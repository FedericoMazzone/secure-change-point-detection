import math
import matplotlib.pyplot as plt
import numpy as np

def gen_time_series(
    n,
    distribution="gaussian",
    phi1 = 0.3,
    phi2 = 0.7,
    tau = 0.5,
    burn = 500,
):
    x = np.zeros(n)
    T1 = int(tau * n)
    sd1 = np.sqrt(1 - phi1**2)
    sd2 = np.sqrt(1 - phi2**2)

    for t in range(burn):
        if distribution == "gaussian":
            innovation = sd1 * np.random.randn()
        elif distribution == "laplace":
            u = np.random.rand() - 0.5
            innovation = -4 * np.sign(u) * np.log(1 - 2 * np.abs(u))  # Laplace(0, 4)
        elif distribution == "t":
            innovation = np.random.standard_t(5)
        else:
            raise ValueError("Unsupported distribution type. Use 'gaussian', 'laplace', or 't'.")
        x[0] = phi1 * x[0] + innovation  # Burn-in period to stabilize the process
    for t in range(1, T1):  # Generate AR(1) process before the change point
        if distribution == "gaussian":
            innovation = sd1 * np.random.randn()
        elif distribution == "laplace":
            u = np.random.rand() - 0.5
            innovation = -4 * np.sign(u) * np.log(1 - 2 * np.abs(u))  # Laplace(0, 4)
        elif distribution == "t":
            innovation = np.random.standard_t(5)
        else:
            raise ValueError("Unsupported distribution type. Use 'gaussian', 'laplace', or 't'.")
        x[t] = phi1 * x[t-1] + innovation
    for t in range(T1, n):  # Generate AR(1) process after the change point
        if distribution == "gaussian":
            innovation = sd2 * np.random.randn()
        elif distribution == "laplace":
            u = np.random.rand() - 0.5
            innovation = -4 * np.sign(u) * np.log(1 - 2 * np.abs(u))
        elif distribution == "t":
            innovation = np.random.standard_t(5)
        else:
            raise ValueError("Unsupported distribution type. Use 'gaussian', 'laplace', or 't'.")
        x[t] = phi2 * x[t-1] + innovation

    return x, T1


def turning_rates(x, m=None):
    """
    Split the time series x into blocks of size m and compute the turning rates
    statistic of each block.

    Parameters:
    x (list or array-like): Time series data.
    m (int): Block length.

    Returns:
    np.array: Turning rates for each block.
    """
    n = len(x)

    if m is None:
        m = int(n ** 0.5)

    turn_rates = []

    for t in range(n // m):
        block = x[t * m : (t + 1) * m]

        non_monotonic_count = sum(
            not (block[j] < block[j + 1] < block[j + 2] or
                block[j] > block[j + 1] > block[j + 2])
            for j in range(len(block) - 2)
        )
    
        turn_rates.append(non_monotonic_count / (len(block) - 2) if len(block) > 2 else 1)

    return np.array(turn_rates)

def clip(x, M):
    # clip x to [-M, M]
    return max(-M, min(M, x))


def run_experiment(
    n = 10000,              # length of time series
    num_experiments = 100,  # number of experiments
    M = 1.0,                # clipping bound
):
    t = n // 2       # change point
    m = int(n**0.5)  # block size

    print(f"Length of time series n = {n}")
    print(f"Number of experiments = {num_experiments}")
    print(f"Clipping bound M = {M}")
    print(f"Change point t = {t}")
    print(f"Using block size m = {m}")

    print("\n=== Generating time series data ===")
    X = []
    for i in range(num_experiments):
        print(f"Generating time series {i+1}/{num_experiments}")
        x_exp, _ = gen_time_series(n)
        x_exp = [clip(val, M) for val in x_exp]
        X.append(x_exp)

    print("\n=== Running local DP change point detection experiments ===")
    rel_err = {}
    eps_list = [1] # list(np.arange(0.5, 30.5, 0.5))

    for eps in eps_list:
        errs = []
        for exp_id in range(num_experiments):

            x = X[exp_id]

            # DP parameters
            delta = 1 / len(x)**2
            L = - math.log(delta)
            Delta = 2 * M
            sigma_DP = Delta * (math.sqrt(L) + math.sqrt(L + eps)) / (math.sqrt(2) * eps)
            
            # Client privatization
            noise = np.random.normal(0, sigma_DP, size=len(x))
            x_priv = np.array(x) # + noise

            # Server detection
            tr = turning_rates(x_priv, m)
            tr_s = sum(tr) / len(tr)
            cs = np.abs(np.cumsum(tr - tr_s))
            idx_epoch = np.argmax(cs)
            idx_time = min(len(x), max(1, idx_epoch * m))

            print(idx_time)

            errs.append(abs(idx_time - t) / t)

        rel_err[eps] = float(np.mean(errs))
        print(f"Relative error for eps = {eps}: {rel_err[eps]}")



if __name__ == "__main__":

    

    # plot relative error vs epsilon
    plt.figure(figsize=(10, 6))
    plt.plot(list(rel_err.keys()), list(rel_err.values()), marker='o')
    plt.xlabel("Privacy Budget ε")
    plt.ylabel("Relative Error in Change Point Detection")
    plt.title("Local DP Change Point Detection Relative Error vs Privacy Budget")
    plt.grid()
    plt.savefig("benchmark-localDP.png")
