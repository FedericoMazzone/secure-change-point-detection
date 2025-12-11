import math
import matplotlib.pyplot as plt
import numpy as np

from data_utils.ts_gen_freq_change import gen_time_series


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
    eps_list = [1.0],       # list of privacy budgets
    M = 1.0,                # clipping bound
):
    t = n // 2       # change point
    m = int(n**0.5)  # block size

    print(f"Length of time series n = {n}")
    print(f"Number of experiments = {num_experiments}")
    print(f"Privacy budgets epsilon = {eps_list}")
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

    for eps in eps_list:
        errs = []
        for exp_id in range(num_experiments):

            x = X[exp_id]

            # DP parameters
            delta = 1 / len(x)**2
            L = -math.log(delta)
            Delta = 2 * M
            sigma_DP = Delta * (math.sqrt(L) + math.sqrt(L + eps)) / (math.sqrt(2) * eps)
            
            # Client privatization
            noise = np.random.normal(0, sigma_DP, size=len(x))
            x_priv = np.array(x) + noise

            # Server detection
            tr = turning_rates(x_priv, m)
            tr_s = sum(tr) / len(tr)
            cs = np.abs(np.cumsum(tr - tr_s))
            idx_epoch = np.argmax(cs)
            idx_time = min(len(x), max(1, idx_epoch * m))

            errs.append(abs(idx_time - t) / t)

        rel_err[eps] = float(np.mean(errs))
        print(f"Relative error for eps = {eps}: {rel_err[eps]}")
    
    return rel_err



if __name__ == "__main__":

    # Setup experiment configurations
    experiment_config = [                       # (n, num_experiments)
        (1000, 10),
        (10000, 1),
        (100000, 1),
        (1000000, 1)
    ]
    eps_list = np.arange(0.5, 30.5, 0.5).tolist()   # list of privacy budgets
    M = 1.0                                         # clipping bound

    # Run experiments
    rel_err = {}
    for n, num_experiments in experiment_config:
        rel_err[n] = run_experiment(
            n = n,
            num_experiments = num_experiments,
            eps_list = eps_list,
            M = M
        )

    # Plot results in one figure
    plt.figure(figsize=(10, 6))
    for n in rel_err:
        epsilons = list(rel_err[n].keys())
        errors = list(rel_err[n].values())
        plt.plot(epsilons, errors, marker='o', label=f'n={n}')
    plt.xlabel('Privacy Budget (epsilon)')
    plt.ylabel('Relative Error')
    plt.title('Local DP Change Point Detection Accuracy')
    plt.legend()
    plt.grid(True)
    plt.savefig('benchmark-localDP.png')