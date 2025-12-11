import numpy as np
import matplotlib.pyplot as plt
import os


def compute_innovations(
    distribution,
    size,
    sd
):
    """
    Generate innovations from specified distribution.

    Parameters:
        distribution (str): Type of distribution ("gaussian", "laplace", "t").
        size (int): Number of samples to generate.
        sd (float): Standard deviation for the innovations.

    Returns:
        np.array: Generated innovations.
    
    Raises:
        ValueError: If an unsupported distribution type is provided.
    
    Example:
        innovations = compute_innovations("gaussian", 1000, 1.0)
    """

    if distribution == "gaussian":
        return sd * np.random.randn(size)
    elif distribution == "laplace":
        u = np.random.rand(size) - 0.5
        return -4 * np.sign(u) * np.log(1 - 2 * np.abs(u))  # Laplace(0, 4)
    elif distribution == "t":
        return np.random.standard_t(5, size=size)
    else:
        raise ValueError("Unsupported distribution type. Use 'gaussian', 'laplace', or 't'.")



def gen_time_series(
    n,
    distribution="gaussian",
    phi1 = 0.3,
    phi2 = 0.7,
    tau = 0.5,
    burn = 500
):
    """
    Generate a time series with a change in AR(1) coefficient at a specified point.

    Parameters:
        n (int): Length of the time series.
        distribution (str): Distribution of innovations ("gaussian", "laplace", "t").
        phi1 (float): AR(1) coefficient before the change.
        phi2 (float): AR(1) coefficient after the change.
        tau (float): Change point location ratio (between 0 and 1).
        burn (int): Burn-in period (initial samples to discard).

    Returns:
        np.array: Generated time series.
        int: Change point index.
    """

    t = int(tau * n)            # change point index
    sd1 = np.sqrt(1 - phi1**2)  # standard deviation for innovations before change
    sd2 = np.sqrt(1 - phi2**2)  # standard deviation for innovations after change

    x = np.zeros(n)

    # Burn-in period (to stabilize the process)
    innovation = compute_innovations(distribution, size=burn, sd=sd1)
    for i in range(burn):
        x[0] = phi1 * x[0] + innovation[i]
    
    # Generate the AR(1) process before the change point
    innovation = compute_innovations(distribution, size=t, sd=sd1)
    for i in range(1, t):
        x[i] = phi1 * x[i-1] + innovation[i]
    
    # Generate the AR(1) process after the change point
    innovation = compute_innovations(distribution, size=n - t, sd=sd2)
    for i in range(t, n):
        x[i] = phi2 * x[i-1] + innovation[i - t]

    return x, t



if __name__ == "__main__":

    folder = "data/generated/frequency-change"   # Folder to save the time series
    plot = True                 # Whether to plot the time series
    distribution = "gaussian"   # Distribution of innovations ["gaussian", "laplace", "t"]
    n = 10000                   # Number of points in each time series
    phi1 = 0.3                  # AR(1) coefficient before the change
    phi2 = 0.7                  # AR(1) coefficient after the change
    tau = 0.5                   # Change point location in percentage
    burn = 500                  # Burn-in period

    # --- Generate time series ---
    x, t = gen_time_series(
        n=n,
        distribution=distribution,
        phi1=phi1,
        phi2=phi2,
        tau=tau,
        burn=burn
    )

    # --- Save time series (x) and label (t) to file ---
    if not os.path.exists(folder):
        os.makedirs(folder)
    filename = f"time_series_{distribution}_n{n}.csv"
    filepath = os.path.join(folder, filename)
    # save time series as x1,...,xn,t with a 3 decimal precision
    x = np.round(x, 3)  # Round to 3 decimal places
    with open(filepath, "w", newline='') as f:
        f.write(",".join(map(str, x)) + f",{t}\n")
    print(f"Time series saved to {filepath}")

    # --- Plotting ---
    if plot:
        plt.figure(figsize=(15, 6))
        plt.plot(x, color='black')
        plt.xlabel("time $t$")
        plt.ylabel("$X_t$")
        plt.title("Time series")
        plt.axvline(t, color='red', linestyle='--', label='Change Point')
        plt.legend()
        plt.grid()
        plt.savefig(filepath.replace('.csv', '.png'))

