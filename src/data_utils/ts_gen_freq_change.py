"""
Generate AR(1) time series with a change in frequency (AR coefficient).
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

from scipy.signal import lfilter


def compute_innovations(
    distribution,
    size,
    sd
):
    """
    Generate innovations from specified distribution.

    Parameters:
        distribution (str): Type of distribution ("normal", "laplace", "t").
        size (int): Number of samples to generate.
        sd (float): Standard deviation for the innovations.

    Returns:
        np.array: Generated innovations.
    
    Raises:
        ValueError: If an unsupported distribution type is provided.
    
    Example:
        innovations = compute_innovations("normal", 1000, 1.0)
    """

    if distribution == "normal":
        return sd * np.random.randn(size)
    elif distribution == "laplace":
        u = np.random.rand(size) - 0.5
        return -4 * np.sign(u) * np.log(1 - 2 * np.abs(u))  # Laplace(0, 4)
    elif distribution == "t":
        return np.random.standard_t(5, size=size)
    else:
        raise ValueError("Unsupported distribution type. Use 'normal', 'laplace', or 't'.")



def gen_time_series(
    n,
    distribution="normal",
    phi1 = 0.3,
    phi2 = 0.7,
    tau = 0.5,
    burn = 500
):
    """
    Generate a time series with a change in AR(1) coefficient at a specified point.

    Parameters:
        n (int): Length of the time series.
        distribution (str): Distribution of innovations ("normal", "laplace", "t").
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
    
    # --- Burn-in period (to stabilize the process) ---
    if burn > 0:
        innovation_burn = compute_innovations(distribution, size=burn, sd=sd1)
        # for i in range(burn):
        #     x[0] = phi1 * x[0] + innovation_burn[i]

        x_burn = lfilter([1.0], [1.0, -phi1], innovation_burn)
        x[0] = x_burn[-1]  # The final burn-in value goes to x[0]
    else:
        x[0] = 0.0
    
    # --- Pre-change segment (indices 1 to t-1) ---
    if t > 1:
        innovation_1 = compute_innovations(distribution, size=t-1, sd=sd1)
        # for i in range(1, t):
        #     x[i] = phi1 * x[i-1] + innovation_1[i - 1]

        # lfilter assumes zero initial condition, but we start from x[0]
        x_pre_zero_init = lfilter([1.0], [1.0, -phi1], innovation_1)
        
        # Account for starting from x[0]: add x[0] * phi1^k for k=1,2,...,t-1
        k = np.arange(1, t)
        carry = x[0] * (phi1 ** k)
        x[1:t] = x_pre_zero_init + carry
    elif t == 1:
        # No pre-change segment (change happens immediately after burn-in)
        pass

    # ---------- Post-change segment (indices t to n-1) ----------
    m = n - t
    if m > 0:
        innovation_2 = compute_innovations(distribution, size=n-t, sd=sd2)
        # for i in range(t, n):
        #     x[i] = phi2 * x[i-1] + innovation_2[i - t]

        x_post_zero_init = lfilter([1.0], [1.0, -phi2], innovation_2)
    
    # Account for starting from x[t-1]
    k = np.arange(1, m + 1)
    carry = x[t-1] * (phi2 ** k)
    x[t:n] = x_post_zero_init + carry

    return x, t



if __name__ == "__main__":

    # --- Argument parser ---
    parser = argparse.ArgumentParser(
        description="Generate AR(1) time series with a frequency change"
    )

    parser.add_argument(
        "--folder",
        type=str,
        default="data/generated/frequency-change",
        help="Folder to save the time series"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Whether to plot the time series"
    )
    parser.add_argument(
        "--distribution",
        type=str,
        choices=["normal", "laplace", "t"],
        default="normal",
        help="Distribution of innovations"
    )
    parser.add_argument(
        "--n",
        type=int,
        default=10000,
        help="Number of points in each time series"
    )
    parser.add_argument(
        "--phi1",
        type=float,
        default=0.3,
        help="AR(1) coefficient before the change"
    )
    parser.add_argument(
        "--phi2",
        type=float,
        default=0.7,
        help="AR(1) coefficient after the change"
    )
    parser.add_argument(
        "--tau",
        type=float,
        default=0.5,
        help="Change point location as a fraction of n"
    )
    parser.add_argument(
        "--burn",
        type=int,
        default=500,
        help="Burn-in period"
    )

    args = parser.parse_args()

    folder = args.folder
    plot = args.plot
    distribution = args.distribution
    n = args.n
    phi1 = args.phi1
    phi2 = args.phi2
    tau = args.tau
    burn = args.burn

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

