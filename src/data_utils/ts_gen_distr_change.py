"""
Generate time series with a change in distribution.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t, chi2
import os


def sample_dist(
    dist,
    mu,
    sigma,
    size
):
    """
    Sample from specified distribution with given mean and standard deviation.

    Parameters:
        dist (str): Type of distribution ("normal", "laplace", "t", "uniform", "chi2").
        mu (float): Mean of the distribution.
        sigma (float): Standard deviation of the distribution.
        size (tuple): Size of the output array.
    
    Returns:
        np.array: Samples from the specified distribution.
    
    Raises:
        ValueError: If an unsupported distribution type is provided.
    
    Example:
        samples = sample_dist("normal", 0, 1, (1000,))
    """

    if dist == 'normal':
        return mu + sigma * np.random.randn(*size)
    
    elif dist == 'laplace':
        u = np.random.rand(size) - 0.5
        v = -np.sign(u) * np.log(1 - 2 * np.abs(u))
        return mu + sigma / np.sqrt(2) * v
    
    elif dist == 't':
        return mu + sigma * t.rvs(df=5, size=size) / np.sqrt(5/3)
    
    elif dist == 'uniform':
        a = mu - np.sqrt(3) * sigma
        b = mu + np.sqrt(3) * sigma
        return np.random.uniform(a, b, size=size)
    
    elif dist == 'chi2':
        k = 5
        raw = chi2.rvs(df=k, size=size)
        standardized = (raw - k) / np.sqrt(2 * k)
        return mu + sigma * standardized
    
    else:
        raise ValueError(f"Unknown distribution: {dist}")



def gen_time_series(
    n,
    dist1="normal",
    dist2="normal",
    mu1=0.0,
    mu2=0.0,
    sigma1=1.0,
    sigma2=1.0,
    tau=0.5
):
    """
    Generate a time series with a change in distribution at a specified point.

    Parameters:
        n (int): Length of the time series.
        dist1 (str): Distribution before the change point.
        dist2 (str): Distribution after the change point.
        mu1 (float): Mean before the change point.
        mu2 (float): Mean after the change point.
        sigma1 (float): Standard deviation before the change point.
        sigma2 (float): Standard deviation after the change point.
        tau (float): Change point location as a fraction of n.

    Returns:
        np.array: Generated time series.
        int: Change point index.
    """

    t = int(tau * n)   # Change point location

    # --- Simulate time series ---
    x1 = sample_dist(dist1, mu1, sigma1, (1, t))
    x2 = sample_dist(dist2, mu2, sigma2, (1, n - t))
    x = np.concatenate([x1, x2], axis=1).flatten()

    return x, t



if __name__ == "__main__":

    # --- Argument parser ---
    parser = argparse.ArgumentParser(
        description="Generate time series with a distribution change"
    )

    parser.add_argument(
        "--folder",
        type=str,
        default="data/generated/distr_change",
        help="Folder to save the generated time series"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Whether to plot the generated time series"
    )
    parser.add_argument(
        "--dist1",
        type=str,
        choices=["normal", "laplace", "t", "uniform", "chi2"],
        default="normal",
        help="Distribution before change point"
    )
    parser.add_argument(
        "--dist2",
        type=str,
        choices=["normal", "laplace", "t", "uniform", "chi2"],
        default="normal",
        help="Distribution after change point"
    )
    parser.add_argument(
        "--n",
        type=int,
        default=10000,
        help="Total number of points in the time series"
    )
    parser.add_argument(
        "--mu1",
        type=float,
        default=0.0,
        help="Mean before change point"
    )
    parser.add_argument(
        "--mu2",
        type=float,
        default=0.0,
        help="Mean after change point"
    )
    parser.add_argument(
        "--sigma1",
        type=float,
        default=1.0,
        help="Standard deviation before change point"
    )
    parser.add_argument(
        "--sigma2",
        type=float,
        default=1.0,
        help="Standard deviation after change point"
    )
    parser.add_argument(
        "--tau",
        type=float,
        default=0.5,
        help="Change point location as a fraction of n"
    )

    args = parser.parse_args()

    folder = args.folder
    plot = args.plot
    dist1 = args.dist1
    dist2 = args.dist2
    n = args.n
    mu1 = args.mu1
    mu2 = args.mu2
    sigma1 = args.sigma1
    sigma2 = args.sigma2
    tau = args.tau

    # --- Generate time series ---
    x, t = gen_time_series(
        n=n,
        dist1=dist1,
        dist2=dist2,
        mu1=mu1,
        mu2=mu2,
        sigma1=sigma1,
        sigma2=sigma2,
        tau=tau
    )

    # --- Save time series (x) and label (t) to file ---
    if not os.path.exists(folder):
        os.makedirs(folder)
    filename = f"time_series_n{n}_{dist1}_{mu1}_{sigma1}_to_{dist2}_{mu2}_{sigma2}.csv"
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
