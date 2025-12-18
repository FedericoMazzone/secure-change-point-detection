"""
Change Point Detection Plotter

This module provides functionality to perform change point detection on time
series data using CUSUM on block statistics, and visualize the results.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys


def cpd(
    x,
    change_type="frequency",
    m=None
):
    """
    Change point detection using CUSUM on block statistics.

    Parameters:
        x (list or np.array): Input time series data.
        change_type (str): Type of change to detect ("mean", "variance", "frequency").
        m (int): Block size for computing statistics. If None, defaults to sqrt(n).
    
    Returns:
        int: Estimated change point index.
        np.array: Block statistic values.
        np.array: CUSUM statistic values.
    
    Raises:
        ValueError: If an unsupported change_type is provided.
    
    Example:
        idx_time = cpd(x, change_type="frequency", m=50)
    """
    n = len(x)
    if m is None:
        m = int(n ** 0.5)
    
    q = []
    for t in range(n // m):
        block = x[t * m : (t + 1) * m]
        if change_type == "mean":
            q.append(np.mean(block))
        elif change_type == "variance":
            q.append(np.var(block))
        elif change_type == "frequency":
            non_monotonic_count = sum(
                not (block[j] < block[j + 1] < block[j + 2] or block[j] > block[j + 1] > block[j + 2])
                for j in range(len(block) - 2)
            )
            q.append(non_monotonic_count / (len(block) - 2) if len(block) > 2 else 1)
        else:
            raise ValueError("Unsupported change_type. Use 'mean', 'variance', or 'frequency'.")
    q = np.array(q)
    
    q_t = np.mean(q)
    cs = np.abs(np.cumsum(q - q_t))
    idx_epoch = np.argmax(cs) + 1
    idx_time = min(n, max(1, idx_epoch * m))

    return idx_time, q, cs



if __name__ == "__main__":

    # --- Argument parser ---
    parser = argparse.ArgumentParser(
        description="Run change-point detection on a time series from CSV"
    )

    parser.add_argument(
        "filepath",
        type=str,
        help="Path to the CSV file containing the time series data"
    )
    parser.add_argument(
        "--change_type",
        type=str,
        choices=["mean", "variance", "frequency"],
        default="frequency",
        help="Type of change to detect"
    )
    parser.add_argument(
        "--m",
        type=int,
        default=None,
        help="Block size for statistics (default: sqrt(n))"
    )

    args = parser.parse_args()

    filepath = args.filepath
    change_type = args.change_type
    m = args.m

    # --- Read the time series data ---
    with open(filepath, 'r') as f:
        data = f.read().strip().split(',')
    x = list(map(float, data[:-1]))     # read the time series values
    t = int(float(data[-1]))           # read the change point index
    if m is None:
        m = int(len(x) ** 0.5)
    print(f"Loaded time series from {filepath}")
    print(f"Size of time series: {len(x)}")
    print(f"Ground truth change point at t={t}")
    print(f"Using block size m={m}")

    # --- Change point detection ---
    computed_t, q, q_cumsum = cpd(x, change_type=change_type, m=m)
    print(f"Computed change point at t={computed_t}")

    # --- Plotting ---
    fig, axs = plt.subplots(3, 1, figsize=(15, 12))
    # Plot the original time series
    axs[0].plot(x, color='black')
    axs[0].set_ylabel("$X_t$")
    axs[0].set_title(f"Time series from {filepath}")
    axs[0].axvline(t, color='red', linestyle='--', label='Change Point')
    axs[0].axvline(computed_t, color='green', linestyle='--', label='Computed Change Point')
    axs[0].legend()
    axs[0].grid()
    # Plot the q values
    axs[1].plot(q, color='blue')
    axs[1].set_ylabel("$q_t$")
    axs[1].set_title("Summarized block statistics")
    axs[1].axvline(t // m - 1, color='red', linestyle='--', label='Change Point')
    axs[1].axvline(computed_t // m - 1, color='green', linestyle='--', label='Computed Change Point')
    axs[1].legend()
    axs[1].grid()
    # Plot the CUSUM statistic
    axs[2].plot(q_cumsum, color='green')
    axs[2].set_ylabel("CUSUM")
    axs[2].set_title("CUSUM statistic")
    axs[2].axvline(t // m - 1, color='red', linestyle='--', label='Change Point')
    axs[2].axvline(computed_t // m - 1, color='green', linestyle='--', label='Computed Change Point')
    axs[2].legend()
    axs[2].grid()
    # Save the plot to a file (change extension to .png)
    imgpath = filepath.replace('.csv', '.png')
    plt.savefig(imgpath)
    print(f"Plot saved to {imgpath}")
