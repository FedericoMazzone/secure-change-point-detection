"""
Benchmarking local differential privacy (DP) change point detection.
Generates synthetic time series data with a change point in frequency, applies
local DP noise, and evaluates change point detection accuracy.
Measures accuracy using relative error between detected and true change points.
Saves results to CSV and generates plots of relative error vs. privacy budget
(epsilon).
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import os
from multiprocessing import Process, Queue, cpu_count

from data_utils.ts_gen_freq_change import gen_time_series



def turning_rates(x, m=None):
    """
    Compute turning rates of time series x using block size m.

    Parameters:
        x (np.array): Input time series data.
        m (int, optional): Block size for computing turning rates.
                           Defaults to sqrt(len(x)).

    Returns:
        np.array: Array of turning rates for each block.
    
    Example:
        x = np.array([1, 3, 2, 4, 5, 3, 6, 7])
        tr = turning_rates(x, m=4)
        print(tr)  # Output: array of turning rates for each block
    """
    n = len(x)
    if m is None:
        m = int(n ** 0.5)
    turn_rates = []
    for t in range(n // m):
        block = x[t * m : (t + 1) * m]
        non_monotonic_count = sum(
            not (block[j] < block[j + 1] < block[j + 2] or block[j] > block[j + 1] > block[j + 2])
            for j in range(len(block) - 2)
        )
        turn_rates.append(non_monotonic_count / (len(block) - 2) if len(block) > 2 else 1)
    return np.array(turn_rates)



def compute_error_for_eps(eps, X, t, m, M, q):
    """
    Compute relative error for a given epsilon over multiple experiments.

    Parameters:
        eps (float): Privacy budget.
        X (list of np.array): List of time series data.
        t (int): True change point time index.
        m (int): Block size for turning rates.
        M (float): Clipping bound for time series values.
        q (multiprocessing.Queue): Queue to send results back to main process.
    
    Returns:
        None: Results are sent back via the queue.
    """
    errs = []
    for x in X:
        # DP parameters
        delta = 1 / len(x) ** 2
        L = -math.log(delta)
        Delta = 2 * M
        sigma_DP = Delta * (math.sqrt(L) + math.sqrt(L + eps)) / (math.sqrt(2) * eps)

        # Client privatization
        noise = np.random.normal(0, sigma_DP, size=len(x)) if not math.isnan(sigma_DP) else np.zeros(len(x))
        x_priv = np.array(x) + noise

        # Server detection
        tr = turning_rates(x_priv, m)
        tr_s = np.mean(tr)
        cs = np.abs(np.cumsum(tr - tr_s))
        idx_epoch = np.argmax(cs) + 1
        idx_time = min(len(x), max(1, idx_epoch * m))

        errs.append(abs(idx_time - t) / t)

    mean_err = float(np.mean(errs))
    std_dev_err = float(np.std(errs))
    q.put((eps, mean_err, std_dev_err))  # send result back to main process



def run_experiment(
    n=10000,
    num_experiments=100,
    eps_list=[1.0],
    M=1.0
):
    """
    Run local DP change point detection experiments.

    Parameters:
        n (int): Length of each time series.
        num_experiments (int): Number of experiments to run.
        eps_list (list of float): List of epsilon values to test.
        M (float): Clipping bound for time series values.
    
    Returns:
        dict: Mapping from epsilon to (mean relative error, std dev of error).
    """

    m = int(n**0.5)     # Block size for turning rates

    print(f"Generating {num_experiments} time series of length {n}...")
    X = []
    for i in range(num_experiments):
        x_exp, t = gen_time_series(n)
        X.append(np.clip(x_exp, -M, M))  # Clip to [-M, M]

    print("Running experiments in parallel for different epsilons...")
    rel_err = {}
    q = Queue()
    processes = []

    n_workers = min(cpu_count(), len(eps_list))

    # Start a process for each epsilon
    for eps in eps_list:
        p = Process(target=compute_error_for_eps, args=(eps, X, t, m, M, q))
        processes.append(p)
        p.start()

        # If we have spawned max workers, wait for them to finish before starting more
        if len(processes) >= n_workers:
            for p_wait in processes:
                p_wait.join()
            processes = []

    # join any remaining processes
    for p in processes:
        p.join()

    # Collect results from queue
    while not q.empty():
        eps, mean_err, std_dev_err = q.get()
        rel_err[eps] = (mean_err, std_dev_err)

    # Sort rel_err by epsilon
    rel_err = dict(sorted(rel_err.items()))
    return rel_err



if __name__ == "__main__":

    # Experiment configuration
    experiment_config = [                           # (n, M, num_experiments)
        (1000, 1.0, 100000),
        (10000, 1.0, 10000),
        (100000, 1.0, 1000),
        (1000000, 1.0, 100),
    ]
    eps_list = np.arange(0.5, 30.5, 0.5).tolist()   # Epsilon values to test

    # Run experiments
    rel_err_all = {}
    for n_val, M, num_exp in experiment_config:
        rel_err_all[(n_val, M, num_exp)] = run_experiment(n=n_val, num_experiments=num_exp, eps_list=eps_list, M=M)
        print(f"Results for n={n_val}, M={M}, num_experiments={num_exp}:")
        for eps in rel_err_all[(n_val, M, num_exp)]:
            mean_err, std_dev_err = rel_err_all[(n_val, M, num_exp)][eps]
            print(f"  Epsilon: {eps:.2f}, Relative Error: {mean_err:.4f}, Std Dev: {std_dev_err:.4f}")
    
    # Run experiments for baseline (no DP)
    print("Running baseline experiments (no DP)...")
    baseline_errors = {}
    M = float("inf")                    # Large M to avoid clipping
    eps_list_baseline = [float("inf")]  # Representing no DP
    for n_val, _, num_exp in experiment_config:
        baseline_errors[(n_val, num_exp)] = run_experiment(n=n_val, num_experiments=num_exp, eps_list=eps_list_baseline, M=M)
        print(f"Baseline Results for n={n_val}, num_experiments={num_exp}:")
        for eps in baseline_errors[(n_val, num_exp)]:
            mean_err, std_dev_err = baseline_errors[(n_val, num_exp)][eps]
            print(f"  Epsilon: {eps}, Relative Error: {mean_err:.4f}, Std Dev: {std_dev_err:.4f}")

    # Save results to csv
    if not os.path.exists("logs/local-dp"):
        os.makedirs("logs/local-dp")
    with open("logs/local-dp/benchmark-localDP.csv", "w") as f:
        f.write("n,M,epsilon,relative_error,std_dev_error\n")
        for n_val, M, num_exp in rel_err_all:
            for eps in rel_err_all[(n_val, M, num_exp)]:
                rel_err, std_dev_err = rel_err_all[(n_val, M, num_exp)][eps]
                f.write(f"{n_val},{M},{eps},{rel_err},{std_dev_err}\n")
        for n_val, num_exp in baseline_errors:
            for eps in baseline_errors[(n_val, num_exp)]:
                rel_err, std_dev_err = baseline_errors[(n_val, num_exp)][eps]
                f.write(f"{n_val},inf,inf,{rel_err},{std_dev_err}\n")
    
    # Plot results
    plt.figure(figsize=(10, 6))
    color_map = {}  # (n_val, num_exp) -> color
    for n_val, M, num_exp in rel_err_all:
        epsilons = list(rel_err_all[(n_val, M, num_exp)].keys())
        errors = [val[0] for val in rel_err_all[(n_val, M, num_exp)].values()]
        std_devs = [val[1] for val in rel_err_all[(n_val, M, num_exp)].values()]
        sems = [std_dev / np.sqrt(num_exp) for std_dev in std_devs]
        line, = plt.plot(
            epsilons,
            errors,
            label=f"n={n_val}, num_exp={num_exp}"
        )
        # store the color used for this curve
        color_map[(n_val, num_exp)] = line.get_color()

        plt.fill_between(
            epsilons,
            np.array(errors) - 2 * np.array(sems),
            np.array(errors) + 2 * np.array(sems),
            alpha=0.2,
            color=line.get_color()
        )
    # baseline lines with matching colors
    for n_val, num_exp in baseline_errors:
        baseline_err, baseline_std = baseline_errors[(n_val, num_exp)][float("inf")]
        plt.hlines(
            baseline_err,
            xmin=0,
            xmax=max(eps_list),
            colors=color_map[(n_val, num_exp)],
            linestyles="dashed",
            label=f"baseline n={n_val}"
        )
    plt.xlabel("Privacy Budget (epsilon)")
    plt.ylabel("Relative Error")
    plt.title("Local DP Change Point Detection Accuracy")
    plt.legend()
    plt.grid(True)
    plt.savefig("logs/local-dp/benchmark-localDP.png")
    print("Plot saved to logs/local-dp/benchmark-localDP.png")
