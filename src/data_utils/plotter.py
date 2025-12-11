import matplotlib.pyplot as plt
import pandas as pd
import sys


filepath = sys.argv[1]  # Path to the CSV file containing the time series data
# read the time series data from the CSV file
# it is in the form x1,...,xn,T1 where T1 is the change point
x = pd.read_csv(filepath, header=None).values.flatten()[:-1]  # Read the time series data
T1 = pd.read_csv(filepath, header=None).values.flatten()[-1]  # Read the change point

# turn from frequency-change to mean-change problem
# block size is sqrt(n)
n = len(x)
block_size = int(n**0.5)
# for each block, take a sliding window of size 3 and count for each block how many triplets are NOT monotonic
q = []
for i in range(0, n - block_size + 1, block_size):
    block = x[i:i + block_size]
    non_monotonic_count = sum(
        not (block[j] < block[j + 1] < block[j + 2] or
             block[j] > block[j + 1] > block[j + 2])
        for j in range(len(block) - 2)
    )
    q.append(non_monotonic_count)
    # normalize by the number of triplets in the block
    q[-1] /= (block_size - 2) if block_size > 2 else 1
# compute the CUSUM statistic of the q values (no Pandas)
q_cumsum = [0] * len(q)
for i in range(len(q)):
    if i == 0:
        q_cumsum[i] = q[i]
    else:
        q_cumsum[i] = q_cumsum[i - 1] + (q[i] - sum(q) / len(q))
for i in range(len(q)):
    q_cumsum[i] = abs(q_cumsum[i])

# --- Plotting ---
# plot in multiple subplots: (1) the original time series, (2) the q values, (3) the CUSUM statistic
# plt.style.use('seaborn-darkgrid')
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))
# # Plot the original time series
# ax1.plot(x, color='black')
# ax1.set_xlabel("time $t$")
# ax1.set_ylabel("$X_t$")
# ax1.set_title(f"Time series from {filepath}")
# ax1.axvline(T1, color='red', linestyle='--', label='Change Point')
# ax1.legend()
# ax1.grid()
# # Plot the CUSUM statistic
# ax2.plot(q_cumsum, color='blue')
# ax2.set_xlabel("Block Index")
# ax2.set_ylabel("CUSUM Statistic")
# ax2.set_title("CUSUM Statistic of Non-Monotonic Triplets")
# ax2.axvline(T1 // block_size, color='red', linestyle='--', label='Change Point')
# ax2.legend()
# ax2.grid()


fig, axs = plt.subplots(3, 1, figsize=(15, 12))
# Plot the original time series
axs[0].plot(x, color='black')
axs[0].set_ylabel("$X_t$")
axs[0].set_title(f"Time series from {filepath}")
axs[0].axvline(T1, color='red', linestyle='--', label='Change Point')
axs[0].legend()
axs[0].grid()
# Plot the q values
axs[1].plot(q, color='blue')
axs[1].set_ylabel("$q_t$")
axs[1].set_title("Non-monotonic triplet count per block")
axs[1].axvline(T1 // block_size, color='red', linestyle='--', label='Change Point')
axs[1].legend()
axs[1].grid()
# Plot the CUSUM statistic
axs[2].plot(q_cumsum, color='green')
axs[2].set_ylabel("CUSUM")
axs[2].set_title("CUSUM statistic of non-monotonic triplet counts")
axs[2].axvline(T1 // block_size, color='red', linestyle='--', label='Change Point')
axs[2].legend()
axs[2].grid()


# plt.figure(figsize=(15, 6))
# plt.plot(x, color='black')
# plt.xlabel("time $t$")
# plt.ylabel("$X_t$")
# plt.title(f"Time series from {filepath}")
# plt.axvline(T1, color='red', linestyle='--', label='Change Point')
# plt.legend()
# plt.grid()
# plt.show()

# Save the plot to a file (change extension to .png)
plt.savefig(filepath.replace('.csv', '.png'))
