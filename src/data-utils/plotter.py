import matplotlib.pyplot as plt
import pandas as pd
import sys


filepath = sys.argv[1]  # Path to the CSV file containing the time series data
# read the time series data from the CSV file
# it is in the form x1,...,xn,T1 where T1 is the change point
x = pd.read_csv(filepath, header=None).values.flatten()[:-1]  # Read the time series data
T1 = pd.read_csv(filepath, header=None).values.flatten()[-1]  # Read the change point

# --- Plotting ---
plt.figure(figsize=(15, 6))
plt.plot(x, color='black')
plt.xlabel("time $t$")
plt.ylabel("$X_t$")
plt.title(f"Time series from {filepath}")
plt.axvline(T1, color='red', linestyle='--', label='Change Point')
plt.legend()
plt.grid()
plt.show()

# Save the plot to a file (change extension to .png)
plt.savefig(filepath.replace('.csv', '.png'))
