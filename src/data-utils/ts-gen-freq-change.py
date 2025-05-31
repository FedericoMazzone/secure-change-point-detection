# This script creates N stationary AR(1) time series, each containing T data points.
# A structural change occurs at time 'tau', where the AR(1) coefficient changes:
# Xt = phi * X_{t-1} + eps_t,   for t = 1, ..., tau
# Xt = phi2 * X_{t-1} + eps_t,  for t = tau+1, ..., T
#
# Key variables:
# - T: Number of points in each time series.
# - eps_t: Innovations (errors) whose distribution can change (e.g., Gaussian, Laplace, etc.).
# - tau: Location of the Change Point, expressed as a fraction in [0,1].
# - phi, phi2: AR(1) coefficients before and after the change point.
# - innovations: This allows the user to modify the distribution of the innovations.
# 
# Outputs:
# - 'data': Self-normalized CUSUM values for each time series.
# - 'v': Frequency (in percentage) of rejecting the null hypothesis ("no change")
#        for different change point locations and heights.
# 
# The results provide insight into the power of detecting changes in AR(1) processes.

import numpy as np
import matplotlib.pyplot as plt
import os

# --- Parameters ---
folder = "data/generated"   # Folder to save the generated time series
distribution = "gaussian"   # Distribution of innovations (can be changed)
n = 20000                   # Number of points in each time series
x = np.zeros(n)             # AR(1) data container
h = 0.4                     # Height of the change in AR(1) coefficient
phi = 0.3                   # AR(1) coefficient before the change
phi2 = phi + h              # AR(1) coefficient after the change
tau = 0.5                   # Change point location in percentage
T1 = int(tau * n)           # Change point location
sd1 = np.sqrt(1 - phi**2)   # Standard deviation for innovations before change
sd2 = np.sqrt(1 - phi2**2)  # Standard deviation for innovations after change
sigma = lambda x: 1 / (1 + np.exp(-x))  # sigmoid function


# --- AR(1) process generation ---
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
    x[t] = phi * x[t-1] + innovation
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
# x = sigma(x)  # Optional: uncomment if needed

# --- Plotting ---
plt.figure(figsize=(15, 6))
plt.plot(x, color='black')
plt.xlabel("time $t$")
plt.ylabel("$X_t$")
plt.title("Time series")
plt.axvline(T1, color='red', linestyle='--', label='Change Point')
plt.legend()
plt.grid()
plt.show()

# --- Save time series and label (T1) to file ---
if not os.path.exists(folder):
    os.makedirs(folder)
filename = f"time_series_{distribution}_n{n}.csv"
filepath = os.path.join(folder, filename)
# save time series as x1,...,xn,T1 with a 3 decimal precision
x = np.round(x, 3)  # Round to 3 decimal places
with open(filepath, "w", newline='') as f:
    f.write(",".join(map(str, x)) + f",{T1}\n")
print(f"Time series saved to {filepath}")

# save the plot
plt.savefig(filepath.replace('.csv', '.png'))
