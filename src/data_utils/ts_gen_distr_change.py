import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t, chi2
import os

# --- Parameters ---
folder = "data/generated"   # Folder to save the generated time series
n = 10000            # Total number of points in the time series
mu_1 = 0            # Mean before change point
mu_2 = 0            # Mean after change point
sigma_1 = 1         # Standard deviation before change point
sigma_2 = 1         # Standard deviation after change point
dist1 = 'normal'    # Distribution before change point
dist2 = 't'   # Distribution after change point
# Available distributions: 'normal', 'laplace', 't', 'uniform', 'chi2'
tau = 0.5           # Change point location in percentage
T1 = int(tau * n)   # Change point location

# Function to sample from standard Laplace(0,1)
def randnlap(size):
    u = np.random.rand(*size) - 0.5
    return -np.sign(u) * np.log(1 - 2 * np.abs(u))

# Function to sample from distribution
def sample_dist(dist, mu, sigma, size):
    if dist == 'normal':
        return mu + sigma * np.random.randn(*size)
    
    elif dist == 'laplace':
        return mu + sigma / np.sqrt(2) * randnlap(size)
    
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



# --- Simulate time series ---
x1 = sample_dist(dist1, mu_1, sigma_1, (1, T1))
x2 = sample_dist(dist2, mu_2, sigma_2, (1, n - T1))
x = np.concatenate([x1, x2], axis=1).flatten()

# --- Plotting ---
plt.figure(figsize=(15, 6))
plt.plot(x, color='black')
plt.xlabel("time $t$")
plt.ylabel("$X_t$")
plt.title(f"Change from {dist1} to {dist2}")
plt.axvline(T1, color='red', linestyle='--', label='Change Point')
plt.legend()
plt.grid()
plt.show()

# --- Save time series and label (T1) to file ---
if not os.path.exists(folder):
    os.makedirs(folder)
filename = f"time_series_n{n}_{dist1}_{mu_1}_{sigma_1}_to_{dist2}_{mu_2}_{sigma_2}.csv"
filepath = os.path.join(folder, filename)
# save time series as x1,...,xn,T1 with a 3 decimal precision
x = np.round(x, 3)  # Round to 3 decimal places
with open(filepath, "w", newline='') as f:
    f.write(",".join(map(str, x)) + f",{T1}\n")
print(f"Time series saved to {filepath}")

# save the plot
plt.savefig(filepath.replace('.csv', '.png'))
