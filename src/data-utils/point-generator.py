import sys
from pathlib import Path

# import matplotlib.pyplot as plt
import numpy as np

# from mpl_toolkits.mplot3d import Axes3D

# Set random seed for reproducibility
np.random.seed(42)

# Define the total number of points and the number of clusters
num_points = int(sys.argv[1])
num_clusters = int(sys.argv[2])
num_dims = int(sys.argv[3]) if len(sys.argv) > 3 else 2
print(f"Generating {num_points} {num_dims}d points around {num_clusters} clusters.")

# Define the actual cluster centers
centers_dict = {
    (2, 2) : np.array([[0.25, 0.25], [0.75, 0.75]]),
    (2, 3) : np.array([[0.2, 0.8], [0.5, 0.5], [0.8, 0.2]]),
    (2, 4) : np.array([[0.2, 0.2], [0.2, 0.8], [0.8, 0.2], [0.8, 0.8]]),
    (2, 5) : np.array([[0.2, 0.2], [0.2, 0.8], [0.5, 0.5], [0.8, 0.2], [0.8, 0.8]]),
    (2, 8) : np.array([[0.4, 0.2], [0.4, 0.6], [0.6, 0.4], [0.6, 0.8], [0.2, 0.8], [0.8, 0.6], [0.2, 0.4], [0.8, 0.2]]),
    (3, 3) : np.array([[0.2, 0.2, 0.2], [0.8, 0.8, 0.2], [0.5,0.5,0.8]]),
    (8, 8) : np.array([
        [0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9],
        [0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9],
        [0.1, 0.9, 0.5, 0.5, 0.5, 0.5, 0.9, 0.1],
        [0.5, 0.5, 0.1, 0.9, 0.9, 0.1, 0.1, 0.5]
        [0.5, 0.5, 0.5, 0.5, 0.9, 0.1, 0.9, 0.1],
        [0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.1, 0.9],
        [0.9, 0.1, 0.9, 0.1, 0.1, 0.1, 0.9, 0.9],
        [0.9, 0.9, 0.1, 0.1, 0.9, 0.1, 0.1, 0.9],
    ])
}
cluster_centers = centers_dict[num_dims, num_clusters]

# Calculate points per cluster by dividing evenly
num_points_per_cluster = [num_points // num_clusters] * num_clusters
num_points_per_cluster[0] += num_points % num_clusters  # Add remainder to the first cluster

# Generate points around each center with a small random noise
points = []
for i, center in enumerate(cluster_centers):
    # Add noise and clip points to stay within [0,1] for both x and y coordinates
    points_cluster = center + 0.03 * np.random.randn(num_points_per_cluster[i], num_dims)
    points_cluster = np.clip(points_cluster, 0, 1)
    points_cluster = np.hstack((points_cluster, np.full((points_cluster.shape[0], 1), i)))
    points.append(points_cluster)

# Combine all points into one array
points = np.vstack(points)

# Create a 'data' directory if it doesn't exist
data_folder = Path("data")
data_folder.mkdir(exist_ok=True)

# Save the points in a file with two columns (x and y)
filepath = data_folder / f"points_d{num_dims}_n{num_points}_k{num_clusters}.csv"
with open(filepath, "w") as f:
    f.writelines(",".join(f"{x:.3f}" for x in point[:-1]) + f",{int(point[-1])}\n" for point in points)

# # Plot the points and centers for visualization
# if num_dims == 2:
#     plt.figure(figsize=(6, 6))
#     plt.scatter(points[:, 0], points[:, 1], c='blue', label='Points')
#     plt.scatter(cluster_centers[:, 0], cluster_centers[:, 1], c='red', marker='x', s=100, label='Centers (Ground Truth)')
#     plt.title('Generated Points and Cluster Centers')
#     plt.xlabel('X')
#     plt.ylabel('Y')
#     plt.legend()
#     plt.grid(True)
#     plt.xlim(0, 1)
#     plt.ylim(0, 1)
#     plt.show()
# elif num_dims == 3:
#     fig = plt.figure(figsize=(8, 8))
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='blue', label='Points')
#     ax.scatter(cluster_centers[:, 0], cluster_centers[:, 1], cluster_centers[:, 2], c='red', marker='x', s=100, label='Centers (Ground Truth)')
#     ax.set_title('Generated Points and Cluster Centers (3D)')
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.legend()
#     ax.grid(True)
#     ax.set_xlim(0, 1)
#     ax.set_ylim(0, 1)
#     ax.set_zlim(0, 1)
#     plt.show()