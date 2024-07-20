import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def generate_random_points(n, seed=None):
    np.random.seed(seed)
    return np.random.rand(n, 2)

def generate_points_boundary(n, margin=0.1, seed=None):
    np.random.seed(seed)
    points = np.concatenate((
        np.random.rand(n // 4, 2) * [1, margin],  # bottom
        np.random.rand(n // 4, 2) * [1, margin] + [0, 1 - margin],  # top
        np.random.rand(n // 4, 2) * [margin, 1],  # left
        np.random.rand(n // 4, 2) * [margin, 1] + [1 - margin, 0]  # right
    ))
    return points

def generate_points_ring(n, radius=0.5, width=0.1, seed=None):
    np.random.seed(seed)
    angles = np.random.rand(n) * 2 * np.pi
    radii = radius + (np.random.rand(n) - 0.5) * width
    points = np.column_stack((radii * np.cos(angles) + 0.5, radii * np.sin(angles) + 0.5))
    return points

def generate_points_clusters(n, num_clusters=4, cluster_std=0.05, seed=None):
    np.random.seed(seed)
    centers = np.random.rand(num_clusters, 2)
    points = np.vstack([np.random.randn(n // num_clusters, 2) * cluster_std + center for center in centers])
    return points

def construct_gabriel_graph(points):
    G = nx.Graph()
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            u, v = points[i], points[j]
            midpoint = (u + v) / 2
            radius = np.linalg.norm(u - v) / 2
            circle_contains_no_other_points = all(np.linalg.norm(points[k] - midpoint) > radius for k in range(n) if k != i and k != j)
            if circle_contains_no_other_points:
                G.add_edge(i, j, weight=1)  # All edges have weight 1 for unweighted graph
    return G

def average_shortest_path_length(G):
    lengths = dict(nx.all_pairs_dijkstra_path_length(G))
    total_length = 0
    count = 0
    for u in lengths:
        for v in lengths[u]:
            if u != v:
                total_length += lengths[u][v]
                count += 1
    return total_length / count

def longest_shortest_path_length(G):
    lengths = dict(nx.all_pairs_dijkstra_path_length(G))
    max_length = 0
    for u in lengths:
        for v in lengths[u]:
            if u != v and lengths[u][v] > max_length:
                max_length = lengths[u][v]
    return max_length

# Parameters
n = 200  # number of points
seed = 42  # random seed for reproducibility

# Generate points using different methods
points_random = generate_random_points(n, seed=seed)
points_boundary = generate_points_boundary(n, seed=seed)
points_ring = generate_points_ring(n, seed=seed)
points_clusters = generate_points_clusters(n, seed=seed)

# Construct Gabriel graphs
G_random = construct_gabriel_graph(points_random)
G_boundary = construct_gabriel_graph(points_boundary)
G_ring = construct_gabriel_graph(points_ring)
G_clusters = construct_gabriel_graph(points_clusters)

# Calculate average shortest path length
avg_length_random = average_shortest_path_length(G_random)
avg_length_boundary = average_shortest_path_length(G_boundary)
avg_length_ring = average_shortest_path_length(G_ring)
avg_length_clusters = average_shortest_path_length(G_clusters)

# Calculate longest shortest path length (diameter)
longest_length_random = longest_shortest_path_length(G_random)
longest_length_boundary = longest_shortest_path_length(G_boundary)
longest_length_ring = longest_shortest_path_length(G_ring)
longest_length_clusters = longest_shortest_path_length(G_clusters)

print(f"Expected shortest path length (random): {avg_length_random}")
print(f"Longest shortest path length (random): {longest_length_random}\n")

print(f"Expected shortest path length (boundary): {avg_length_boundary}")
print(f"Longest shortest path length (boundary): {longest_length_boundary}\n")

print(f"Expected shortest path length (ring): {avg_length_ring}")
print(f"Longest shortest path length (ring): {longest_length_ring}\n")

print(f"Expected shortest path length (clusters): {avg_length_clusters}")
print(f"Longest shortest path length (clusters): {longest_length_clusters}\n")

# Plot the graphs for visualization
def plot_graph(points, G, title):
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos={i: points[i] for i in range(len(points))}, with_labels=False, node_size=20, node_color='blue')
    plt.title(title)
    plt.show()

plot_graph(points_random, G_random, "Gabriel Graph (Random)")
plot_graph(points_boundary, G_boundary, "Gabriel Graph (Boundary)")
plot_graph(points_ring, G_ring, "Gabriel Graph (Ring)")
plot_graph(points_clusters, G_clusters, "Gabriel Graph (Clusters)")
