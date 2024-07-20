import networkx as nx
import matplotlib.pyplot as plt

# Edge list provided
edges = [
    (7, 9), (25, 27), (10, 13), (26, 29), (13, 26), (9, 23), (22, 25), (5, 26), (13, 29), (12, 27),
    (2, 23), (9, 23), (13, 30), (24, 25), (6, 13), (12, 17), (1, 19), (2, 9), (12, 22), (10, 26),
    (8, 24), (21, 28), (20, 21), (1, 21), (18, 21), (18, 28), (5, 18), (14, 18), (20, 27), (20, 23),
    (20, 24), (8, 20), (1, 20), (1, 19), (1, 8), (1, 28), (19, 28), (11, 28), (5, 28), (8, 24),
    (8, 19), (5, 15), (5, 11), (5, 26), (5, 14), (14, 15), (19, 30), (11, 19), (19, 24), (24, 30),
    (9, 24), (23, 24), (11, 30), (11, 26), (23, 27), (7, 23), (22, 23), (9, 23), (26, 30), (4, 26),
    (10, 26), (15, 26), (9, 29), (9, 30), (9, 22), (7, 22), (22, 25), (22, 29), (13, 30), (29, 30),
    (4, 30), (7, 25), (7, 27), (3, 27), (25, 27), (4, 13), (4, 10), (16, 29), (2, 29), (13, 29),
    (25, 29), (6, 15), (10, 15), (3, 25), (17, 25), (16, 25), (10, 13), (6, 10), (16, 17), (12, 16),
    (2, 16), (2, 6), (2, 12), (2, 13), (6, 13), (12, 17), (6, 12), (3, 17), (6, 17), (3, 6)
]
#15 10 27 40 40 1
path = '46 28 25 24 33 30 12 1 72 84 74 10 19 68 78 89 99 15 76 71 50 21 41 37 38 54 61'.split(' ')
path = list(map(int, path))
path = [i - 1 for i in path]
print(path)


# Create a graph
G = nx.Graph()
G.add_edges_from(edges)
edge_color = []
for u,v in G.edges():
    if edges.count((u,v)) > 0:
        if path.count(edges.index((u,v))) > 0:
            edge_color.append('blue')
            continue
    if edges.count((v,u)) > 0:
        if path.count(edges.index((v,u))) > 0:
            edge_color.append('blue')
            continue
        
    edge_color.append('gray')




node_colors = ['red' if node == 15 or node == 10 else 'blue' for node in G.nodes()]
# Draw the graph
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, seed=42)  # positions for all nodes
nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=700, edge_color=edge_color)
plt.title("Graph Visualization")
plt.show()


print(edges[46-1])
print(edges[28-1])
print(edges[25-1])