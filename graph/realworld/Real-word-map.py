import osmnx as ox
import networkx as nx
import random

# Specify the city or region
city_name = "Warsaw Poland"
place_query = f"{city_name}"

# Download the OSM data
graph = ox.graph_from_place(place_query, network_type="drive")

# Visualize the graph
ox.plot_graph(graph, bgcolor='white', edge_color='black', 
              node_color='red', edge_linewidth=0.1, 
              node_size=1, figsize=(10, 10), dpi=300, 
              save=True, filepath=f"{city_name}.png")

# Convert OSM graph to NetworkX graph
G = nx.Graph(graph)

# Add distance as edge weight
for u, v, data in G.edges(data=True):
    data['weight'] = ox.distance.great_circle(G.nodes[u]['y'], G.nodes[u]['x'],
                                                 G.nodes[v]['y'], G.nodes[v]['x'])

# Get graph information
num_vertices = G.number_of_nodes()
num_edges = G.number_of_edges()

print(f"Number of Vertices: {num_vertices}")
print(f"Number of Edges: {num_edges}")

mapVertices = {}
i = 0
for v in G.nodes():
    mapVertices[v] = i
    i += 1

# Save the graph with number of vertices and edges and edge weights
f = open(city_name + ".txt", "w")
f.write(f"{num_vertices}\n")
for u, v, data in G.edges(data=True):
    f.write(f"{mapVertices[u]} {mapVertices[v]} {int(data['weight'])}\n")

# Save graph with negative edge weights
f = open(city_name + "N.txt", "w")
f.write(f"{num_vertices}\n")
phi = [random.randint(-200, 200) for _ in range(num_vertices)]
for u, v, data in G.edges(data=True):
    f.write(f"{mapVertices[u]} {mapVertices[v]} {int(data['weight']) - phi[mapVertices[u]] + phi[mapVertices[v]]}\n")