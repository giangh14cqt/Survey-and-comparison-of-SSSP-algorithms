import random
import networkx as nx
import matplotlib.pyplot as plt
import os.path

def generate_random_directed_weighted_graph(n, density_ratio):
    # Create an empty directed graph
    G = nx.DiGraph()

    # Add vertices
    G.add_nodes_from(range(n))
    edge_max = n * (n - 1)  # Maximum number of edges possible
    density = n * density_ratio / edge_max  # Density of the graph
    # Add random weighted edges based on density
    for i in range(n):
        print(i)
        for j in range(n):
            if (i != j and random.random() < density) or (i == j+1):
                # Generate a random weight for the edge
                weight = random.randint(n, 2*n)  # You can adjust the range as needed
                if i == j+1:
                    weight = -1
                # Add the directed edge with the random weight
                G.add_edge(i, j, weight=weight)

    return G

def graph_gen(n_vertices, density, file_name, file_dir):
    if not os.path.isdir(file_dir):
        os.mkdir(file_dir)
    if not os.path.isdir(os.path.join(file_dir, str(n_vertices))):
        os.mkdir(os.path.join(file_dir, str(n_vertices)))
    random_graph = generate_random_directed_weighted_graph(n_vertices, density)
    file_path = os.path.join(file_dir, str(n_vertices), file_name)
    outFile=open(file_path, "w")
    outFile.write(str(n_vertices) + "\n")
    phi = [random.randint(-200, 200) for _ in range(n_vertices)]
    for edge in random_graph.edges():
        u = edge[0]
        v = edge[1]
        w = random_graph.get_edge_data(u, v)['weight']
        # if file_dir.find("negative") != -1:
        #     w += phi[u] - phi[v]
        outFile.write(str(u) + " " + str(v) + " " + str(w) + "\n")

def manual():
    n_vertices = int(input("Enter the number of vertices: "))
    directory = input("Enter the directory: ")
    file_name = input("Enter the file name: ")
    density = float(input("Enter the density (0, 1]: "))
    display_graph = bool(input("Display graph? (skip for no, any character for yes): ") != '')
    if not os.path.isdir(directory):
        os.mkdir(directory)
    random_graph = generate_random_directed_weighted_graph(n_vertices, density)
    file_path = os.path.join(directory, file_name)
    outFile=open(file_path, "w")
    outFile.write(str(n_vertices) + "\n")
    for edge in random_graph.edges():
        outFile.write(str(edge[0]) + " " + str(edge[1]) + " " + str(random_graph.get_edge_data(edge[0], edge[1])['weight']) + "\n")

    if display_graph:
        # Plotting the graph
        pos = nx.circular_layout(random_graph)
        labels = nx.get_edge_attributes(random_graph, 'weight')
        nx.draw(random_graph, pos, with_labels=True)
        nx.draw_networkx_edge_labels(random_graph, pos, edge_labels=labels)
        plt.show()

if __name__ == '__main__':
    # n_vertices = [100, 500, 1000, 2000, 4000, 6000, 8000, 10000]
    # n_vertices = [15000, 20000, 25000, 30000, 35000, 40000]
    # n_vertices = [50000, 55000, 60000, 65000, 70000]
    n_vertices = [10000]
    # densities = [10, 25]#, 50, 75, 100]
    densities = [25]
    # file_dir = ["non-negative-weighted2", "negative-weighted2"]
    file_dir = ["negative-weighted-worst-case"]
    for fd in file_dir:
        for n in n_vertices:
            for d in densities:
                graph_gen(n, d, str(n) + "_" + str(d) + ".txt", fd)

