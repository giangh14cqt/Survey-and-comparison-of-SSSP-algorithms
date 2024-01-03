import subprocess
import os.path

# binary_file_path = "../src/LasVegas/build/SSSP_New"
binary_file_path = {
    'Dijkstra': "../src/Dijkstra/build/Dijkstra",
    'BellmanFord': "../src/BellmanFord/build/BellmanFord",
    'FloydWarshall': "../src/Floyd-Warshall/build/FloydWarshall",
}


def benchmark(algorithm, n_vertices, density, f):
    file_name = str(n_vertices) + "_" + str(density) + ".txt"
    file_path = os.path.join("../graph/non-negative-weighted2", str(n_vertices), file_name)
    f.write(str(n_vertices) + " " + str(density) + "\n")
    # Run the binary file and capture the output
    running_file = binary_file_path.get(algorithm)
    print(algorithm, n_vertices, density)
    result = subprocess.run(
        [running_file, file_path], stdout=subprocess.PIPE, text=True
    )

    # Access the output
    output = result.stdout
    # write the output to a file
    f.write(output)


if __name__ == "__main__":
    n_vertices = [100, 500, 1000, 2000, 4000, 6000, 8000, 10000]
    densities = [10, 25, 50, 75, 100]
    file_dir = ["non-negative-weighted"]
    for k, v in binary_file_path.items():
        f = open(k + "Performance2.txt", "w")
        for n in n_vertices:
            for d in densities:
                benchmark(k, n, d, f)

        f.close()
