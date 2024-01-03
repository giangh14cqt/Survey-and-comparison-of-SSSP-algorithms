import subprocess
import os.path

binary_file_path = {
    'LasVegas': "../src/LasVegas/build/SSSP_New",
    'BellmanFord': "../src/BellmanFord/build/BellmanFord"
}

def benchmark(algorithm, n_vertices, density, f):
    file_name = str(n_vertices) + "_" + str(density) + ".txt"
    file_path = os.path.join("../graph/negative-weighted2", str(n_vertices), file_name)
    f.write(str(n_vertices) + " " + str(density) + "\n")
    # Run the binary file and capture the output
    running_file = binary_file_path.get(algorithm)
    print([running_file, file_path])
    result = subprocess.run(
        [running_file, file_path], stdout=subprocess.PIPE, text=True
    )
    # Access the output
    output = result.stdout
    print(output)
    # write the output to a file
    f.write(output)


if __name__ == "__main__":
    # n_vertices = [100, 500, 1000, 2000, 4000, 6000, 8000, 10000, 15000, 20000, 25000, 30000, 35000, 40000]
    # n_vertices = [15000, 20000, 25000, 30000, 35000, 40000]
    n_vertices = [45000, 50000, 55000, 60000, 65000, 70000]
    densities = [10, 25]#, 50, 75, 100]
    for k, v in binary_file_path.items():
        f = open(k + "Performance.txt", "w")
        for n in n_vertices:
            for d in densities:
                benchmark(k, n, d, f)

        f.close()
