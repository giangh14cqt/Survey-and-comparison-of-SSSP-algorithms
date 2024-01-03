import matplotlib.pyplot as plt

sources = [
    {
        "algorithm": "Dijkstra",
        "benchmark_file": "./running-result/DijkstraPerformance2.txt",
    },
    {
        "algorithm": "Bellman-Ford",
        "benchmark_file": "./running-result/BellmanFordPerformance.txt",
    },
    {
        "algorithm": "Floyd-Warshall",
        "benchmark_file": "./running-result/FloydWarshallPerformance2.txt",
    },
    {
        "algorithm": "Wulff-Nilsen",
        "benchmark_file": "./running-result/LasVegasPerformance.txt",
    },
    # {
    #     "algorithm": "Wulff-Nilsen",
    #     "benchmark_file": "./running-result/LasVegasPerformance3.txt",
    # },
]

def compare_each_algorithm_by_density():
    for source in sources:
        f = open(source["benchmark_file"], "r")
        lines = f.readlines()
        f.close()
        x = []
        y = []
        last_line = ''
        data = {}
        for line in lines:
            if len(line.split(' ')) == 1:
                n, m = last_line.split(' ')
                if data.get(n) is None:
                    data[n] = {}
                data[n][m] = float(line)
            last_line = line
        for n, m in data.items():
            x.append(int(n))
            y.append(list(m.values()))
        plt.plot(x, y, label=source["algorithm"])
        # plt.legend(data['100'].keys())
        plt.legend([('density = ' + str(m)) for m in data['100'].keys()])
        plt.xlabel('Number of vertices')
        plt.ylabel('Time (ms)')
        plt.title('Performance of shortest path algorithms ' + source["algorithm"])
        plt.show()

def compare_each_density_by_algorith():
    file_readers = [open(source["benchmark_file"], "r") for source in sources]
    lines_list = [f.readlines() for f in file_readers]
    [f.close() for f in file_readers]
    data = []
    for lines in lines_list:
        data.append({})
        last_line = ''
        for line in lines:
            if len(line.split(' ')) == 1:
                n, m = last_line.split(' ')
                m = int(m)
                n = int(n)
                
                if data[-1].get(m) is None:
                    data[-1][m] = {}
                data[-1][m][n] = float(line)
            last_line = line

    for m in data[0].keys():
        y = []
        x = list(data[0][m].keys())
        # for i in range(len(sources)):
        #     y.append(list(data[i][m].values()))
        for n in data[0][m].keys():
            y.append([data[i][m][n] for i in range(len(sources))])
        plt.plot(x, y, label=m)
        plt.legend([source["algorithm"] for source in sources])
        plt.xlabel('Number of vertices')
        plt.ylabel('Time (ms)')
        plt.title('Performance of shortest path algorithms for density ' + str(m))
        plt.show()

# compare_each_density_by_algorith()
compare_each_algorithm_by_density()