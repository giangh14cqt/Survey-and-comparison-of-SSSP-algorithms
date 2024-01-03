#include <iostream>
#include <limits>
#include <vector>
#include <fstream>

using namespace std;

const int INF = numeric_limits<int>::max(); // Set initial distance to infinity

class Graph
{
private:
    int V; // Number of vertices
    vector<vector<int>> dist;

public:
    // Constructor
    Graph(int V)
    {
        this->V = V;
        dist.resize(V, vector<int>(V, INF));
        // Initialize diagonal elements to 0
        for (int i = 0; i < V; i++)
        {
            dist[i][i] = 0;
        }
    }

    // Add a directed edge from node u to node v with weight w
    void add_edge(int u, int v, int w) { dist[u][v] = w; }

    Graph(ifstream &inputFile)
    {
        inputFile >> V;
        dist.resize(V, vector<int>(V, INF));
        int u, v, w;
        while (inputFile >> u >> v >> w)
        {
            add_edge(u, v, w);
        }
    }

    // Execute the Floyd-Warshall algorithm
    void floyd_warshall()
    {
        for (int k = 0; k < V; k++)
            for (int i = 0; i < V; i++)
                if (dist[i][k] != INF)
                    for (int j = 0; j < V; j++)
                        if (dist[k][j] != INF)
                            dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
    }

    // Return shortest distance between two vertices
    int get_shortest_distance(int u, int v) { return dist[u][v]; }
};
