#include <iostream>
#include <limits>
#include <vector>

using namespace std;

const int INF = numeric_limits<int>::max(); // Set initial distance to infinity

class Graph
{
private:
    int V;                                   // Number of vertices
    vector<vector<pair<int, int>>> adj_list; // Adjacency list to store the graph
    vector<int> dist;

public:
    // Constructor
    Graph(int V)
    {
        this->V = V;
        adj_list.resize(V);
    }

    // Add a directed edge from node u to node v with weight w
    void add_edge(int u, int v, int w) { adj_list[u].push_back({v, w}); }

    Graph(ifstream &inputFile)
    {
        inputFile >> V;
        adj_list.resize(V);
        int u, v, w;
        while (inputFile >> u >> v >> w)
        {
            add_edge(u, v, w);
        }
    }

    // Execute the Bellman-Ford algorithm from source
    bool bellman_ford(int source)
    {
        dist.assign(V, INF); // Distance from source to each vertex
        dist[source] = 0;
        for (int i = 1; i < V; i++)
        {
            for (int u = 0; u < V; u++)
            {
                for (auto iter : adj_list[u])
                {
                    int v = iter.first;
                    int w = iter.second;
                    if (dist[u] != INF && dist[u] + w < dist[v])
                    {
                        dist[v] = dist[u] + w;
                    }
                }
            }
        }
        // Check for negative-weight cycles
        for (int u = 0; u < V; u++)
        {
            for (auto iter : adj_list[u])
            {
                int v = iter.first;
                int w = iter.second;
                if (dist[u] != INF && dist[u] + w < dist[v])
                {
                    return false; // Negative-weight cycle exists
                }
            }
        }
        return true; // No negative-weight cycle found
    }

    // Return shortest distance to a vertex
    int get_shortest_distance(int vertex) { return dist[vertex]; }
};
