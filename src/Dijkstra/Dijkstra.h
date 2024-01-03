#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <iostream>
#include <queue>
#include <vector>
#include <fstream>

using namespace std;

const int INF = INT_MAX; // Set initial distance to infinity

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

    // Execute the Dijkstra algorithm from source
    void dijkstra(int source)
    {
        dist.assign(V, INF);        // Distance from source to each vertex
        vector<bool> vis(V, false); // Visited flag for each vertex
        dist[source] = 0;
        priority_queue<pair<int, int>, vector<pair<int, int>>,
                       greater<pair<int, int>>>
            pq;
        pq.push(make_pair(0, source));
        while (!pq.empty())
        {
            int u = pq.top().second;
            pq.pop();
            if (vis[u])
                continue;
            vis[u] = true;
            for (auto iter : adj_list[u])
            {
                auto v = iter.first;
                auto w = iter.second;
                if (dist[u] + w < dist[v])
                {
                    dist[v] = dist[u] + w;
                    pq.push(make_pair(dist[v], v));
                }
            }
        }
    }

    // Return shortest distace to a vertex
    int get_shortest_distance(int vertex) { return dist[vertex]; }
};

#endif