//
// Created by Truong Giang Do on 28/11/2023.
//

#ifndef SSSP_NEW_GRAPH_H
#define SSSP_NEW_GRAPH_H


#include <iostream>
#include <queue>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <algorithm>
#include <map>
#include <stack>
#include <set>
#include <fstream>
#include <unordered_set>

#include "Timer.h"

using namespace std;

class Node {
public:
    int node;
    int cost;

    Node(int n, int c) {
        node = n;
        cost = c;
    }

    bool operator<(const Node &n) const {
        return cost > n.cost;
    }

    bool operator>(const Node &n) const {
        return cost < n.cost;
    }

    bool operator==(const Node &n) const {
        return cost == n.cost;
    }
};

class Graph {
public:
    int v_max;            // max number of vertices
    vector<int> vertices; // specifies the vertices in the graph
    vector<bool> containsVertex;

    int n;                             // number of vertices
    vector<vector<int>> adjacencyList; // adjacency list of the graph
    vector<vector<int>> weights;       // weights of the edges

    int time; // Used for SCCs by Tarjan's algorithm

    Graph(const Graph &g, bool addDummy = false) {
        if (!addDummy) {
            v_max = g.v_max;
            vertices = g.vertices;
            containsVertex = g.containsVertex;
            n = g.n;
            adjacencyList = g.adjacencyList;
            weights = g.weights;
            time = g.time;
        } else {
            v_max = g.v_max + 1;
            vertices = g.vertices;
            containsVertex = g.containsVertex;
            n = g.n + 1;
            adjacencyList = g.adjacencyList;
            weights = g.weights;
            time = 0;

            vertices.push_back(g.v_max);
            containsVertex.push_back(true);
            adjacencyList[g.v_max].reserve(g.n);
            weights[g.v_max].reserve(g.n);
            for (int i = 0; i < g.n; ++i) {
                adjacencyList[g.v_max].emplace_back(g.vertices[i]);
                weights[g.v_max].emplace_back(0);
            }
        }
    }

    Graph(int numVertices, bool withAllVertices) {
        v_max = numVertices;
        vertices = vector<int>();
        containsVertex = vector<bool>(v_max, false);

        if (withAllVertices) {
            vertices = vector<int>(numVertices);
            for (int v = 0; v < numVertices; v++) {
                vertices[v] = v;
                containsVertex[v] = true;
            }
            n = numVertices;
        } else
            n = 0;
        adjacencyList = vector<vector<int>>(numVertices);
        weights = vector<vector<int>>(numVertices);

        time = 0;
    }

    void addVertices(vector<int> &vertices_) {
        for (int &v: vertices_) {

            if (0 <= v && v < v_max && !containsVertex[v]) {
                containsVertex[v] = true;
            } else
                throw_with_nested("Vertex out of bounds or already in graph");
        }
        vertices.reserve(vertices.size() + vertices_.size());
        vertices.insert(vertices.end(), vertices_.begin(), vertices_.end());
        n += vertices_.size();
    }

    void addVertex(int v) {
        if (0 <= v && v < v_max && !containsVertex[v]) {
            vertices.push_back(v);
            containsVertex[v] = true;
            n++;
        } else
            throw_with_nested("Vertex out of bounds or already in graph");
    }

    void addEdges(int v, vector<int> &outVertices, vector<int> &outWeights) {
        if (outVertices.size() != outWeights.size())
            throw_with_nested("Number of out vertices and weights must be equal");
        adjacencyList[v] = vector<int>(outVertices.begin(), outVertices.end());
        weights[v] = vector<int>(outWeights.begin(), outWeights.end());
    }

    void removeVertices(vector<int> &ball) {
        for (int v: ball) {
            if (containsVertex[v]) {
                containsVertex[v] = false;
                n--;
            }
        }
        for (int v: vertices)
            if (!containsVertex[v]) {
                adjacencyList[v].clear();
                weights[v].clear();
            } else {
                vector<int> newOutVertices;
                vector<int> newOutWeights;
                for (int i = 0; i < adjacencyList[v].size(); ++i) {
                    int outV = adjacencyList[v][i];
                    if (containsVertex[outV]) {
                        newOutVertices.push_back(outV);
                        newOutWeights.push_back(weights[v][i]);
                    }
                }
                adjacencyList[v] = newOutVertices;
                weights[v] = newOutWeights;
            }
        vector<int> newVertices(n);
        int i = 0;
        for (int v: vertices) {
            if (containsVertex[v]) {
                newVertices[i] = v;
                i++;
            }
        }
        vertices = newVertices;
    }

    void displayGraph() {
        cout << "Displaying graph:" << endl;
//        for (int i = 0; i < v_max; i++) {
//            if (containsVertex[i]) {
//                cout << "Vertex " << i << " has out edges: ";
//                for (unsigned long j = 0; j < adjacencyList[i].size(); j++) {
//                    cout << adjacencyList[i][j] << " with weight " << weights[i][j] << ", ";
//                }
//                cout << endl;
//            }
//        }
        for (int v: vertices) {
            for (unsigned long j = 0; j < adjacencyList[v].size(); j++) {
                cout << v << ' ' << adjacencyList[v][j] << endl;
            }
        }
    }

    // Run Tarjan's algorithm to find SCCs
    vector<vector<int>> SCC() {
        vector<vector<int>> SCCverts;

        vector<int> vertsToN = vector<int>(v_max, 0);
        vector<int> NtoVerts = vector<int>(n, 0);

        for (unsigned long i = 0; i < vertices.size(); ++i) {
            vertsToN[vertices[i]] = i;
            NtoVerts[i] = vertices[i];
        }

        vector<int> disc = vector<int>(n, 0);
        vector<int> low = vector<int>(n, 0);
        time = 0;

        vector<bool> stackMember = vector<bool>(n, false);
        stack<int> st;

        for (int i = 0; i < n; ++i) {
            if (disc[i] == -1) {
//                vector<vector<int>> SCCTmp = SCCUtil(i, low, disc, stackMember, st, vertsToN, NtoVerts);
//                for (const vector<int> &SCC: SCCTmp)
//                    SCCverts.push_back(SCC);
                dfs(i, low, disc, stackMember, st, vertsToN, NtoVerts, SCCverts);
            }
        }

        return SCCverts;
    }

    void dfs(int u, vector<int> &low, vector<int> &disc, vector<bool> &stackMember, stack<int> &st,
             vector<int> &vertsToN, vector<int> &NtoVerts, vector<vector<int>> &SCCverts) {
        disc[u] = low[u] = ++time;
        st.push(u);
        for (int v_true: adjacencyList[NtoVerts[u]]) {
            int v = vertsToN[v_true];
            if (stackMember[v]) continue;
            if (!disc[v]) {
                dfs(v, low, disc, stackMember, st, vertsToN, NtoVerts, SCCverts);
                low[u] = min(low[u], low[v]);
            } else low[u] = min(low[u], disc[v]);
        }
        if (low[u] == disc[u]) {
            vector<int> SCC;
            int v;
            do {
                v = st.top();
                st.pop();
                stackMember[v] = true;
                SCC.push_back(v);
            } while (v != u);
            SCCverts.push_back(SCC);
        }
    }

    vector<vector<int>> SCCUtil(int u, vector<int> &low, vector<int> &disc, vector<bool> &stackMember, stack<int> &st,
                                vector<int> &vertsToN, vector<int> &NtoVerts) {

//        vector<vector<int>> SCCverts;
//        disc[u] = low[u] = time;
//        time++;
//        stackMember[u] = true;
//        st.push(u);
//
//        for (int v_true: adjacencyList[NtoVerts[u]]) {
//            int v = vertsToN[v_true];
//            if (disc[v] == -1) {
//                vector<vector<int>> SCCTmp = SCCUtil(v, low, disc, stackMember, st, vertsToN, NtoVerts);
//                SCCverts.insert(SCCverts.end(), SCCTmp.begin(), SCCTmp.end());
//                low[u] = min(low[u], low[v]);
//            } else if (stackMember[v])
//                low[u] = min(low[u], disc[v]);
//        }
//
//        int w = -1;
//        if (low[u] == disc[u]) {
//            vector<int> SCC = vector<int>();
//            while (w != u) {
//                w = st.top();
//                st.pop();
//                SCC.push_back(NtoVerts[w]);
//                stackMember[w] = false;
//            }
//            SCCverts.push_back(SCC);
//        }
//        return SCCverts;
    }

    // Returns true if the graph has a negative weight cycle.
    // Assumes every vertex is being used.
    bool hasNegCycle() {
        vector<int> vertsToN = vector<int>(v_max, -1);
        vector<int> NtoVerts = vector<int>(n, -1);

        for (unsigned long i = 0; i < vertices.size(); ++i) {
            vertsToN[vertices[i]] = i;
            NtoVerts[i] = vertices[i];
        }

        vector<bool> visited = vector<bool>(n, false);
        vector<int> dist = vector<int>(n, INT_MAX);

        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                if (isNegCycleBellmanFord(i, dist, vertsToN, NtoVerts))
                    return true;
            }

            for (int j = 0; j < n; ++j)
                if (dist[j] != INT_MAX)
                    visited[j] = true;
        }

        return false;
    }

    bool isNegCycleBellmanFord(int src, vector<int> &dist, vector<int> &vertsToN, vector<int> &NtoVerts) {
        for (int i = 0; i < n; ++i)
            dist[i] = INT_MAX;
        dist[src] = 0;

        for (int i = 1; i <= n - 1; ++i) {
            for (int u = 0; u < n; ++u) {
                int u_real = NtoVerts[u];
                for (unsigned long j = 0; j < adjacencyList[u_real].size(); ++j) {
                    int v_real = adjacencyList[u_real][j];
                    int v = vertsToN[v_real];
                    int weight = weights[u_real][j];
                    if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
                        dist[v] = dist[u] + weight;
                }
            }
        }

        for (int u = 0; u < n; ++u) {
            int u_real = NtoVerts[u];
            for (unsigned long j = 0; j < adjacencyList[u_real].size(); ++j) {
                int v_real = adjacencyList[u_real][j];
                int v = vertsToN[v_real];
                int weight = weights[u_real][j];
                if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
                    return true;
            }
        }

        return false;
    }

    void BellmanFord(int src_real) {
        vector<int> vertsToN = vector<int>(v_max, -1);
        vector<int> NtoVerts = vector<int>(n, -1);

        for (unsigned long i = 0; i < vertices.size(); ++i) {
            vertsToN[vertices[i]] = i;
            NtoVerts[i] = vertices[i];
        }

        vector<int> dist = vector<int>(n, INT_MAX);

        int src = vertsToN[src_real];
        dist[src] = 0;

        for (int i = 1; i < n; ++i) {
            for (int u = 0; u < n; ++u) {
                int u_real = NtoVerts[u];
                for (unsigned long j = 0; j < adjacencyList[u_real].size(); ++j) {
                    int v_real = adjacencyList[u_real][j];
                    int v = vertsToN[v_real];
                    int weight = weights[u_real][j];
                    if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
                        dist[v] = dist[u] + weight;
                }
            }
        }
    }

    bool hasNoNegativeEdgeWeights() {
        for (int v: vertices) {
            for (unsigned long i = 0; i < adjacencyList[v].size(); ++i) {
                if (weights[v][i] < 0)
                    return false;
            }
        }
        return true;
    }
};

#endif //SSSP_NEW_GRAPH_H
