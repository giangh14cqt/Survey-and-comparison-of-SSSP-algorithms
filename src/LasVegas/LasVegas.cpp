//
// Created by Truong Giang Do on 01/12/2023.
//
#include "LasVegas.h"

bool CHECKS = false;
int SRC = 0;
bool WITH_LDD = true;

Graph readInput(ifstream &inputFile) {
    int g_size;
    inputFile >> g_size; // >> maxWeight;
    ++g_size;

    Graph g(g_size, false);
    vector<vector<int>> edges(g_size);
    vector<vector<int>> weights(g_size);

    vector<vector<bool>> edge_exists(g_size, vector<bool>(g_size, false));
    int u, v, w;
    while (inputFile >> u >> v >> w) {
        if (u > g_size || v > g_size) {
            exit(1);
        }
        if (!g.containsVertex[u])
            g.addVertex(u);

        if (!g.containsVertex[v])
            g.addVertex(v);

        if (!edge_exists[u][v]) {
            edges[u].push_back(v);
            weights[u].push_back(w);
            edge_exists[u][v] = true;
        }
    }

    for (int i = 0; i < g_size; ++i)
        g.addEdges(i, edges[i], weights[i]);
    return g;
}

vector<int> lasVegas(Graph &g) {
    gn_global = g.n;
    Timer::startTimer();
    int minWeight = INT_MAX;
    for (int u: g.vertices)
        for (unsigned long i = 0; i < g.adjacencyList[u].size(); ++i)
            minWeight = min(minWeight, g.weights[u][i]);

    if (minWeight >= 0) {
        vector<int> dist = Dijkstra(g, SRC);
        cout << Timer::getDuration() << endl;
        return dist;
    }
    vector<int> phi(g.v_max);
    while (minWeight < -1) {
        Graph gScaledS(g, false);

        for (int u: g.vertices) {
            for (unsigned long i = 0; i < g.adjacencyList[u].size(); ++i) {
                gScaledS.weights[u][i] = phi[u] - phi[g.adjacencyList[u][i]] + g.weights[u][i];
            }
        }

        vector<int> phi_i = ScaleDown(gScaledS, gScaledS.n, -1 * minWeight / 2);

        for (int u: g.vertices)
            phi[u] += phi_i[u];

        minWeight = INT_MAX;
        for (int u: g.vertices)
            for (unsigned long i = 0; i < g.adjacencyList[u].size(); ++i)
                minWeight = min(minWeight, g.weights[u][i] + phi[u] - phi[g.adjacencyList[u][i]]);
    }

    Graph gFinal(g);

    vector<int> tree = SPMain(gFinal, SRC);

    vector<int> dist = getDistFromTree(g, tree);
    cout << Timer::getDuration() << endl;
    return dist;
}

vector<int> SPMain(Graph &g_in, int s) {
    int scaleFactor = 2 * g_in.n;
    Graph g = getScaledGraph(g_in, scaleFactor);
    int B = roundPower2(scaleFactor);
    vector<int> phi(g.v_max);
    set<vector<int>> nullErem;
    for (int i = 1; i <= logBase2(B); i++) {
        nullErem.clear();
        Graph g_phi = createModifiedGB(g, 0, false, nullErem, phi);
        vector<int> phi_i = ScaleDown(g_phi, g.n, B / (int) pow(2, i));

        if (CHECKS && hasNegativeEdges(g_phi, phi_i, B / (int) pow(2, i))) {
            throw_with_nested("ScaleDown failed.");
        }

        phi = addPhi(phi, phi_i);
    }

    // create G^*
    for (int u: g.vertices) {
        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            g.weights[u][i] += phi[u] - phi[g.adjacencyList[u][i]] + 1;

            if (CHECKS && g.weights[u][i] < 0) {
                throw_with_nested(
                        "After applying the phi outputted from SPMain, there exists an edge with negative weight.");
            }
        }
    }

    vector<int> tree = getShortestPathTree(g, s);

    return tree;
}

Graph getScaledGraph(Graph &g_in, int scaleFactor) {
    Graph g = Graph(g_in.v_max, false);
    g.addVertices(g_in.vertices);

    for (int u: g_in.vertices) {
        vector<int> edges(g_in.adjacencyList[u].size());
        vector<int> weights(g_in.adjacencyList[u].size());

        for (unsigned long i = 0; i < g_in.adjacencyList[u].size(); ++i) {
            edges[i] = g_in.adjacencyList[u][i];
            weights[i] = g_in.weights[u][i] * scaleFactor;
        }
        g.addEdges(u, edges, weights);
    }

    return g;
}

int roundPower2(int n) {
    return int(pow(2, ceil(logBase2(n))));
}

double logBase2(int n) {
    return log(n) / log(2);
}

/*
 * 1. INPUT REQUIREMENTS:
 * (a) B is positive integer, w is integral, and w(e) ≥ −2B for all e ∈ E
 * (b) If the graph G does not contain a negative-weight cycle then the input
 * must satisfy η(GB) ≤ ∆; that is, for every v ∈ V there is a shortest
 * sv-path in GBs with at most ∆ negative edges
 * (c) All vertices in G have constant out-degree
 * 2. OUTPUT: If it terminates, the algorithm returns an integral price function
 * φ
 * such that wφ(e) ≥ −B for all e ∈ E
 * 3. RUNNING TIME: If G does not contain a negative-weight cycle, then the
 * algorithm has epected runtime O(m log3(n) log(∆)).
 * Remark: If G contains a negative-weight cycle, there is no guarantee
 * on the runtime, and the algorithm might not even terminate; but if the
 * algorithm does terminate, it always produces a correct output.
 */
vector<int> ScaleDown(Graph &g, int delta, int B) {
    vector<int> phi_2(g.v_max);
    vector<int> emptyPhi;
    set<vector<int>> emptyRem;

    if (delta > 2) {
        double d = delta / 2.0;
        Graph g_b_nneg = createModifiedGB(g, B, true, emptyRem, emptyPhi);

        // phase 0
        vector<vector<int>> E_sep;
        if (WITH_LDD)
            E_sep = LDDRework(g_b_nneg, int(d * B));

        set<vector<int>> E_sep_hash(E_sep.begin(), E_sep.end());
        Graph g_B_Esep = createModifiedGB(g, B, false, E_sep_hash, emptyPhi);
        vector<vector<int>> SCCs = g_B_Esep.SCC();

        // phase 1
        vector<int> vertexToSCCMap = getVertexToSCCMap(SCCs, g.v_max);
        set<vector<int>> edgesBetweenSCCs = getEdgesBetweenSCCs(g, vertexToSCCMap);
        Graph H = createModifiedGB2(g, 0, false, edgesBetweenSCCs, emptyPhi);
        vector<int> phi_1 = ScaleDown(H, delta / 2, B);
        // phase 2
        Graph g_B_E_sep_phi1 = createModifiedGB(g, B, false, E_sep_hash, phi_1);
        vector<int> phi = FixDAGEdges(g_B_E_sep_phi1, SCCs, vertexToSCCMap, edgesBetweenSCCs);
        phi_2 = addPhi(phi_1, phi);

        if (CHECKS && hasNegativeEdges(g_B_Esep, phi_2, 0))
            throw_with_nested("FixDAGEdges failed.");
    }
    // phase 3
    Graph g_B_phi2 = createModifiedGB(g, B, false, emptyRem, phi_2);
    vector<int> phi_prime = ElimNeg(g_B_phi2);
    vector<int> phi_3 = addPhi(phi_2, phi_prime);

    if (CHECKS) {
        if (hasNegativeEdges(g, phi_3, B))
            throw_with_nested("ElimNeg failed.");
    }

    return phi_3;
}

bool hasNegativeEdges(Graph &g, vector<int> &phi, int B) {
    for (int u: g.vertices) {
        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            if (g.weights[u][i] + phi[u] - phi[g.adjacencyList[u][i]] < -1 * B) {
                return true;
            }
        }
    }

    return false;
}

/*
 * Creates G^B_phi = (V, E, w^B_phi), where
 * w^B_phi(e) = w(e) + phi(u) - phi(v) if w(e) >= 0,
 * and w(e) + B + phi(u) - phi(v) if w(e) < 0.
 * If nneg == true, w(e) = max{0, w^B(e)}.
 * Removes all the edges in remEdges.
 */
Graph createModifiedGB(Graph &g, int B, bool nneg, set<vector<int>> &remEdges, vector<int> &phi) {
    Graph modG(g.v_max, false);
    modG.addVertices(g.vertices);

    for (int u: g.vertices) {
        vector<int> edges, weights;
        edges.reserve(g.adjacencyList[u].size());
        weights.reserve(g.adjacencyList[u].size());

        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            vector<int> edge = {u, v};
            if ((remEdges.empty()) || (remEdges.find(edge) == remEdges.end())) {
                int weight = g.weights[u][i];

                if (weight < 0) {
                    weight += B;
                }
                if (nneg) {
                    weight = max(0, weight);
                }
                if (!phi.empty()) {
                    weight += phi[u] - phi[v];
                }

                edges.push_back(v);
                weights.push_back(weight);
            }
        }

        modG.addEdges(u, edges, weights);
    }

    return modG;
}

Graph createModifiedGB2(Graph &g, int B, bool nneg, set<vector<int>> &remEdges, vector<int> &phi) {
    Graph modG(g.v_max, false);
    modG.addVertices(g.vertices);

    for (int u: g.vertices) {
        vector<int> edges, weights;
        edges.reserve(g.adjacencyList[u].size());
        weights.reserve(g.adjacencyList[u].size());

        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            vector<int> edge = {u, v};
            if (remEdges.find(edge) != remEdges.end()) {
                int weight = g.weights[u][i];

                if (weight < 0) {
                    weight += B;
                }
                if (nneg) {
                    weight = max(0, weight);
                }
                if (!phi.empty()) {
                    weight += phi[u] - phi[v];
                }

                edges.push_back(v);
                weights.push_back(weight);
            }
        }

        modG.addEdges(u, edges, weights);
    }

    return modG;
}

set<vector<int>> getEdgesBetweenSCCs(Graph &g, vector<int> &vertexToSCCMap) {
    set<vector<int>> edgesBetweenSCCs;
    for (int u: g.vertices) {
        for (int v: g.adjacencyList[u]) {
            if (vertexToSCCMap[u] != vertexToSCCMap[v]) {
                edgesBetweenSCCs.insert({vertexToSCCMap[u], vertexToSCCMap[v]});
            }
        }
    }

    return edgesBetweenSCCs;
}

vector<int> getVertexToSCCMap(vector<vector<int>> &SCCs, int numVertices) {
    vector<int> vertexToSCCMap(numVertices, 0);
    for (unsigned long i = 0; i < SCCs.size(); i++) {
        for (int v: SCCs[i]) {
            vertexToSCCMap[v] = i;
        }
    }

    return vertexToSCCMap;
}

vector<int> addPhi(vector<int> &phi_1, vector<int> &phi_2) {
    if (phi_1.size() != phi_2.size())
        throw_with_nested("addPhi: phi_1 and phi_2 must have the same length.");

    vector<int> phi(phi_1.size());
    for (unsigned long i = 0; i < phi_1.size(); i++) {
        phi[i] = phi_1[i] + phi_2[i];
    }
    return phi;
}

vector<int> FixDAGEdges(Graph &g,
                        vector<vector<int>> &SCCs,
                        vector<int> &vertexToSCCMap,
                        set<vector<int>> &edgesBetweenSCCs) {
    int n = SCCs.size();
    vector<vector<int>> SCCAdjList = createSCCAdjList(SCCs, vertexToSCCMap, edgesBetweenSCCs);
    vector<int> topOrdering = topSort(n, SCCAdjList);

    vector<int> mu(n, 0); // indices are in topological order (e.g., index 0 corresponds to the first SCC
    // in topological order)
    for (int u: g.vertices) {
        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            int SCCu = vertexToSCCMap[u];
            int SCCv = vertexToSCCMap[v];
            int edgeWeight = g.weights[u][i];

            if ((SCCu != SCCv) && edgeWeight < mu[topOrdering[SCCv]]) {
                mu[topOrdering[SCCv]] = edgeWeight;
            }
        }
    }

    if (CHECKS) {
        for (vector<int> edge: edgesBetweenSCCs) {
            int u = edge[0];
            int v = edge[1];
            if (topOrdering[u] > topOrdering[v]) {
                throw_with_nested("FixDAGEdges went wrong.");
            }
        }
    }

    vector<int> phi(g.v_max);
    int m = 0;
    for (int j = 1; j < n; j++) {
        m += mu[j];
        for (int v: SCCs[topOrdering[j + n]]) {
            phi[v] = m;
        }
    }

    return phi;
}

// returns the adjacency list for the DAG where every SCC is viewed as a single vertex
vector<vector<int>> createSCCAdjList(vector<vector<int>> &SCCs,
                                     vector<int> &vertexToSCCMap,
                                     set<vector<int>> &edgesBetweenSCCs) {
    vector<vector<int>> SCCAdjList(SCCs.size());

    for (vector<int> edge: edgesBetweenSCCs) {
        int u = edge[0];
        int v = edge[1];
        SCCAdjList[u].push_back(v);
    }

    return SCCAdjList;
}

/*
 * Input: DAG
 * Calculates a topological ordering of the n vertices in the DAG
 * Returns a map of vertex v to its index i in the ordering
 * The index in the order i is also mapped to the vertex v (bijective mapping);
 * however, to prevent duplicate keys, instead of using the key i we use the
 * key i + n.
 */
vector<int> topSort(int n, vector<vector<int>> &adjList) {
    vector<int> topOrdering(2 * n);
    stack<int> stack;
    vector<bool> visited(n, false);

    vector<int> indegree(n, 0);
    queue<int> listSource;

    for (int u = 0; u < n; ++u)
        for (auto v: adjList[u])
            indegree[v]++;

    for (int u = 0; u < n; ++u)
        if (!indegree[u]) listSource.push(u);

    while (!listSource.empty()) {
        int u = listSource.front();
        listSource.pop();
        stack.push(u);
        for (auto v: adjList[u])
            if (--indegree[v] == 0)
                listSource.push(v);
    }

    if (stack.size() < n) {
        throw_with_nested("topSort went wrong: contains cycles");
    }

    int i = 0;
    while (!stack.empty()) {
        int v = stack.top();
        stack.pop();
        topOrdering[v] = i;
        topOrdering[i + n] = v;
        i++;
    }

    return topOrdering;
}

/*
 * ElimNeg takes as input a graph G = (V,E,w) in which all vertices have
 * constant out-degree.
 * The algorithm outputs a price function φ such that w_φ(e) ≥ 0 for all e ∈ E
 * and has running time O(log(n)*(n + sum_v∈V ηG(v)));
 * note that if G contains a negative-weight cycle then v∈V ηG(v) = ∞
 * so the algorithm will never terminate and hence not produce any output.
 */
vector<int> ElimNeg(Graph &g) {
    Graph Gs = createGs(g);
    vector<int> dist(Gs.v_max);
    int s = Gs.v_max - 1;

    dist[s] = 0;
    for (int v = 0; v < s; v++) {
        dist[v] = INT_MAX;
    }

    custom_priority_queue<Node> pq(Gs.v_max);
    vector<bool> inPQ(Gs.v_max);
    pq.emplace(s, dist[s]);
    inPQ[s] = true;

    set<int> marked;
    while (!pq.empty()) {
        // Dijkstra Phase
        while (!pq.empty()) {
            int v = pq.top().node;
            int w = pq.top().cost;
            pq.pop();
            if (w > dist[v]) {
                continue;
            }

            inPQ[v] = false;

            marked.emplace(v);

            for (unsigned long i = 0; i < Gs.adjacencyList[v].size(); i++) {
                int x = Gs.adjacencyList[v][i];
                int edgeWeight = Gs.weights[v][i];

                if (edgeWeight >= 0 && (dist[v] + edgeWeight < dist[x])) {
                    marked.emplace(x);

                    dist[x] = dist[v] + edgeWeight;
                    pq.emplace(Node(x, dist[x]));
                    inPQ[x] = true;
                }
            }
        }

        // Bellman-Ford Phase
        for (int v: marked) {
            for (unsigned long i = 0; i < Gs.adjacencyList[v].size(); i++) {
                int x = Gs.adjacencyList[v][i];
                int edgeWeight = Gs.weights[v][i];

                if (edgeWeight < 0 && (dist[v] + edgeWeight < dist[x])) {
                    dist[x] = dist[v] + edgeWeight;
                    pq.emplace(Node(x, dist[x]));
                    inPQ[x] = true;
                }
            }
        }
        marked.clear();
    }

    vector<int> phi(g.v_max);
    for (int v: g.vertices) {
        phi[v] = dist[v];
    }
    return phi;
}

/*
 * Returns the graph that is g with an added dummy vertex s
 * and edges of weight 0 connecting s to every vertex in G.
 */
Graph createGs(Graph &g) {
    int s = g.v_max;
    Graph Gs(g.v_max + 1, false);
    Gs.addVertices(g.vertices);
    Gs.addVertex(s);

    for (int u: g.vertices) {
        Gs.addEdges(u, g.adjacencyList[u], g.weights[u]);
    }

    vector<int> edges(g.n);
    vector<int> weights(g.n);
    for (int i = 0; i < g.vertices.size(); i++) {
            int v = g.vertices[i];
            edges[i] = v;
            weights[i] = 0;
        }
    Gs.addEdges(s, edges, weights);

    return Gs;
}

vector<int> getShortestPathTree(Graph &g, int s) {
    vector<bool> settled(g.v_max);
    custom_priority_queue<Node> pq(g.v_max);
    vector<int> dist(g.v_max, INT_MAX);
    vector<int> tree(g.v_max, -1);

    pq.emplace(s, 0);
    dist[s] = 0;

    while (!pq.empty()) {
        int u = pq.top().node;
        pq.pop();

        if (settled[u]) {
            continue;
        }

        settled[u] = true;
        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];
            if (!settled[v]) {
                int weight = g.weights[u][i];
                int newDistance = dist[u] + weight;
                if (newDistance < dist[v]) {
                    dist[v] = newDistance;
                    tree[v] = u;
                    pq.emplace(v, dist[v]);
                }
            }
        }
    }

    return tree;
}

vector<int> bellmanFord(Graph &g) {
    Timer::startTimer();
    vector<int> dist(g.v_max);
    for (int i = 0; i < g.v_max; i++) {
        dist[i] = INT_MAX;
    }
    dist[SRC] = 0;
    for (int i = 1; i < g.v_max; i++) {
        for (int u: g.vertices) {
            for (unsigned long j = 0; j < g.adjacencyList[u].size(); j++) {
                int v = g.adjacencyList[u][j];
                int weight = g.weights[u][j];

                if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                }
            }
        }
    }

    for (int u: g.vertices) {
        for (unsigned long j = 0; j < g.adjacencyList[u].size(); j++) {
            int v = g.adjacencyList[u][j];
            int weight = g.weights[u][j];

            if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                throw_with_nested("Negative cycle detected");
            }
        }
    }
    cout << "Bellman-Ford took: " << Timer::getDuration() << endl;
    return dist;
}

vector<int> getDistFromTree(Graph &g, vector<int> &tree) {
    vector<int> dist(g.v_max, INT_MAX);
    dist[SRC] = 0;

    updateDistFromTree(g, tree, dist, SRC);

    return dist;
}

void updateDistFromTree(Graph &g, vector<int> &tree, vector<int> &dist, int u) {
    for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
        int v = g.adjacencyList[u][i];
        int weight = g.weights[u][i];

        if (tree[v] == u) {
            dist[v] = dist[u] + weight;
            updateDistFromTree(g, tree, dist, v);
        }
    }
}