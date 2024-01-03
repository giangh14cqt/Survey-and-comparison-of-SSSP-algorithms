//
// Created by Truong Giang Do on 28/11/2023.
//
#include "LDD.h"

const int LDD_BASE_CASE = 10;
const int constant_c = 1;
int gn_global = 0;


vector<vector<int>> preLDD(Graph &g, int d) {
    return LDD(g, d);
}

// OUTPUT: A set of edges e_sep with the following guarantees:
// – each SCC of G\e_sep has weak diameter at most D; that is, if u,v are in the
// same SCC,
// then dist_G(u, v) ≤ D and dist_G(v, u) ≤ D.
// – For every e ∈ E, Pr[e ∈ e_sep] = O(w(e)·(logn)^2/D +n^−10). These
// probabilities are not
// guaranteed to be independent.
// Each vector<int> in the output ArrayList has size two and represents an edge
// (int[0], int[1])
vector<vector<int>> LDD(Graph &g, int d) {
//        if (g.n <= max(1, LDD_BASE_CASE)) {
//            return {};
//        }

    Graph g_rev = createGRev(g);

    int s = -1;
    bool foundGoodS = false;
    for (int v: g.vertices) {
        if (!g.adjacencyList[v].empty() || !g_rev.adjacencyList[v].empty()) {
            s = v;
            foundGoodS = true;
            break;
        }
    }

    if (!foundGoodS)
        return {};

    vector<int> condAndi_max = CoreOrLayerRange(g, g_rev, s, d);

    if (condAndi_max[0] == 1)
        return RandomTrim(g, g_rev, s, d);

    int r = (int) ceil(d / (3.0 * log(g.n)));
    int i_min = condAndi_max[1] - r;

    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<int> distr(calculateGeoProb(g.n, r));

    int i_tilda = distr(gen);
    int i_rnd = i_min + min(i_tilda, r);

    if (condAndi_max[0] == 2) {
        vector<int> ball = volume(g, s, i_rnd);
        Graph subGraph = getSubGraph(g, ball, false);
        Graph minusSubGraph = getSubGraph(g, ball, true);

        vector<vector<int>> layer_g = layer(g, ball);
        vector<vector<int>> preLDD_subGraph = preLDD(subGraph, d);
        vector<vector<int>> preLDD_minusSubGraph = preLDD(minusSubGraph, d);
        return edgeUnion(layer_g, preLDD_subGraph, preLDD_minusSubGraph);
    }

    if (condAndi_max[0] == 3) {
        vector<int> ball = volume(g_rev, s, i_rnd);
        Graph subGraph = getSubGraph(g_rev, ball, false);
        Graph minusSubGraph = getSubGraph(g_rev, ball, true);

        vector<vector<int>> layer_g_rev = layer(g_rev, ball);
        vector<vector<int>> preLDD_subGraph = preLDD(subGraph, d);
        vector<vector<int>> preLDD_minusSubGraph = preLDD(minusSubGraph, d);
        vector<vector<int>> rev = edgeUnion(layer_g_rev, preLDD_subGraph, preLDD_minusSubGraph);
        return revEdges(rev);
    }

    throw_with_nested("LowDiamDecomposition failed.");
}


vector<vector<int>> LDDRework(Graph &g0, int d) {
    Graph g(g0, false);
    Graph g_rev = createGRev(g);
    set<vector<int>> e_sep;

    int k = constant_c * int(log2(gn_global));
    vector<int> randomVetices = Random::Get().randomSubset(g.vertices, k);

    vector<vector<int>> ball_g_in(randomVetices.size());
    for (int i = 0; i < randomVetices.size(); ++i) {
        ball_g_in[i] = volume(g, randomVetices[i], d / 4);
        ball_g_in[i] = vertexUnion(ball_g_in[i], randomVetices);
    }
    vector<vector<int>> ball_g_out(randomVetices.size());
    for (int i = 0; i < randomVetices.size(); ++i) {
        ball_g_out[i] = volume(g_rev, randomVetices[i], d / 4);
        ball_g_out[i] = vertexUnion(ball_g_out[i], randomVetices);
    }

    queue<int> lightNodeV;
    for (int i = 0; i < ball_g_in.size(); ++i)
        if (ball_g_in[i].size() < 0.5 * g.n && ball_g_out[i].size() < 0.5 * g.n)
            lightNodeV.push(i);

    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<int> distr(min(1.0, 80.0 * log2(g.n) / d));
    while (!lightNodeV.empty()) {
        int v = g0.vertices[lightNodeV.front()];
        lightNodeV.pop();

        if (!g.containsVertex[v])
            continue;

        int rv = distr(gen);
        vector<int> ball_in = volume(g, v, rv);
        vector<int> ball_out = volume(g_rev, v, rv);
        vector<vector<int>> boun = boundary(g, g_rev, ball_in, ball_out);
        if (rv > d / 4 || ball_in.size() > 0.7 * g.n || ball_out.size() > 0.7 * g.n) {
//            cout << "Rare case" << endl;
            return g.adjacencyList;
        }
        vector<int> mergeBall = vertexUnion(ball_in, ball_out);
        set<int> mergeBallSet(mergeBall.begin(), mergeBall.end());
        mergeBall = vector<int>(mergeBallSet.begin(), mergeBallSet.end());
        Graph subGraph = getSubGraph(g, mergeBall, false);
        vector<vector<int>> e_recurse = LDDRework(subGraph, d);
        e_sep.insert(boun.begin(), boun.end());
        e_sep.insert(e_recurse.begin(), e_recurse.end());
        g.removeVertices(mergeBall);
    }

    if (!g.vertices.empty()) {
        int arbitraryVertex = g0.vertices[0];
        vector<int> ball_in = volume(g0, arbitraryVertex, d / 2);
        vector<int> ball_out = volume(g_rev, arbitraryVertex, d / 2);
        set<int> vertices(g.vertices.begin(), g.vertices.end());
        set<int> ball_in_set(ball_in.begin(), ball_in.end());
        set<int> ball_out_set(ball_out.begin(), ball_out.end());
        if (!includes(vertices.begin(), vertices.end(), ball_in_set.begin(), ball_in_set.end()) ||
            !includes(vertices.begin(), vertices.end(), ball_out_set.begin(), ball_out_set.end())) {
//            cout << "Rare case 2" << endl;
            return g.adjacencyList;
        }
    }

    return vector<vector<int>>(e_sep.begin(), e_sep.end());
}

vector<vector<int>> boundary(Graph &g, Graph &g_rev, vector<int> &ball_in, vector<int> &ball_out) {
    set<vector<int>> edges;
    for (int u: ball_out) {
        for (int v: g.adjacencyList[u]) {
            if (find(ball_out.begin(), ball_out.end(), v) == ball_out.end())
                edges.insert({u, v});
        }
    }
    for (int u: ball_in) {
        for (int v: g_rev.adjacencyList[u]) {
            if (find(ball_in.begin(), ball_in.end(), v) == ball_in.end())
                edges.insert({v, u});
        }
    }

    return vector<vector<int>>(edges.begin(), edges.end());
}

vector<vector<int>> revEdges(vector<vector<int>> &edges) {
    vector<vector<int>> revEdgeSet(edges.size());
    for (unsigned long i = 0; i < edges.size(); i++) {
        revEdgeSet[i] = {edges[i][1], edges[i][0]};
    }
    return revEdgeSet;
}

double calculateGeoProb(int n, int r) {
    double prob = pow(log(n), 2) / r;
    return min(1.0, prob);
}

vector<vector<int>> RandomTrim(Graph &g, Graph &g_rev, int s, int d) {
    vector<vector<int>> e_sep;

    if (INT_MAX / 4 <= d)
        return e_sep;

    vector<int> dist = Dijkstra(g, s);
    vector<int> dist_rev = Dijkstra(g_rev, s);

    int numDistOver4D = 0;
    int over4DVert = -1;

    vector<int> v_far;
    for (int v: g.vertices) {
        if (max(dist[v], dist_rev[v]) > 2 * d) {
            v_far.push_back(v);
            if (max(dist[v], dist_rev[v]) > 4 * d) {
                numDistOver4D++;
                over4DVert = v;
            }
        }
    }

    // base case
    if (numDistOver4D <= 1) {
        if (numDistOver4D == 1) {
            for (int v: g.adjacencyList[over4DVert]) {
                e_sep.push_back({over4DVert, v});
            }
        }
        return e_sep;
    }

    vector<int> m; // marked vertices
    int i_max = d;
    int r = (int) ceil(d / (3.0 * log(g.n)));
    int i_min = i_max - r;

    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<int> distr(calculateGeoProb(g.n, r));

    int v = diffVertex(v_far, m, g.v_max);
    while (v != -1) {
        int i_rnd = i_min + min(distr(gen), r);

        if (dist_rev[v] > 2 * d) {
            Graph gVMinusM = getSubGraph(g, m, true);
            vector<int> ball = volume(gVMinusM, v, i_rnd);
            Graph GVMinusMSubGraph = getSubGraph(gVMinusM, ball, false);
//                e_sep = edgeUnion(e_sep, layer(gVMinusM, ball), preLDD(GVMinusMSubGraph, d));
            vector<vector<int>> layer_gVMinusM = layer(gVMinusM, ball);
            vector<vector<int>> preLDD_GVMinusMSubGraph = preLDD(GVMinusMSubGraph, d);
            e_sep = edgeUnion(e_sep, layer_gVMinusM, preLDD_GVMinusMSubGraph);
            m = vertexUnion(m, ball);
        } else if (dist[v] > 2 * d) {
            Graph gVMinusM_rev = getSubGraph(g_rev, m, true);
            vector<int> ball_rev = volume(gVMinusM_rev, v, i_rnd);
            Graph GVMinusMSubGraph_rev = getSubGraph(gVMinusM_rev, ball_rev, false);
//                e_sep = edgeUnion(e_sep, revEdges(layer(gVMinusM_rev, ball_rev)),
//                                  revEdges(preLDD(GVMinusMSubGraph_rev, d)));
            vector<vector<int>> layer_gVminusM_rev = layer(gVMinusM_rev, ball_rev);
            vector<vector<int>> preLDD_GVMinusMSubGraph_rev = preLDD(GVMinusMSubGraph_rev, d);
            layer_gVminusM_rev = revEdges(layer_gVminusM_rev);
            preLDD_GVMinusMSubGraph_rev = revEdges(preLDD_GVMinusMSubGraph_rev);
            e_sep = edgeUnion(e_sep, layer_gVminusM_rev, preLDD_GVMinusMSubGraph_rev);
            m = vertexUnion(m, ball_rev);
        } else {
            throw_with_nested("RandomTrim failed.");
        }

        v = diffVertex(v_far, m, g.v_max);
    }
    return e_sep;
}

// returns the subgraph of g containing only the vertices in ball
// if setMinus is true, the function returns the subgraph of g containing only
// the vertices outside of the ball
Graph getSubGraph(Graph &g, vector<int> &ball, bool setMinus) {
    vector<bool> contains(g.v_max, false);
    for (int i: ball)
        if (g.containsVertex[i])
            contains[i] = true;

    vector<int> vert;
    for (int v: g.vertices)
        if (!setMinus && contains[v])
            vert.push_back(v);
        else if (setMinus && !contains[v])
            vert.push_back(v);

    Graph subGraph(g.v_max, false);
    subGraph.addVertices(vert);

    for (int u: vert) {
        vector<int> edges;
        vector<int> weights;

        for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            if (subGraph.containsVertex[v]) {
                edges.push_back(v);
                weights.push_back(g.weights[u][i]);
            }
        }

        subGraph.addEdges(u, edges, weights);
    }
    return subGraph;
}

// returns the union of two vertex sets
vector<int> vertexUnion(vector<int> &set1, vector<int> &set2) {
    set<int> set;
    set.insert(set1.begin(), set1.end());
    set.insert(set2.begin(), set2.end());

    return vector<int>(set.begin(), set.end());
}

vector<vector<int>> edgeUnion(vector<vector<int>> &set1,
                              vector<vector<int>> &set2,
                              vector<vector<int>> &set3) {
    set<vector<int>> set(set1.begin(), set1.end());
    set.insert(set2.begin(), set2.end());
    set.insert(set3.begin(), set3.end());
    return {set.begin(), set.end()};
}

int diffVertex(vector<int> &set1, vector<int> &set2, int v_max) {
    vector<bool> contains(v_max, false);
    for (int i: set2)
        contains[i] = true;

    for (int i: set1)
        if (!contains[i])
            return i;

    return -1;
}

// OUTPUT: a pair (Condition,i) where Condition ∈ {1,2,3} and i ≤ D is a
// non-negative integer such that
// – if Condition = 1 then n_G(s,i) > 2n and n_G_rev(s,i) > 2n,
// - if Condition = 2 then n_G(s, i) ≤ 2n and Vol_G(s, i) and Vol_G(s, i −
// ⌈D/(3lgn)⌉) are the same canonical range,
// – if Condition = 3 then n_G_rev(s, i) ≤ 2n and Vol_G_rev(s, i) and
// Vol_G_rev(s, i − ⌈D/(3lgn)⌉) are in the same canonical range.
// If Condition ∈ {2, 3} then i ≥ D/(3lgn).
// Runs LayerRange on G and G_rev in parallel.
vector<int> CoreOrLayerRange(Graph &g, Graph &g_rev, int s, int d) {
    vector<vector<int>> farthestDistancesSeen;
    vector<vector<int>> farthestDistancesSeen_rev;
    double constant = d / (3.0 * log(g.n));
    vector<bool> settled(g.v_max, false);
    vector<bool> settled_rev(g_rev.v_max, false);
    int numSettled = 0;
    int numSettled_rev = 0;
    custom_priority_queue<Node> pq(g.v_max);
    custom_priority_queue<Node> pq_rev(g.v_max);
    vector<int> dist(g.v_max, INT_MAX);
    vector<int> dist_rev(g.v_max, INT_MAX);
    init(g, pq, dist, s);
    init(g_rev, pq_rev, dist_rev, s);
    bool finished = false;
    bool finished_rev = false;
    int j = -1;
    int j_rev = -1;

    while (true) {
        if (numSettled == g.n)
            finished = true;

        if (numSettled_rev == g.n)
            finished_rev = true;

        if (finished && finished_rev) {
            // case 1
            if (j == -1 || j_rev == -1)
                return {2, (int) ceil(constant)};

            return {1, max(j, j_rev)};
        }

        if (!finished) {
            vector<int> result = oneIterationLayerRange(g, pq, settled, numSettled, farthestDistancesSeen, constant,
                                                        dist, d);
            if (!result.empty()) {
                if (result[0] == 1) {
                    j = result[1];
                    finished = true;
                } else if (result[0] == 2) {
                    // case 2
                    return {2, result[1]};
                }
            }
            numSettled++;
        }

        if (!finished_rev) {
            vector<int> result_rev = oneIterationLayerRange(g_rev, pq_rev, settled_rev, numSettled_rev,
                                                            farthestDistancesSeen_rev, constant, dist_rev, d);
            if (!result_rev.empty()) {
                if (result_rev[0] == 1) {
                    j_rev = result_rev[1];
                    finished_rev = true;
                } else if (result_rev[0] == 2) {
                    // case 3
                    return {3, result_rev[1]};
                }
            }
            numSettled_rev++;
        }
    }
}

// returns a copy of g but with edges reversed
Graph createGRev(Graph &g) {
    Graph g_rev(g.v_max, false);
    g_rev.addVertices(g.vertices);

    vector<vector<int>> edges(g.v_max);
    vector<vector<int>> weights(g.v_max);

    for (int v: g.vertices) {
        for (unsigned long i = 0; i < g.adjacencyList[v].size(); i++) {
            edges[g.adjacencyList[v][i]].push_back(v);
            weights[g.adjacencyList[v][i]].push_back(g.weights[v][i]);
        }
    }

    for (int i = 0; i < g.v_max; i++) {
        g_rev.addEdges(i, edges[i], weights[i]);
    }
    return g_rev;
}

vector<int> oneIterationLayerRange(Graph &g, custom_priority_queue<Node> &pq, vector<bool> &settled,
                                   int numSettled, vector<vector<int>> &farthestDistancesSeen, double constant,
                                   vector<int> &dist, int d) {
    if (pq.empty()) {
        /*
         * Nothing left to search.
         * g is disconnected, since n_G(s,i) <= 2n/3.
         * Will never be the case that i ≥ D/(3lgn).
         * Return i_big = i + ceil(D/(3lgn)).
         * Guaranteed that i_big will satisfy i_big >= D/(3lgn) and the canonical
         * ranges.
         */
        int farthestDistanceSeen = farthestDistancesSeen[farthestDistancesSeen.size() - 1][0];
        int i_big = min(d, farthestDistanceSeen + (int) ceil(constant));
        return {2, i_big};
    }
    int u = pq.top().node;
    pq.pop();

    if (settled[u])
        return {};

    settled[u] = true;

    if (farthestDistancesSeen.empty() || dist[u] > farthestDistancesSeen[farthestDistancesSeen.size() - 1][0]) {
        farthestDistancesSeen.push_back({dist[u], numSettled + 1});
    }

    int farthestDistanceSeen = farthestDistancesSeen[farthestDistancesSeen.size() - 1][0];

    // case 1
    if (numSettled + 1 > 2.0 * g.n / 3.0)
        return {1, farthestDistanceSeen};

    // case 2
    if ((farthestDistanceSeen >= constant) && sameCanonicalRange(farthestDistancesSeen, constant))
        return {2, farthestDistanceSeen};

    updateNeighbors(g, u, settled, pq, dist, d);

    return {};
}

// Checks whether Vol_G(s, i - ceil[D/(3logn)]) and Vol_G(s, i) are in the same
// canonical range.
// Two numbers are in the same canonical range if they lie in the same half-open
// interval
// [2^j, 2^{j+1}), where j is a non-negative integer.
bool sameCanonicalRange(vector<vector<int>> &farthestDistancesSeen, double constant) {
    int i = farthestDistancesSeen[farthestDistancesSeen.size() - 1][0];
    int vol1 = farthestDistancesSeen[farthestDistancesSeen.size() - 1][1];

    for (unsigned long j = farthestDistancesSeen.size() - 2; j >= 0; j--)
        if (farthestDistancesSeen[j][0] <= i - ceil(constant)) {
            int vol2 = farthestDistancesSeen[j][1];
            if (floor(log(vol1) / log(2)) == floor(log(vol2) / log(2)))
                return true;
            break;
        }
    return false;
}

// {(u,v) in E_H | u in V_H(s,r) and v not in V_H(s,r)}

vector<vector<int>> layer(Graph &g, vector<int> &ball) {
    vector<bool> contains(g.v_max, false);
    for (int i: ball)
        contains[i] = true;
    vector<vector<int>> edges;
    for (int u: ball) {
        for (int v: g.adjacencyList[u]) {
            if (!contains[v])
                edges.push_back({u, v});
        }
    }
    return edges;
}

// returns all the vertices in g within a distance of r from source vertex s
// using Dijkstra's
vector<int> volume(Graph &g, int s, int r) {
    vector<int> output;
    vector<bool> settled(g.v_max, false);
    custom_priority_queue<Node> pq(g.n);
    vector<int> dist(g.v_max, INT_MAX);
    init(g, pq, dist, s);
    while (!pq.empty()) {
        int u = pq.top().node;
        pq.pop();

        if (settled[u] || dist[u] > r)
            continue;

        output.push_back(u);
        settled[u] = true;

        updateNeighbors(g, u, settled, pq, dist, r);
    }
    return output;
}

// returns all the vertices in g within a distance of r from source vertex s
// using Dijkstra's
vector<int> volume2(Graph &g, int s, int r) {
    vector<int> output;
    vector<bool> settled(g.n, false);
    custom_priority_queue<Node> pq(g.n);
    vector<int> VerToN(g.v_max, -1);
    for (int i = 0; i < g.vertices.size(); ++i) {
        VerToN[g.vertices[i]] = i;
    }
    vector<int> dist(g.n, INT_MAX);
    s = VerToN[s];
    init(g, pq, dist, s);
    while (!pq.empty()) {
        int u = pq.top().node;
        pq.pop();

        if (settled[u] || dist[u] > r)
            continue;

        output.push_back(g.vertices[u]);
        settled[u] = true;

        updateNeighbors2(g, u, VerToN, settled, pq, dist, r);
    }
    return output;
}

void updateNeighbors2(Graph &g, int u,
                      vector<int> &VerToN, vector<bool> &settled,
                      custom_priority_queue<Node> &pq, vector<int> &dist,
                      int d) {
    int u_real = g.vertices[u];
    for (unsigned long i = 0; i < g.adjacencyList[u_real].size(); i++) {
        int v_real = g.adjacencyList[u_real][i];
        int v = VerToN[v_real];
        if (!settled[v]) {
            dist[v] = min(dist[v], dist[u] + g.weights[u_real][i]);

            // only want to process nodes within a distance of d from the source
            if (dist[v] <= d)
                pq.push(Node(v, dist[v]));
        }
    }
}

vector<int> Dijkstra(Graph &g, int s) {
    vector<bool> settled(g.v_max, false);
    custom_priority_queue<Node> pq(g.v_max);
    vector<int> dist(g.v_max, INT_MAX);
    pq.push(Node(s, 0));
    dist[s] = 0;

    while (!pq.empty()) {
        int u = pq.top().node;
        pq.pop();

        if (settled[u])
            continue;

        settled[u] = true;

        updateNeighbors(g, u, settled, pq, dist, INT_MAX);
    }
    return dist;
}


void init(Graph &g, custom_priority_queue<Node> &pq, vector<int> &dist, int s) {
    for (int i = 0; i < g.v_max; i++) {
        dist[i] = INT_MAX;
    }
    pq.push(Node(s, 0));
    dist[s] = 0;
}

void updateNeighbors(Graph &g, int u,
                     vector<bool> &settled,
                     custom_priority_queue<Node> &pq,
                     vector<int> &dist, int d) {
    for (unsigned long i = 0; i < g.adjacencyList[u].size(); i++) {
        int v = g.adjacencyList[u][i];
        if (!settled[v]) {
            dist[v] = min(dist[v], dist[u] + g.weights[u][i]);

            // only want to process nodes within a distance of d from the source
            if (dist[v] <= d)
                pq.push(Node(v, dist[v]));
        }
    }
}