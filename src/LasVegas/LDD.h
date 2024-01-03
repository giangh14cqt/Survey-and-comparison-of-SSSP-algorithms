//
// Created by Truong Giang Do on 28/11/2023.
//

#ifndef SSSP_NEW_LDD_H
#define SSSP_NEW_LDD_H

#include "Randd.h"
#include "CustomPriorityQueue.h"

extern int gn_global;

vector<vector<int>> preLDD(Graph &g, int d);

bool hasLargeDiameter(Graph &g, int s, int diameter);

vector<vector<int>> LDD(Graph &g, int d);

vector<vector<int>> LDDRework(Graph &g, int d);

vector<vector<int>> boundary(Graph &g, Graph &g_rev, vector<int> &ball_in, vector<int> &ball_out);

vector<vector<int>> revEdges(vector<vector<int>> &edges);

double calculateGeoProb(int n, int r);

vector<vector<int>> RandomTrim(Graph &g, Graph &g_rev, int s, int d);

Graph getSubGraph(Graph &g, vector<int> &ball, bool setMinus);

vector<int> vertexUnion(vector<int> &set1, vector<int> &set2);

void addVerticesToSet(set<int> &set, vector<int> &vertices);

vector<vector<int>> edgeUnion(vector<vector<int>> &set1,
                              vector<vector<int>> &set2,
                              vector<vector<int>> &set3);

int diffVertex(vector<int> &set1, vector<int> &set2, int v_max);

vector<int> CoreOrLayerRange(Graph &g, Graph &g_rev, int s, int d);

Graph createGRev(Graph &g);

vector<int> oneIterationLayerRange(Graph &g, custom_priority_queue<Node> &pq, vector<bool> &settled,
                                   int numSettled, vector<vector<int>> &farthestDistancesSeen, double constant,
                                   vector<int> &dist, int d);

bool sameCanonicalRange(vector<vector<int>> &farthestDistancesSeen, double constant);

vector<vector<int>> layer(Graph &g, vector<int> &ball);

vector<int> volume(Graph &g, int s, int r);

vector<int> volume2(Graph &g, int s, int r);

vector<int> Dijkstra(Graph &g, int s);

void updateNeighbors(Graph &g, int u, vector<bool> &settled, custom_priority_queue<Node> &pq, vector<int> &dist, int d);

void updateNeighbors2(Graph &g, int u, vector<int> &VerToN, vector<bool> &settled, custom_priority_queue<Node> &pq, vector<int> &dist, int d);

void init(Graph &g, custom_priority_queue<Node> &pq, vector<int> &dist, int s);

#endif //SSSP_NEW_LDD_H
