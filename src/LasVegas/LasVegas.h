//
// Created by Truong Giang Do on 28/11/2023.
//

#ifndef SSSP_NEW_LASVEGAS_H
#define SSSP_NEW_LASVEGAS_H

#include "LDD.h"

extern int SRC;
extern bool WITH_LDD;

vector<int> lasVegas(Graph &g);

Graph readInput(ifstream &inputFile);

void findReachable(Graph &g, int s, vector<bool> &reachable);

vector<int> bitScaling(Graph &g);

void verifyTree(Graph &g, vector<int> &tree, vector<int> &dist, int src);

bool containsCycles(Graph &g, vector<vector<bool>> &adjList, int src, vector<bool> &visited);

void getDistances(Graph &g, vector<int> &tree, vector<int> &dist, int curDis, int curVertex,
                  vector<bool> &visited);

vector<int> SPMain(Graph &g_in, int s);

Graph getScaledGraph(Graph &g_in, int scaleFactor);

bool invalidTree(Graph &g, int s, vector<int> &tree);

int roundPower2(int n);

double logBase2(int n);

vector<int> ScaleDown(Graph &g, int delta, int B);

vector<vector<int>> SPmainLDD(Graph &g, int diameter);

bool hasNegativeEdges(Graph &g, vector<int> &phi, int B);

Graph createModifiedGB(Graph &g, int B, bool nneg, set<vector<int>> &remEdges, vector<int> &phi);
Graph createModifiedGB2(Graph &g, int B, bool nneg, set<vector<int>> &remEdges, vector<int> &phi);

set<vector<int>> getEdgesBetweenSCCs(Graph &g, vector<int> &vertexToSCCMap);

vector<int> getVertexToSCCMap(vector<vector<int>> &SCCs, int numVertices);

vector<int> addPhi(vector<int> &phi_1, vector<int> &phi_2);

vector<int>
FixDAGEdges(Graph &g, vector<vector<int>> &SCCs, vector<int> &vertexToSCCMap, set<vector<int>> &edgesBetweenSCCs);

vector<vector<int>>
createSCCAdjList(vector<vector<int>> &SCCs, vector<int> &vertexToSCCMap, set<vector<int>> &edgesBetweenSCCs);

vector<int> topSort(int n, vector<vector<int>> &adjList);

void topSortUtil(int u, vector<bool> &visited, stack<int> &stack, vector<vector<int>> &adjList);

vector<int> ElimNeg(Graph &g);

Graph createGs(Graph &g);

vector<int> getShortestPathTree(Graph &g, int s);

vector<int> bellmanFord(Graph &g);

vector<int> getDistFromTree(Graph &g, vector<int> &tree);

void updateDistFromTree(Graph &g, vector<int> &tree, vector<int> &dist, int u);

#endif //SSSP_NEW_LASVEGAS_H
