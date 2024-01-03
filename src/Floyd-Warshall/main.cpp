#include <fstream>
#include <iostream>

#include "Floyd-Warshall.h"
#include "Timer.h"

using namespace std;

int main(int argc, char *argv[])
{
    string filename = argv[1];
    ifstream inputFile(filename);
    Graph g(inputFile);
    Timer::startTimer();
    g.floyd_warshall();
    cout << Timer::getDuration() << endl;
}