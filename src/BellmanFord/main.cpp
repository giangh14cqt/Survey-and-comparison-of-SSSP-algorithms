#include <fstream>
#include <iostream>

#include "Bellman-Ford.h"
#include "Timer.h"

using namespace std;

int main(int argc, char *argv[])
{
    string filename = argv[1];
    int SRC = 0;
    if (argc > 2)
        SRC = stoi(argv[3]);
    ifstream inputFile(filename);
    Graph g(inputFile);
    Timer::startTimer();
    g.bellman_ford(SRC);
    cout << Timer::getDuration() << endl;
}