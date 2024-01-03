#include "LasVegas.h"

bool CHECK_RESULT = false;
const char *optionalFlag = "--source=";

int main(int argc, char *argv[]) {
//    string filename = "../graph/40000_25.txt";
    string filename = argv[1];

    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-c") == 0)
            CHECK_RESULT = true;
        if (strncmp(argv[i], optionalFlag, strlen(optionalFlag)) == 0) {
            SRC = stoi(argv[i] + strlen(optionalFlag));
        }
    }

    Random::Get().Seed();
    ifstream inputFile(filename);
    Graph g = readInput(inputFile);

    vector<int> LasVegasReal;
    try {
        LasVegasReal = lasVegas(g);
    } catch (const char *msg) {
        cout << msg << endl;
        return 0;
    }

    if (CHECK_RESULT) {
        vector<int> BellmanFord = bellmanFord(g);

        for (unsigned long i = 0; i < BellmanFord.size(); i++) {
            if (LasVegasReal[i] != BellmanFord[i]) {
                cout << "Failed: Bellman-Ford and Las Vegas are not equal" << endl;
                return 0;
            }
        }
        cout << "Success: Bellman-Ford and Las Vegas are equal" << endl;
    }
    return 0;
}
