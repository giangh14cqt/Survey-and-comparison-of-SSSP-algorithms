//
// Created by Truong Giang Do on 28/11/2023.
//

#ifndef SSSP_NEW_RANDD_H
#define SSSP_NEW_RANDD_H

#include "Graph.h"

class Random {
public:
private:
    // Make the default constructor private.
    Random() {}

public:
    // Delete the copy constructor function.
    Random(const Random &) = delete;

    // Delete the overloading of assignment operator
    Random &operator=(const Random &) = delete;

    static Random &Get() {
        static Random inst;
        return inst;
    }

    vector<int> randomSubset(vector<int> &set, int k) {
        vector<int> subset(k);
        for (int i = 0; i < k; i++) {
            int index = GenInt(0, set.size() - 1);
            subset[i] = set[index];
        }
        return subset;
    }

    // Seed the random number generator.
    static void Seed() {
        srand(time(nullptr));
    }

    static int GenInt() {
        return rand();
    }

    static int GenInt(int min, int max) {
        return rand() % (max - min + 1) + min;
    }

    static double GenDouble() {
        return (double) rand() / RAND_MAX;
    }
};

#endif //SSSP_NEW_RANDD_H
