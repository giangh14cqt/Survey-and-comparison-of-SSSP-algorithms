//
// Created by Truong Giang Do on 05/12/2023.
//

#ifndef SSSP_NEW_TIMER_H
#define SSSP_NEW_TIMER_H

#include <chrono>

class Timer {
public:
    static void startTimer() {
        getInstance().startTime = std::chrono::high_resolution_clock::now();
    }

    static double getDuration() {
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - getInstance().startTime);
        return duration.count() * 1e-3; // Convert to milliseconds
    }

    static void startDebugTimer() {
        getInstance().debugTime -= getDuration();
    }

    static void stopDebugTimer() {
        getInstance().debugTime += getDuration();
    }

    static void resetDebugTimer() {
        getInstance().debugTime = 0;
    }

    static double getDebugDuration() {
        return getInstance().debugTime;
    }

private:
    Timer() {} // Private constructor to ensure only one instance can be created

    static Timer& getInstance() {
        static Timer instance; // Guaranteed to be destroyed and instantiated on first use.
        return instance;
    }

    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    double debugTime;
};

#endif //SSSP_NEW_TIMER_H
