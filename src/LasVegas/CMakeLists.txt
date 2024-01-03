cmake_minimum_required(VERSION 3.26)
project(SSSP_New)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

add_executable(SSSP_New main.cpp
        LDD.h
        Graph.h
        Randd.h
        LasVegas.h
        LDD.cpp
        LasVegas.cpp
        CustomPriorityQueue.h
        Timer.h
)
