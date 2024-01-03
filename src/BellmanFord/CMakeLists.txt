cmake_minimum_required(VERSION 3.26)
project(BellmanFord)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

add_executable(BellmanFord main.cpp
        Bellman-Ford.h
        Timer.h
)
