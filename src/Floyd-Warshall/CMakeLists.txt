cmake_minimum_required(VERSION 3.26)
project(FloydWarshall)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

add_executable(FloydWarshall main.cpp
        Floyd-Warshall.h
        Timer.h
)
