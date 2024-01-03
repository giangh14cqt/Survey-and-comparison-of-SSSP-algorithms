cmake_minimum_required(VERSION 3.26)
project(Dijkstra)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

add_executable(Dijkstra main.cpp
        Dijkstra.h
        Timer.h
)
