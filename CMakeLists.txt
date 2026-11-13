# Build from the project root: cmake -S . -B build && cmake --build build
cmake_minimum_required(VERSION 4.3)
project(SuperFV LANGUAGES CXX)

add_executable(superfv src/superfv.cpp src/array.cpp)
