# from the project root:
#   build: cmake -S . -B build && cmake --build build
#   clean: cmake --build build --target clean

cmake_minimum_required(VERSION 4.3)
project(SuperFV LANGUAGES CXX)

add_executable(superfv src/superfv.cpp src/array.cpp)
