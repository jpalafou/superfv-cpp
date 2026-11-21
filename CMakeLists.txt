# from the project root:
#   build:  cmake -S . -B build && cmake --build build
#   run:    ./build/superfv
#   clean:  cmake --build build --target clean

cmake_minimum_required(VERSION 4.3)
project(SuperFV LANGUAGES CXX)

add_executable(superfv
    src/superfv.cpp
    src/array.cpp
    src/boundary_conditions.cpp
    src/hydro.cpp
    src/initial_condition_defs.cpp
    src/initial_conditions.cpp
    src/mesh.cpp
    src/riemann_solver.cpp
    src/stencils.cpp
    src/sweep.cpp
    src/update.cpp
)
set(SUPERFV_REAL_TYPE "double" CACHE STRING "Floating-point type")
target_compile_definitions(superfv PRIVATE
    SUPERFV_REAL_TYPE=${SUPERFV_REAL_TYPE}
)
target_compile_features(superfv PUBLIC cxx_std_20)
