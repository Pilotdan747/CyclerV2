cmake_minimum_required(VERSION 3.14)
project(Mars_Cycle)

set(CMAKE_CXX_COMPILER "g++-11")
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-fopenmp -L /usr/local/lib -lnlopt -lm -DARMA_DONT_USE_WRAPPER -framework accelerate")

# S1L1
add_executable(TestS1L1 src/TestS1L1.cpp src/cycler.cpp src/Lambert_Multi.cpp src/Lambert_Battin.cpp src/helperFuncs.cpp)

# main
add_executable(Cycle src/main.cpp src/cycler.cpp src/Lambert_Multi.cpp src/helperFuncs.cpp)

# optimize
add_executable(Optimize src/Optimize.cpp src/cycler.cpp src/Lambert_Multi.cpp src/helperFuncs.cpp)

# Derivative
add_executable(Derivative src/Derivative.cpp src/cycler.cpp src/Lambert_Multi.cpp src/helperFuncs.cpp)

# Derivative
add_executable(SingleCycle src/SingleCycler.cpp src/cycler.cpp src/Lambert_Multi.cpp src/helperFuncs.cpp)