cmake_minimum_required(VERSION 3.24)
project(solve4)

set(CMAKE_CXX_STANDARD 17)

add_executable(solve4 main.cpp Solve.h Solve.cpp)
add_subdirectory(tests)