cmake_minimum_required(VERSION 3.9.1)

project(testcomputeRHS)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

add_executable(testcomputeRHS
  test_computeRHS.cpp
  computeErrorRHS.cpp
  ${CMAKE_SOURCE_DIR}/../src/compute_RHS.cpp
)
target_include_directories(testcomputeRHS PRIVATE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/src
) 
