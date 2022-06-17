cmake_minimum_required(VERSION 3.9.1)

#project(2DSolverIncompressibleNS LANGUAGES CXX CUDA)
project(2DSolverIncompressibleNS LANGUAGES CUDA CXX)

find_package(CUDA  REQUIRED)
include_directories("${CUDA_INCLUDE_DIRS}")

# Setup CUDA project
include_directories(${CMAKE_CUDA_TOOLKIT_DIRECTORIES})
#include_directories(${CMAKE_CUDA_TOOLKIT_DIRECTORIES} ${CUDA_all_LIBRARIES})
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

add_executable(2dSolverNS
  main.cpp
  ${CMAKE_SOURCE_DIR}/src/compute_RHS.cpp
  ${CMAKE_SOURCE_DIR}/src/compute_curl.cpp
  ${CMAKE_SOURCE_DIR}/src/compute_projection_step.cpp
  ${CMAKE_SOURCE_DIR}/src/compute_RHS.cu
  ${CMAKE_SOURCE_DIR}/src/compute_projection_step.cu
  ${CMAKE_SOURCE_DIR}/src/time_advance_RK3.cu
  ${CMAKE_SOURCE_DIR}/src/helper_functions.cpp
  ${CMAKE_SOURCE_DIR}/src/time_advance_RK3.cpp
  ${CMAKE_SOURCE_DIR}/tests/computeErrorRHS.cpp
)
target_include_directories(2dSolverNS PRIVATE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/tests
) 

set_target_properties(2dSolverNS PROPERTIES CUDA_ARCHITECTURES "35;50;72")
#target_link_libraries(2dSolverNS ${CUDA_LIBRARIES})

# Find cufft
find_library(CUDA_CUFFT_LIBRARY cufft PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
target_link_libraries(2dSolverNS
    ${PROJECT_LINK_LIBS}
    ${CUDA_CUFFT_LIBRARY}
    ${CUDA_CUBLAS_LIBRARIES}
)