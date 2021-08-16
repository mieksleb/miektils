cmake_minimum_required(VERSION 3.15)
project(miekTEP)

set(CMAKE_CXX_STANDARD 14)

set(EIGEN_UNSUPP_DIR "~/c++/eigen/")
include_directories(${EIGEN_UNSUPP_DIR})

add_executable(miektep main.cpp tensor_tools.cpp tensor_tools.h forces.cpp forces.h model.h)

