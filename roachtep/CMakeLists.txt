cmake_minimum_required(VERSION 3.15)
project(miekTEP)

set(CMAKE_CXX_STANDARD 14)

set(EIGEN_UNSUPP_DIR "~/c++/eigen/")
include_directories(${EIGEN_UNSUPP_DIR})

add_executable(roachtep main.cpp forces.cpp forces.h)

