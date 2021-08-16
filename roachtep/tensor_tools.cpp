//
// Created by Michael Selby on 02/08/2021.
//
// Tensor related functions
//

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;
#include "tensor_tools.h"

void print_rank_4_tensor(Tensor<double,4> tens) {
    Eigen::array<Index, 4> dims = tens.dimensions();
    for(int i = 0; i < dims[0]; i++) {
        for(int j = 0; j < dims[1]; j++) {
            for (int alpha = 0; alpha < dims[2]; alpha++) {
                for (int beta = 0; beta < dims[3]; beta++) {
                    cout << tens(i,j,alpha,beta) << endl;
                }
            }
        }
    }
}

void print_rank_3_tensor(Tensor<double,3> tens) {
    Eigen::array<Index, 3> dims = tens.dimensions();
    for(int i = 0; i < dims[0]; i++) {
        for(int j = 0; j < dims[1]; j++) {
            for (int alpha = 0; alpha < dims[2]; alpha++) {
                    cout << tens(i,j,alpha) << endl;
            }
        }
    }
}


void print_rank_2_tensor(Tensor<double,2> tens) {
    // Prints rank 2 tensor
    Eigen::array<Index, 2> dims = tens.dimensions();
    for(int i = 0; i < dims[0]; i++) {
        for (int alpha = 0; alpha < dims[1]; alpha++) {
            cout << tens(i,alpha) << " ";
        }
        cout << "\n";
    }
}


Tensor<double,1> get_slice(Tensor<double,2> input, int i) {
    // This function slices a the ith column vector from a rank 2 tensor
    // Input: Rank 2 Tensor (n x m matrix)
    // Output: Rank 1 tensor (n dimensional vector)
    Eigen::array<Eigen::Index, 2> dims = input.dimensions();
    Eigen::array<long, 2> offsets = {i, 0};
    Eigen::array<long, 2> extents = {1, dims[1]};
    Tensor<double, 1> slice = input.slice(offsets, extents).reshape(Eigen::array<long, 1>{3});
    return slice;
}


double norm2(Tensor<double,1> vec) {
    // Returns the squared norm of a vector (rank 1 tensor)
    Eigen::array<Eigen::Index, 1> dims = vec.dimensions();
    double val = 0;
    for (int i = 0; i < dims[0]; i++) {
        val += vec(i) * vec(i);
    }
    return val;
}

Tensor<double,2> get_noise(int np) {
    // this function returns a tensor who's rows are gaussian random variables with zero mean and unit variance
    unsigned seed;
    Tensor<double,2> noise(np,3);
    seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    for (int i = 0; i < np; i++) {
        for (int alpha = 0; alpha < 3; alpha++) {
            normal_distribution<double> distribution(0.0, 1.0);
            noise(i, alpha) = distribution(generator);
        }
    }
    return noise;
}