//
// Created by Michael Selby on 02/08/2021.
//

#ifndef MIEKTEP_TENSOR_TOOLS_H
#define MIEKTEP_TENSOR_TOOLS_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace std;
using namespace Eigen;

void print_rank_4_tensor(Tensor<double,4> tens);
void print_rank_3_tensor(Tensor<double,3> tens);
void print_rank_2_tensor(Tensor<double,2> tens);
void print_rank_1_tensor(Tensor<double,1> tens);
Tensor<double,1> get_slice_rank2(Tensor<double,2> input, int i);
Tensor<double,1> get_slice_rank3(Tensor<double,3> input, int i, int j);
double dot(Tensor<double,1> vec1,Tensor<double,1> vec2);
double norm2(Tensor<double,1> vec);
Tensor<double,2> get_noise(int np);

#endif //MIEKTEP_TENSOR_TOOLS_H
