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
Tensor<double,1> get_slice(Tensor<double,2> input, int i);
double norm2(Tensor<double,1> vec);
Tensor<double,2> get_noise(int np);

#endif //MIEKTEP_TENSOR_TOOLS_H
