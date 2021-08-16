//
// Created by Michael Selby on 04/12/2020.
//

#ifndef MIEKTEP_FORCES_H
#define MIEKTEP_FORCES_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "tensor_tools.h"


using namespace std;
using namespace Eigen;


Tensor<double,2> FENE_force(int np, Tensor<double,3> &diffvectens,Tensor<double,2> &diff2vals, double d, double kappa, double boltzmann, double r0);

class forces {

};


#endif //MIEKTEP_FORCES_H
