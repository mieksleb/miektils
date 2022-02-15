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


void FENE_force(Tensor<double,2> &force, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                double &d21, double &d22, double kappa, double r0, bool left, bool right);

void KRATKY_force( Tensor<double,2> &force, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                   double &d21, double &d22, double A);


class forces {

};


#endif //MIEKTEP_FORCES_H
