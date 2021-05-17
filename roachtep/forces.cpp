//
// Created by Michael Selby on 04/12/2020.
// This is where all external force functions are stored, these include:
// - FENE spring
// - Bending force
// - Dihedral force
//
//

#include "forces.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <random>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


// The FENE force is simply a function of the difference between two particles i and j.
// diffvec is the rank 2 tensor (np x np) which contains all values of the  distance of separation, i.e. difftens(i,j) = |ri-rj|.
// diff2 is the rank 3 (np x np x 3) tensor which contains all values of the separation, i.e. diffvectens(i,j,alpha) = ri,alpha - rj,alpha.
Tensor<double,1> FENE_force(int np, Tensor<double,1> diffvec,Tensor<double,1> diff2, double d, double kappa, double boltzmann, double r0) {
    double weekscut = 1.12246204831 * d; // 1.122... is 2^(1/6), this is the Weeks-Chandler cutoff distance.
    double diff = pow(diff2(0),0.5);
    double rat0 = diff/r0; // this is (rij)^2/R0, which is the FENE spring part of the force
    if (diff < weekscut) {
        double ratio6 = pow(d / diff2(0),3); // this is (d/rij)^6.
        FENE_force = diffvec / diff * ((24 * boltzmann / diff * ratio6 * (1 - 2 * ratio6)) + kappa*r0*rat0/(1-pow(rat0,2)));

    }
    else {
        FENE_force = (diffvec / diff)*kappa*r0*rat0/(1-pow(rat0,2));
    }
}