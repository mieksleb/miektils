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
#include "tensor_tools.h"


using namespace std;
using namespace Eigen;


// The FENE force is simply a function of the difference between two particles i and j.
// diffvec is the rank 2 tensor (np x np) which contains all values of the  distance of separation, i.e. difftens(i,j) = |ri-rj|.
// diff2 is the rank 3 (np x np x 3) tensor which contains all values of the separation, i.e. diffvectens(i,j,alpha) = ri,alpha - rj,alpha.


/*Tensor<double,1> FENE_force(int np, Tensor<double,1> diffvec,Tensor<double,1> diff2, double d, double kappa, double boltzmann, double r0) {
    double weekscut = 1.12246204831 * d; // 1.122... is 2^(1/6), this is the Weeks-Chandler cutoff distance.
    double diff = pow(diff2(0),0.5);
    double rat0 = diff/r0; // this is (rij)^2/R0, which is the FENE spring part of the force
    if (diff < weekscut) {
        double ratio6 = pow(d / diff2(0),3); // this is (d/rij)^6.
        FENE_force = diffvec / (diff * ((24 * boltzmann / diff * ratio6 * (1 - 2 * ratio6)) + kappa*r0*rat0/(1-pow(rat0,2))));

    }
    else {
        FENE_force = (diffvec / diff)*kappa*r0*rat0/(1-pow(rat0,2));
    }
}*/

Tensor<double,2> FENE_force(int np, Tensor<double,3> &diffvectens,Tensor<double,2> &diff2vals, double d, double kappa, double boltzmann, double r0) {

//****************************************************************************80
//
//  Purpose:
//
//    FENE_force calculates the force contribution due to a finitely extensible elastic spring
//

//  Parameters:
//
//    Input, int np, the number of particles.
//
//    Input, Tensor<double,3> diffvectens,  diffvectens(i,j,alpha) = r(i,alpha)-r(j,alpha)

//    Input, Tensor<double,3> diffvec2, diffvec2(i,j,alpha) = ||r(i,alpha)-r(j,alpha)||^2
//
//    Output, Tensor<double,2> pos(nd,np), the positions.
//
//    Output, Tensor<double,2> vel(nd,np), the velocities.
//
//    Output, Tensor<double,2> acc(nd,np), the accelerations.
//

    Tensor<double,2> fene(np,3);
    fene.setZero();

    double weekscut2 = pow(1.12246204831 * d,2); // 1.122... is 2^(1/6), this is the Weeks-Chandler cutoff distance squared
    double r02 = pow(r0,2);

    for (int i = 0; i < np; i++) {
        for (int alpha = 0; alpha < 3; alpha++){
            if (i==0) {

                double rat0 = pow(diff2vals(i,i+1),0.5)/r02; // this is (ri,i+1/R0)^2, which is the FENE spring part of the force
                if (diff2vals(i,i+1) < weekscut2) {
                    double ratio6 = pow(d / rat0, 3); // this is (d/rij)^6.
                    fene(i, alpha) =
                            diffvectens(i,i+1,alpha) / diff2vals(i,i+1) * ((24 * boltzmann /
                            diff2vals(i,i+1) * ratio6 * (1 - 2*ratio6)) + kappa*r0*rat0/(1-pow(rat0,2)));

                }
                else {
                    fene(i, alpha) = kappa*r0*rat0/(1-pow(rat0,2));
                }

            }
            else {
                double rat1 = pow(diff2vals(i,i-1),0.5)/r0; // this is (ri,i-1/R0)^2, which is the FENE spring part of the force
                double rat2 = pow(diff2vals(i,i+1),0.5)/r0; // this is (ri,i+1/R0)^2, which is the FENE spring part of the force
                if (diff2vals(i,i+1) < weekscut2) {
                    double ratio6 = pow(d / diff2vals(i, i + 1), 6); // this is (d/rij)^6.
                    fene(i, alpha) =
                            diffvectens(i,i+1,alpha) / diff2vals(i,i+1) * ((24 * boltzmann /
                                    diff2vals(i, i + 1) * ratio6 * (1 - 2 * ratio6)) + kappa*r0*rat1/(1-pow(rat1,2)));
                }
                else {
                    fene(i, alpha) = kappa*r0*rat1/(1-pow(rat1,2));
                }
                if (diff2vals(i,i-1) < weekscut2) {
                    double ratio6 = pow(d / diff2vals(i, i - 1), 6); // this is (d/rij)^6.
                    fene(i, alpha) = fene(i, alpha) +
                                           diffvectens(i,i-1,alpha) / diff2vals(i,i-1) * ((24 * boltzmann /
                                                   diff2vals(i,i-1) * ratio6 * (1 - 2 * ratio6)) + kappa*r0*rat2/(1-pow(rat2,2)));
                }
                else {
                    fene(i, alpha) = kappa*r0*rat1/(1-pow(rat2,2));
                }
            }
        }
    }
    return fene;
}