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




void FENE_force( Tensor<double,2> &force, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                            double &d21, double &d22, double kappa, double r0, bool left, bool right) {

//****************************************************************************80
//
//  Purpose:
//
//    FENE_force calculates the force contribution due to a finitely extensible non-elastic spring
//
//  Parameters:
//
//    Input, int np, the number of particles.
//
//    Input, Tensor<double,1> force
//
//    Input, boolean left, true if system is connected on the left (i-1 <-> i)
//    Input, boolean right, true if system is connected on the left (i <-> i+1)

    int alpha;

    double r02 = pow(r0, 2);                         // R0^2

    double rat1 = d21 / r02;
    double rat2 = d22 / r02;
    double ratio1 = 1 / d21;
    double rat61 = pow(ratio1, 3);
    double ratio2 = 1 / d22;
    double rat62 = pow(ratio2, 3);

    if (left == TRUE and right == FALSE) {
        for (alpha = 0; alpha < 3; alpha++) {
            force(i, alpha) =   r1(alpha) * ( 24 * rat61 * ratio1 * ( 1 - 2 * rat61 ) + kappa / (1 - rat1));
        }
    }
    else if (left == FALSE and right == TRUE) {
        for (alpha = 0; alpha < 3; alpha++) {
            force(i, alpha) = r2(alpha) * ( 24 * rat62 * ratio2 * (1 -  2 * rat62 ) + kappa / (1 - rat2));
        }

    }
    else {
        for (alpha = 0; alpha < 3; alpha++) {
            force(i, alpha) =   r1(alpha) * ( 24 * rat61 * ratio1 * ( 1 - 2 * rat61 ) + kappa / (1 - rat1)) + r2(alpha) * ( 24 * rat62 * ratio2 * (1 - 2 * rat62 ) + kappa / (1 - rat2));
        }

    }
}


void KRATKY_force( Tensor<double,2> &force, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                 double &d21, double &d22, double A) {
    // Force on particle i due to Kratky-porod hinge where r1 and r2 point away from particle i

    // create u1 and u2 which are normalized versions of r1 and r2
    Tensor<double,1> u1(3);
    Tensor<double,1> u2(3);

    double d1 = pow(d21,0.5);
    u1 = r1/d1;
    double d2 = pow(d22,0.5);
    u2 = r2/d2;
    double dott = dot(u1,u2);

    for (int alpha = 0; alpha < 3; alpha++) {
        force(i, alpha) = A * (u1(alpha) * (1/d2 - dott/d1)+ u2(alpha) * (1/d1 - dott/d2));
    }

}


