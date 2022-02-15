//
// Created by Michael Selby on 10/09/2021.
//

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
#include "potentials.h"


using namespace std;
using namespace Eigen;


void KRATKY_pot( double &pot, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                   double &d21, double &d22, double A) {
    // potential of particle i due to Kratky-porod hinge where r1 and r2 point away from particle i

    // create u1 and u2 which are normalized versions of r1 and r2
    Tensor<double,1> u1(3);
    Tensor<double,1> u2(3);

    double d1 = pow(d21,0.5);
    u1 = r1/d1;
    double d2 = pow(d22,0.5);
    u2 = r2/d2;
    double dott = dot(u1,u2);

    pot += A * (1-dott);

}

void FENE_pot( double &pot, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                 double &d21, double &d22, double kappa, double r0, bool left, bool right) {

//****************************************************************************80
//
//  Purpose:
//
//    FENE_pot calculates the potential contribution due to a finitely extensible non-elastic spring and nearest neighbour WCA
//
//  Parameters:
//
//    Input, int np, the number of particles.
//
//    Input, double, pot
//
//    Input, boolean left, true if system is connected on the left (i-1 <-> i)
//    Input, boolean right, true if system is connected on the left (i <-> i+1)


    double r02 = pow(r0, 2);     // R0^2
    double k =  kappa * r02 / 2;
    double ratio1 = 1 / d21;
    double rat61 = pow(ratio1, 3);
    double ratio2 = 1 / d22;
    double rat62 = pow(ratio2, 3);

    double rat1 = d21 / r02;
    double rat2 = d22 / r02;
    if (left == TRUE and right == FALSE) {
        pot += 4 * rat61 * (rat61 - 1) + 1 - k * log(1 - rat1);
    }
    else if (left == FALSE and right == TRUE) {
        pot += 4 * rat62 * (rat62 - 1) + 1 - k * log(1 - rat2);
    }
    else {
        pot += 4 * rat61 * (rat61 - 1) + 4 * rat62 * (rat62 - 1) + 2 - k * ( log(1 - rat1) + log(1 - rat2) );
    }
}

void WCA_pot(double &pot, Tensor<double,1> &r, double &d2) {
    /// calculates the Weeks-chander-anderson potential between two particles with separation vector r and d2=|r|^2

    double ratio = 1 / d2;               // (d/rij)^2
    double rat6 = pow(ratio, 3); // (d/rij)^6
    pot += 4 * rat6 * (rat6 - 1)+ 1;

}