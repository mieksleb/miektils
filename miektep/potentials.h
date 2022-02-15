//
// Created by Michael Selby on 10/09/2021.
//

#ifndef MIEKTEP_POTENTIALS_H
#define MIEKTEP_POTENTIALS_H

void KRATKY_pot(double &pot, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
                 double &d21, double &d22, double A);

void FENE_pot( double &pot, int i, Tensor<double,1> &r1, Tensor<double,1> &r2,
               double &d21, double &d22, double kappa, double r0, bool left, bool right);

void WCA_pot(double &pot, Tensor<double,1> &r, double &d2);

class potentials {

};


#endif //MIEKTEP_POTENTIALS_H
