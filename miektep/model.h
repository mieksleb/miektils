//
// Created by Michael Selby on 11/08/2021.
//
// All the model constants

#ifndef MIEKTEP_MODEL_H
#define MIEKTEP_MODEL_H


int n = 7;                                               // number of particles which determine core particle
double PI = 3.141592653589793;
double d = n * 3.4 * pow(10,-10);                                       // bead diameter
double gam = 1;                                         // linear drag coefficient
double kb = 1.38064852*pow(10,-23);     // Boltzmann's constant in SI units
double mass = n * 1.66 * pow(10,-27);                                  // mass of a bead in SI units
double kappa = 30;                                  // fene kappa in simulation units
double r0 = 1.6;                                    // r0 = 1.6*d
double A = 150/n;                               // Kratky-Porod stiffness in simulation units
double tau = pow(mass * d * d / kb * 300 ,0.5);

#endif //MIEKTEP_MODEL_H
