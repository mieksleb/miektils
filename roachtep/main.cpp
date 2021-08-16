#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <random>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include "forces.h"
#include "tensor_tools.h"
#include "model.h"

using namespace std;
using namespace Eigen;




void initialize ( int np , string conf_path, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &acc);

void compute ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel,
              Tensor<double,2> &force, double &pot, double &kin );

//void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
//Tensor<double,2> FENE_force(int np, Tensor<double,2> diffvectens,Tensor<double,2> difftens, double d, double kappa, double boltzmann, double r0);


double cpu_time ( );

void timestamp ( );
void update ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &force,Tensor<double,2> &old_force,
              Tensor<double,2> &acc, double dt, double boltzmann);




//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Usage:
//
//    md nd np step_num dt
//
//    where
//
//    * nd is the spatial dimension (2 or 3);
//    * np is the number of particles (500, for instance);
//    * step_num is the number of time steps (500, for instance).
//    * dt is the time step (0.1 for instance)
//
{
    int np;
    double ctime;
    double dt;
    double e0;
    double kinetic;
    double potential;
    double temp;
    int step;
    int step_num;
    int step_print;
    int step_print_index;
    int step_print_num;
    string top_file_name;
    string conf_file_name;

    cout << argv[0] << endl;

    timestamp ( );
    cout << "\n";
    cout << "  Miektep\n";
    cout << "  C++ version\n";
    cout << "  A molecular dynamics program for basic roaches :-) \n";
//


    string line;
    ifstream input_file;
    string inp_delim = "=";
    input_file.open("/Users/michaelselby/Documents/DPhil/miektils/roachtep/input.in");
    if (input_file.is_open()) {   //checking whether the file is open
        while (getline(input_file, line)) { //read data from file object and put it into string
            auto start = 0U;
            auto end_str = line.find(inp_delim);
            while (end_str != string::npos)
            {
                string search = line.substr(start, end_str - start);
                start = end_str + inp_delim.length();
                end_str = line.find(inp_delim, start);
                string val = line.substr(start, end_str);
                search.erase( std::remove_if( search.begin(), search.end(), ::isspace ), search.end() );
                val.erase( std::remove_if( val.begin(), val.end(), ::isspace ), val.end() );
               if ( search == "steps" )
                {
                    step_num = stoi(val);
                }
                else if ( search == "dt" )
                {
                    dt = stod(val);
                }
                else if ( search == "top" )
                {
                    top_file_name = val;
                }
                else if ( search == "conf" )
                {
                    conf_file_name = val;
                }
                else if ( search == "T" )
                {
                    temp = stod(val);
                }


            }
        }
    }
    input_file.close(); //close the file object.

    double boltzmann = kb*temp;

    string top_path = "/Users/michaelselby/Documents/DPhil/miektils/roachtep/"+top_file_name;
    string conf_path = "/Users/michaelselby/Documents/DPhil/miektils/roachtep/"+conf_file_name;

//  We now open the topology file and obtain np and circular;

    ifstream top_file;
    top_file.open(top_path);
    if (top_file.is_open()) {   //checking whether the file is open
        getline(top_file, line);
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        np = stoi(line);
        string s;
        while (getline(top_file, line)) { //read data from file object and put it into string
            istringstream line0(line);
            while (getline(line0, s,' ')) {

                cout << s << endl;
            }

        }
    }
    top_file.close(); //close the file object.



//
//  Report.
//
    cout << "\n";
    cout << "  np, the number of particles in the simulation is " << np << endl;
    cout << "  step_num, the number of time steps, is " << step_num << endl;
    cout << "  dt, the size of each time step, is " << dt << endl;
    cout << "  T, the temperature of each time step, is " << temp << endl;
//    cout << "  The input configuration is " + top_file << endl;
//    cout << "  The input topology is " << conf_file << endl;


//
//  This is the main time stepping loop:
//    Compute forces and energies,
//    Update positions, velocities, accelerations.
//
    cout << "\n";
    cout << "  At each step, we report the potential and kinetic energies.\n";
    cout << "  The sum of these energies should be a constant.\n";
    cout << "  As an accuracy check, we also print the relative error\n";
    cout << "  in the total energy.\n";
    cout << "\n";
    cout << "      Step      Potential       Kinetic        (P+K-E0)/E0\n";
    cout << "                Energy P        Energy K       Relative Energy Error\n";
    cout << "\n";

    step_print = 10;
    step_print_index = 0;
    step_print_num = 10;

    ctime = cpu_time ( );


    // Construct the positions, velocities, accelerations and forces all as rank 2 tensors
    Tensor<double,2> pos(np,3);
    Tensor<double,2> acc(np,3);
    Tensor<double,2> vel(np,3);
    Tensor<double,2> force(np,3);
    Tensor<double,2> old_force(np,3);

    // Initialize with zeros
    pos.setZero();
    acc.setZero();
    vel.setZero();
    force.setZero();
    old_force.setZero();



    for ( step = 0; step <= step_num; step++ ) {

        if ( step == 0 )
        {
            initialize(np,conf_path,pos,vel,acc);
        }
        else
        {
            update(np,pos,vel,force,old_force,acc,dt,boltzmann);
        }
        compute(np,pos,vel,force,potential,kinetic );

        if ( step == 0 )
        {
            e0 = potential + kinetic;
        }

        if ( step == step_print )
        {
            cout << "  " << setw(8) << step
                 << "  " << setw(14) << potential
                 << "  " << setw(14) << kinetic
                 << "  " << setw(14) << ( potential + kinetic - e0 ) / e0 << "\n";
            step_print_index = step_print_index + 1;
            step_print = ( step_print_index * step_num ) / step_print_num;
        }

    }
//
//  Report timing.
//
    ctime = cpu_time ( ) - ctime;
    cout << "\n";
    cout << "  Elapsed cpu time " << ctime << " seconds.\n";
//

//  Terminate.
//
    cout << "\n";
    cout << "MD\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );

    return 0;
}




//****************************************************************************80

void compute ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel,
               Tensor<double,2> &force, double &pot, double &kin ) {

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE computes the forces,
//    force : total force, f(i,alpha) is alpha'th component of the vectorial force acting on particle i
//

    double d;
    double d2;
    double PI2 = PI / 2.0;
    Tensor<double, 2> fene(np,3);
    Tensor<double, 2> bend(np,3);
    Tensor<double, 1> rij(3);       // matrix ri-rj
    Tensor<double, 0> diff2;
    Tensor<double, 3> diffvectens(np,np,3);
    Tensor<double, 2> diff2vals(np,np);

    diffvectens.setZero();
    diff2vals.setZero();
    fene.setZero();
    d2 = diff2();

    double boltzmann = 1;
    double kappa = 1;
    double r0 = 2;

//  calculate the separation tensors
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            // if i=j then separation vector and scalar is zero
            if (i == j) {
                diff2vals(i, j) = 0;
                for (int alpha = 0; alpha < 3; alpha++) {
                    diffvectens(i, j, alpha) = 0;
                    d2 = 0;
                }
            }
            else {
                // if i=!j then we must calculate
                rij = get_slice(pos, i) - get_slice(pos, j); // vector ri-rj
                Eigen::array<IndexPair<int>, 1> product_dims = {IndexPair<int>{0, 0}};
                diff2 = rij.contract(rij, product_dims); // norm ||ri-rj||^2
                //diff2 = norm2(rij);
                //cout << d2 << endl;
                d2 = diff2();
                diff2vals(i,j) = d2;

                for (int alpha = 0; alpha < 3; alpha++) {
                    diffvectens(i, j, alpha) = rij(alpha);

                }
            }
        }
    }
    fene = FENE_force(np, diffvectens, diff2vals, d, kappa, boltzmann, r0);
    force = fene;


//
//  Compute the kinetic energy.
//

    kin = kin * 0.5 * mass;

    return;

}

//****************************************************************************80


double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
    double value;

    value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

    return value;
}
//****************************************************************************80



void initialize ( int np , string conf_path, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &acc)

//****************************************************************************80
//
//  Purpose:
//
//    INITIALIZE initializes the positions, velocities, and accelerations.
//

//  Parameters:
//
//    Input, int np, the number of particles.
//
//    Output, Tensor<double,2> pos(np,3), the positions.
//
//    Output, Tensor<double,2> vel(np,3), the velocities.
//
//    Output, Tensor<double,2> acc(np,3), the accelerations.
//
{

    int i;
    int alpha;
    int seed;
    double d = 1;
    string line;
//
//  Set the positions, velocities and accelerations
//  We now open the configuration file

    ifstream conf_file;
    conf_file.open(conf_path);
    if (conf_file.is_open()) {   //checking whether the file is open
        getline(conf_file, line);
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        np = stoi(line);
        string s;
        while (getline(conf_file, line)) { //read data from file object and put it into string
            istringstream line0(line);
            while (getline(line0, s, ' ')) {

                cout << s << endl;
            }

        }
    }
    conf_file.close(); //close the file object.
}

//****************************************************************************80


void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    cout << time_buffer << "\n";

    return;
# undef TIME_SIZE
}
//****************************************************************************80





void update ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &force,  Tensor<double,2> &old_force,
                            Tensor<double,2> &acc, double dt, double boltzmann)

//****************************************************************************80
//
//  Parameters:
//
//    Input, int np, the number of particles.
//
//    Output, Tensor<double,2> pos(np,3), the positions.
//
//    Output, Tensor<double,2> vel(np,3), the velocities.
//
//   Output, Tensor<double,2> acc(np,3), the accelerations.
//
//    Output, Tensor<double,2> force(np,3), the forces.
//
//     Output, Tensor<double,2> old_force(np,3), the forces at the previous timestep
//
//    Input, double mass, the mass of each particle.
//
//    Input, double dt, the time step.
//
{
    int i;
    int alpha;
    Tensor<double,2> noise(np,3); // gaussian white noise with unit variance and zero mean
    noise = get_noise(np);
    noise *= noise.constant(2*gam*boltzmann);
    double factor;  //square root of boltzmann like factor 2*k_B*T
    double dt2 = pow(dt,2);
    double fact1 = 0.5*dt2/mass;
    double fact2 = dt/2*mass;
    double b = 1/(1+gam*fact1);
    double a = (1+gam*fact1)*b;

    noise = get_noise(np);


    for ( i = 0; i < np; i++ ) {
        for ( alpha = 0; alpha < 3; alpha++ ) {
            pos(i,alpha) += b*gam*dt*vel(i,alpha)+ b*fact1*force(i,alpha)+ b*fact2*noise(i,alpha);
            vel(i,alpha) += fact2*( a * old_force(i,alpha)+force(i,alpha)) + b*noise(i,alpha) / mass;
            old_force(i,alpha) = force(i,alpha);



        }
    }

    return;
}




