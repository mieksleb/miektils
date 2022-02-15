#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <random>
#include <vector>
#define _POSIX_SOURCE
#include <unistd.h>
#undef _POSIX_SOURCE
#include <stdio.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include "forces.h"
#include "potentials.h"
#include "tensor_tools.h"
#include "model.h"

using namespace std;
using namespace Eigen;




void initialize ( int np , string conf_path, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &acc);

//void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
//Tensor<double,2> FENE_force(int np, Tensor<double,2> diffvectens,Tensor<double,2> difftens, double d, double kappa, double boltzmann, double r0);


double cpu_time ( );

void timestamp ( );
void update ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &force,Tensor<double,2> &old_force,
              Tensor<double,2> &acc, double dt, double boltzmann);

void compute ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel,
               Tensor<double,2> &force, double &pot, double &kin, double boltzmann, bool circular );




//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Usage:
//
//
//    where
//
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
    bool circular;
    int step;
    int step_num;
    int step_print;
    int step_print_index;
    int step_print_num;
    string top_file_name;
    string conf_file_name;
    string traj_file_name;
    string last_conf_file_name;

    // get the current working directory, where the inout file is situated
    char cwd_char[1024];
    getcwd(cwd_char, sizeof(cwd_char));
    string cwd = static_cast<string>(cwd_char);
    cwd += static_cast<string>("/");

    string input_file_name = argv[1]; //input file is the only system argument
    string input_file_path = cwd + input_file_name;


    timestamp ( );
    cout << "\n";
    cout << "  Miektep" << endl;
    cout << "  C++ version" << endl;
    cout << "  A molecular dynamics program for coarse grained DNA" << endl;
    cout << "  Current working directory is " << cwd << endl;
//

    string line;
    ifstream input_file;
    string inp_delim = "=";
    input_file.open(input_file_path);
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
               else if ( search == "traj" )
               {
                   traj_file_name = val;
               }
               else if ( search == "last_conf" )
               {
                   last_conf_file_name = val;
               }


            }
        }
    }

    input_file.close(); // close the file object.

    double boltzmann = kb*temp; // kT boltzmann factor

    // assume topology and configuration files are in same directory as input file
    string top_path = cwd + top_file_name;
    string conf_path = cwd + conf_file_name;

//  We now open the topology file and obtain np and circular;

    ifstream top_file;
    top_file.open(top_path);
    if (top_file.is_open()) {   //checking whether the file is open
        getline(top_file, line); // read the first line of the topology, which should contain the number of particles
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

        np = stoi(line);
        string s;
        getline(top_file, line); // read the second line
        istringstream line0(line);

        getline(line0, s,' ');
        getline(line0, s,' ');
        getline(line0, s,' ');

        int num = stoi(s);
        if (num==-1) {
            circular = false;
        }
        else {
            circular = true;
        }

        while (getline(top_file, line)) { //read data from file object and put it into string
            istringstream line0(line);
            while (getline(line0, s,' ')) {
            }

        }
    }
    top_file.close(); // close the file object.



//
//  Report.
//
    cout << "\n";
    cout << "  np, the number of particles in the simulation is " << np << endl;
    cout << "  step_num, the number of time steps, is " << step_num << endl;
    cout << "  dt, the size of each time step, is " << dt << endl;
    cout << "  T, the temperature of each time step, is " << temp << endl;
    if (circular == true) {
        cout << "  circular DNA detected" << endl;
    }
//    cout << "  The input configuration is " + top_file << endl;
//    cout << "  The input topology is " << conf_file << endl;


//
//  This is the main time stepping loop:
//    Compute forces and energies,
//    Update positions, velocities, accelerations.
//

    string energy_file_name = "energy.out";
    // open all files to be written to
    ofstream traj_file;
    ofstream last_conf_file;
    ofstream energy_file;
    traj_file.open(traj_file_name);
    traj_file << np << endl;
    traj_file << endl;
    last_conf_file.open(last_conf_file_name);
    last_conf_file << np << endl;
    last_conf_file << endl;
    energy_file.open(energy_file_name);

    energy_file << "\n";
    energy_file << "  At each step, we report the potential and kinetic energies.\n";
    energy_file << "  The sum of these energies should be a constant.\n";
    energy_file << "  As an accuracy check, we also print the relative error\n";
    energy_file << "  in the total energy.\n";
    energy_file << "\n";
    energy_file << "      Step      Potential       Kinetic        (P+K-E0)/E0\n";
    energy_file << "                Energy P        Energy K       Relative Energy Error\n";
    energy_file << "\n";

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
        potential = 0;
        kinetic = 0;

        if ( step == 0 )
        {
            initialize(np,conf_path,pos,vel,acc);
        }
        else
        {
            update(np,pos,vel,force,old_force,acc, dt, boltzmann);
        }

        compute(np,pos,vel,force,potential,kinetic, boltzmann,circular );


        //print_rank_2_tensor(force);

        if ( step == 0 )
        {
            e0 = potential + kinetic;
        }

        if ( step % step_print == 0 )
        {
            for ( int i = 0; i < np; i++) {
                traj_file << "C " << pos(i,0) << " " << pos(i,1) << " " << pos(i,2) << endl;
            }


            energy_file << "  " << setw(8) << step
                 << "  " << setw(14) << potential
                 << "  " << setw(14) << kinetic
                 << "  " << setw(14) << ( potential + kinetic - e0 ) / e0 << endl;
            step_print_index = step_print_index + 1;
            step_print = ( step_print_index * step_num ) / step_print_num;
        }
        if ( step == step_num ) {
            for ( int i = 0; i < np; i++) {
                last_conf_file << "C " << pos(i,0) << " " << pos(i,1) << " " << pos(i,2) << endl;
            }

        }

    }
    traj_file.close();
    last_conf_file.close();
    energy_file.close();
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
               Tensor<double,2> &force, double &pot, double &kin, double boltzmann, bool circular ) {

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
    Tensor<double, 2> kratky(np,3);
    Tensor<double, 1> ri(3);
    Tensor<double, 1> r1(3);
    Tensor<double, 1> r2(3);
    Tensor<double, 1> rij(3);
    Tensor<double, 0> diff2;
    Tensor<double, 3> diffvectens(np,np,3);
    Tensor<double, 2> diff2vals(np,np);
    Tensor<double,1> vel1;
    double d22;
    double d21;


    diffvectens.setZero();
    diff2vals.setZero();
    kratky.setZero();
    fene.setZero();
    vel1.setZero();
    kin = 0;
    pot = 0;



//  calculate the separation tensors, all of which are symmetric for i <-> j
//  all diagonal elements are zero and since tensors are initialized to zero we can ignore them
    for (int i = 0; i < np; i++) {
        ri = get_slice_rank2(pos, i);

        for (int j = 0; j < np; j++) {

            if (i < j) {
                // if i=!j then we must calculate
                rij = get_slice_rank2(pos, i) - get_slice_rank2(pos, j); // vector ri-rj
                Eigen::array<IndexPair<int>, 1> product_dims = {IndexPair<int>{0, 0}};
                diff2 = rij.contract(rij, product_dims); // norm ||ri-rj||^2
                d2 = diff2();
                diff2vals(i,j) = d2;

                for (int alpha = 0; alpha < 3; alpha++) {
                    diffvectens(i, j, alpha) = rij(alpha);

                }
            }

            else {
                // tensors have symmetry/anti-symmetry!
                diff2vals(i,j) = diff2vals(j,i);
                for (int alpha = 0; alpha < 3; alpha++) {
                    diffvectens(i, j, alpha) =  - diffvectens(j, i, alpha);

                }
            }

        }
    }

    // now we can calculate all quantities we wish from the separation tensors

    for (int i = 0; i < np; i++) {

        vel1 = get_slice_rank2(vel,i);
        kin += dot(vel1,vel1);

        // now we must account for when we are at the end of the chains as some potentials are "one-sided" if system is not circular
        if (i == 0) {
            r1 = get_slice_rank3(diffvectens,0,np-1); // vector r(np-1)-r0
            r2 = get_slice_rank3(diffvectens,0,1);
            d21 = diff2vals(np-1,0);
            d22 = diff2vals(0,1);


            if (circular == TRUE) {
                FENE_force(fene, i, r1, r2, d21, d22, kappa, r0, TRUE,
                           TRUE);   // adds fene force due to particle 1 and np
                FENE_pot(pot, i, r1, r2, d21, d22, kappa, r0, TRUE,
                           TRUE);   // adds fene force due to particle 1 ONLY

                //KRATKY_force(kratky, i, r1, r2, d21, d22,
                             //A);                // adds kratky-porod force due to particles 1 and np

                //KRATKY_pot(pot, i, r1, r2, d21, d22,A);


            }
            else {
                FENE_force(fene, i, r1, r2, d21, d22, kappa, r0, FALSE,
                           TRUE);   // adds fene force due to particle 1 ONLY
                FENE_pot(pot, i, r1, r2, d21, d22, kappa, r0, FALSE,
                           TRUE);   // adds fene force due to particle 1 ONLY
            }
        }
        else if (i == np-1) {
            r1 = get_slice_rank3(diffvectens,np-1,np-2);
            r2 = get_slice_rank3(diffvectens,np-1,0);
            d21 = diff2vals(np-2,np-1);
            d22 = diff2vals(np-1,0);


            if (circular == TRUE) {

                FENE_force(fene, i, r1, r2, d21, d22, kappa, r0, TRUE,
                           TRUE);   // adds fene force due to particle 1 and np
                FENE_pot(pot, i, r1, r2, d21, d22, kappa, r0, TRUE,
                           TRUE);   // adds fene force due to particle 1 and np


                //KRATKY_force(kratky, i, r1, r2, d21, d22,
                             //A);                // adds kratky-porod force due to particles 1 and np
                //KRATKY_pot(pot, i, r1, r2, d21, d22,A);
            }
            else {

                FENE_force(fene, i, r1, r2, d21, d22, kappa, r0, TRUE,
                           FALSE);   // adds fene force due to particle np-1 ONLY
                FENE_pot(pot, i, r1, r2, d21, d22, kappa, r0, TRUE,
                           FALSE);   // adds fene force due to particle np-1 ONLY


            }
        }
        else {
            r1 = get_slice_rank3(diffvectens,i,i-1);
            r2 = get_slice_rank3(diffvectens,i,i+1);
            d21 = diff2vals(i,i-1);
            d22 = diff2vals(i,i+1);


            FENE_force(fene, i, r1, r2, d21, d22, kappa, r0, TRUE, TRUE);
            FENE_pot(pot, i, r1, r2, d21, d22, kappa, r0, TRUE, TRUE);

            //KRATKY_force(kratky, i, r1, r2, d21, d22, A);
            //KRATKY_pot(pot, i, r1, r2, d21, d22, A);

        }
    }




    //wca = WCA_force(np, diffvectens, diff2vals, d, kappa, boltzmann, r0);
    force = fene ;
    cout << "force is" << endl;
    print_rank_2_tensor(fene);
    cout << endl;


//
//  Compute the kinetic energy.
//

    kin *= 0.5 * mass;

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
        string s;
        int part_counter = 0; // this is the counter for the positional part of pos
        while (getline(conf_file, line)) { //read data from file object and put it into string
            istringstream line0(line);
            int counter = 0;
            while (getline(line0, s, ' ')) {
                pos(part_counter,counter) = stod(s);
                counter += 1;  // increment the positional counter by one
            }
            part_counter += 1; // increment particle counter by one

        }
    }
    conf_file.close(); //close the file object.
    //vel = get_noise(np);

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





void update ( int np, Tensor<double,2> &pos, Tensor<double,2> &vel, Tensor<double,2> &force, Tensor<double,2> &old_force,
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
    // input timestep is in ns so must convert to simulation units
    //dt *= pow(10,-9)/tau;
    dt = 0.02;

    Tensor<double,2> noise(np,3); // gaussian white noise with unit variance and zero mean
    double halfdt = dt/2;
    double gamt = pow(gam*dt,0.5);
    Tensor<double,2> t(np,3);

    // velocity verlet algorithm with langevin thermostat

    //generate the noise for the half time-step interval
    noise = get_noise(np);

    noise *= t.constant(gamt);

    cout << "noise is" << endl;
    print_rank_2_tensor(noise);
    cout << endl;

    // first update velocity for half time-step
    vel += halfdt * ( force - vel * t.constant(gam) ) + noise;

    //generate the noise for the full time-step interval
    noise = get_noise(np);
    noise *= t.constant(gamt);
    // then we update after whole time-step
    pos += vel * t.constant(dt);
    vel += halfdt * ( force - vel * t.constant(gam) ) + noise;

    cout << "velocity is" << endl;
    print_rank_2_tensor(vel);
    cout << endl;


    old_force = force;


    return;
}


