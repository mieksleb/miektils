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

int main ( int argc, char *argv[] );
void compute ( int np, Tensor<double,2> pos, Tensor<double,2> vel,
               Tensor<double,2> force, double &pot, double &kin );
double cpu_time ( );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
void timestamp ( );
void update ( int np, double pos[], double vel[], double f[],
              double acc[], double mass, double dt );


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
//    double *acc;
    double ctime;
    double dt;
    double e0;
    Tensor<double,2> force;
    double kinetic;
    double mass = 1.0;
    int np;
    Tensor<double,2> pos;
    double potential;
    int step;
    int step_num;
    int step_print;
    int step_print_index;
    int step_print_num;
    Tensor<double,2> vel;

    timestamp ( );
    cout << "\n";
    cout << "Roachtep\n";
    cout << "  C++ version\n";
    cout << "  A molecular dynamics program for basic roaches :-) \n";
//

//  Get the number of particles.
//
    if ( 1 < argc )
    {
        np = stoi ( argv[2] );
    }
    else
    {
        cout << "\n";
        cout << "  Enter np, the number of particles (500, for instance).\n";
        cin >> np;
    }
//
//  Get the number of time steps.
//
    if ( 2 < argc )
    {
        step_num = atoi ( argv[3] );
    }
    else
    {
        cout << "\n";
        cout << "  Enter step_num, the number of time steps.\n";
        cin >> step_num;
    }
//
//  Get the time step.
//
    if ( 3 < argc )
    {
        dt = atof ( argv[4] );
    }
    else
    {
        cout << "\n";
        cout << "  Enter dt, the time step size.\n";
        cin >> dt;
    }
//
//  Report.
//
    cout << "\n";
    cout << "  NP, the number of particles in the simulation is " << np << "\n";
    cout << "  step_num, the number of time steps, is " << step_num << "\n";
    cout << "  dt, the size of each time step, is " << dt << "\n";


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

    step_print = 0;
    step_print_index = 0;
    step_print_num = 10;

    ctime = cpu_time ( );

    for ( step = 0; step <= step_num; step++ )
    {
        if ( step == 0 )
        {
            initialize(np,pos,vel);
        }
        else
        {
            update(np,pos,vel,force,mass,dt);
        }

        compute(np,pos,vel,mass,force,potential,kinetic );

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

// This function slices a the ith column vector from a rank 2 tensor
// Input: Rank 2 Tensor (n x m matrix)
// Output: Rank 1 tensor (n dimensional vector)

Tensor<double,1> get_slice(Tensor<double,2> input, int i) {
    Eigen::array<Eigen::Index, 2> dims = input.dimensions();
    Eigen::array<long, 2> offsets = {i, 0};
    Eigen::array<long, 2> extents = {1, dims[1]};
    Tensor<double, 1> slice = input.slice(offsets, extents).reshape(Eigen::array<long, 1>{3});
    return slice;
}


//****************************************************************************80

void compute ( int np, Tensor<double,2> pos, Tensor<double,2> vel, double mass,
               Tensor<double,2> force, double &pot, double &kin ) {

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE computes the forces,
//    force : total force, f(i,alpha) is ath component of the vectorial force acting on particle i
//

    double d;
    double d2;
    double PI2 = 3.141592653589793 / 2.0;
    Tensor<double,2> FENE_force;
    Tensor<double,2> BEND_force;
    Tensor<double,1> diffvec;       // matrix ri-rj
    Tensor<double,1> diff2;



    double pot = 0.0;
    double kin = 0.0;

    for ( int i = 0; i < np; i++ ) {
        for ( int j = 0; j < np; j++ ) {
            double alpha;
            if (i == j){
            }
            else{
                diffvec = get_slice(pos,i)-get_slice(pos,j); // vector ri-rj
                Eigen::array<Eigen::IndexPair<int>, 1> product_dims = {Eigen::IndexPair<int>{0, 0}};
                Tensor<double,0> r2;
                diff2 = diffvec.contract(diffvec, product_dims); // norm |ri-rj|^2
                force(i,alpha) = force(i,alpha) + FENE_force(int np, Tensor<double,1> diffvec,Tensor<double,2> diff2, double d, double kappa, double boltzmann, double r0);
            }





            }
        }
//
//  Compute the kinetic energy.
//
        for ( i = 0; i < 3; i++ )
        {
            kin = kin + vel[] * vel[];
        }
    }

    kin = kin * 0.5 * mass;

    return;

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




void initialize ( int np, int nd, Tensor<double,2> pos, Tensor<double,2> vel, Tensor<double,2> acc )

//****************************************************************************80
//
//  Purpose:
//
//    INITIALIZE initializes the positions, velocities, and accelerations.
//

//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Output, double POS[ND*NP], the positions.
//
//    Output, double VEL[ND*NP], the velocities.
//
//    Output, double ACC[ND*NP], the accelerations.
//
{
    int i;
    int alpha;
    int seed;
//
//  Set the positions.
//
    seed = 42069;
    r8mat_uniform_ab ( nd, np, 0.0, 10.0, seed, pos );
//
//  Set the velocities.
//
    for ( i = 0; i < np; i++ )
    {
        for ( alpha = 0; alpha < nd; alpha++ )
        {
            vel(i,alpha) = 0.0;
        }
    }
//
//  Set the accelerations.
//
    for ( i = 0; i < np; i++ )
    {
        for ( alpha = 0; alpha < nd; alpha++ )
        {
            acc(i,alpha) = 0.0;
        }
    }
    return;
}
//****************************************************************************80

void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A, B, the limits of the pseudorandom values.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
    int i;
    const int i4_huge = 2147483647;
    int j;
    int k;

    if ( seed == 0 )
    {
        cerr << "\n";
        cerr << "R8MAT_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit ( 1 );
    }

    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            k = seed / 127773;

            seed = 16807 * ( seed - k * 127773 ) - k * 2836;

            if ( seed < 0 )
            {
                seed = seed + i4_huge;
            }

            r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
        }
    }

    return;
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

vector<double> get_noise() {
    // this function returns a gaussian random variable with zero mean and unit variance
    unsigned seed;
    seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(0.0, 1.0);
    vector<double> noise = {distribution(generator),distribution(generator),distribution(generator)};
    return noise;
}



void update( int np, int nd, Tensor<double,2> pos, Tensor<double,2> vel, Tensor<double,2> force, vector<double> mob, vector<double> mobroot,
              vector<double> acc, double mass, double dt)

//****************************************************************************80
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input/output, double POS[ND*NP], the positions.
//
//    Input/output, double VEL[ND*NP], the velocities.
//
//    Input, double F[ND*NP], the forces.
//
//    Input/output, double ACC[ND*NP], the accelerations.
//
//    Input, double MASS, the mass of each particle.
//
//    Input, double DT, the time step.
//
{
    int i;
    int alpha;
    Tensor<double,2> noise; // gaussian white noise with unit variance and zero mean
    double factor;  //square root of boltzmann like 2*k_B*T

    noise = get_noise;


    for ( i = 0; i < np; i++ ) {
        for ( alpha = 0; alpha < nd; alpha++ ) {
            pos(i,alpha) = pos(i,alpha) + dt*vel(i,alpha)+ (0.5*pow(dt,2) / mass)*force(i,alpha);
            vel(i,alpha) = vel(i,alpha) + (dt/2*mass)*(force(i,alpha)-oldforce(i,alpha));



        }
    }

    return;
}

