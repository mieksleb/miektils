# Paktep

PAKTEP is a fortran90 code which calculates several geometric properties of double stranded DNA. Currently it only reads oxDNA trajectories but extensions to other files are in the works.

Quantities to calculate include Twist, Writhe, Persistence Length and Plectoneme position. All scripts utilize the Dierckx spline fitting routines with the plectoneme position algorithm using a planar self-intersection based method. 

## Compilation
If you haven't already you will need to set your $MIEKTILS path.

First you will need to compile the static library for the Dierckx spline fitting routines, to do this simply enter the directory "direckx" and type in 
```bash
cd direckx
gfortran -c *.f
ar cr libdierckx.a *.o
```

This should create the static library file libdierckx.a. Here gfortran is my choice of compiler, another will suffice.

Finally we return to the PAKTEP main directory "paktep" and enter 
```bash
make
```

If the compilation fails or you wish to retry, simply enter:
```bash
make clean
```
to your terminal to remove all compiled objects.


## Usage
To choose which observables you wish to calculate, simply set the boolean/logical variable inside processor.f90 to be .True. if you wish to calculate it. This will require a recompilation however. For example if you wish to calcuate the twist and writhe, but not the plectoneme position or persistance length. You would set:

```fortran
twist_writhe_logical=.True.
plectoneme_position_logical=.False.
persistence_length_logical=.False.
```

To run PAKTEP simply use the follwing commands:

```bash
./processor.exe
trajectory.dat
conf.conf
top.top
```
This should create up to 2 files: twist_writhe.dat and plectoneme.dat. If you set persistance_length_logical=.True., then the result is printed in the terminal.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

