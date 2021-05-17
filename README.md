# Miektils

Welcome to Miektils. Miektils is a collection of utilities for DNA simulations, designed specifically for oxDNA coarse-grained simulations and all-atom systems from Amber. Submodules include PAKTEP, Ambermiek, oxdna and WLC.


## PAKTEP

PAKTEP is FORTRAN90 code which calculates the Twist, Writhe, Plectoneme position and Persistence Length of DNA (and potentially other macromolecules). PAKTEP utilises spline fitting routines and boasts competitive performance for long trajectory processing.



## Ambermiek

Ambermiek is a file generation system designed specifically for creating linear and circular dsDNA simulation files. Circular DNA can be created with an arbirary Linking number defecit, a toplogical constraint on the DNA related to its "excess-twist".

## oxdna

Contains python and bash scripts which generate and post-process oxdna simulations including minicircle generation and twist/writhe calculations.

## WLC

This directory includes python scripts which generate piecewise elatica configurations and energies, as well as numerical scripts for 3D elasticae.
