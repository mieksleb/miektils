#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:31:31  2023

@author: michaelselby
"""
from tools import disc_curvature
# from tools_fast_math import disc_curvature
import numpy as np
import sys
import matplotlib.pyplot as plt


file = sys.argv[1]
circ = sys.argv[2]

output_file = "curvature.dat"

if int(circ) == 1:
    circular = True
else:
    circular = False


# file = "/Users/michaelselby/Documents/DPhil/miektils/ambermiek/minicircle/Hussain/MDNA/molecular_contour_MDNA.xyz"
# output_file = "/Users/michaelselby/Documents/DPhil/miektils/ambermiek/minicircle/Hussain/MDNA/curvature.dat"
# circular = True


with open(file,"r") as file:
    bp = int(file.readline())
    r = np.empty((bp,3))
    curvature = np.zeros((bp))
    file.readline()
    i = 0
    steps = 0
    
    for line in file:
        split_line = line.split()
        if len(split_line) > 0 and split_line[0] == "C":
            r[i,:] = np.array([float(split_line[1]),float(split_line[2]),float(split_line[3])])
            i += 1
        else:
            file.readline()
            i = 0
            curv_vals = disc_curvature(r, circular)
            curvature += curv_vals
            steps += 1


    curvature /= steps
    
with open(output_file,"w") as file:
    for curv in curvature:
        file.write(str(curv) + "\n")
    