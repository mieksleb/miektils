#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:44:16 2022

@author: michaelselby
"""
import sys
import numpy as np
from molecular_contour import get_molecular_contour


conf = sys.argv[1]
top = sys.argv[2]
xyz = sys.argv[3]
circ = sys.argv[4]

# conf = "/Users/michaelselby/Documents/DPhil/miektils/oxdna/minicircles/circle_336/deltaLk0/MDNA/relax_trajectory_336_0.dat"
# top = "/Users/michaelselby/Documents/DPhil/miektils/oxdna/minicircles/circle_336/deltaLk0/MDNA/circ_336_0.top"
# xyz = "/Users/michaelselby/Documents/DPhil/miektils/oxdna/minicircles/circle_336/deltaLk0/MDNA/circ_336_0.xyz"
# circ = 1



if int(circ) == 1:
    circular = True
else:
    circular = False



toplines = open(top,"r").readlines()

line0 = toplines[0].split()
n = int(line0[0])
strands = int(line0[1])

bp = int(n/strands)

if int(toplines[1].split()[2]) == -1:
    circular = False
else:
    circular = True


strand1pos = np.zeros((bp,3))
strand2pos = np.zeros((bp,3))
step = 0
with open(xyz,"w") as mol_cont_file:
    mol_cont_file.write(str(bp)+ "\n")
    mol_cont_file.write("\n")
    with open(conf,"r") as conf:
        count = 0
        for line in conf:
            split_line = line.split()
            if split_line[0]=="t" or split_line[0]=="E" or split_line[0]=="b":
                continue
            else:
                x,y,z = float(split_line[0]),float(split_line[1]),float(split_line[2])
                if count < bp:
                    strand1pos[count,:] = np.array([x,y,z])
                else:
                    strand2pos[count-bp-1,:] = np.array([x,y,z])
                    
                if count == 2*bp -1:
                    strand2pos = strand2pos[::-1,:]
                    r_long = get_molecular_contour(bp, strand1pos, strand2pos, circular=circular)
                    for r1 in r_long:
                        line = "C "+ str(r1[0])+" "+str(r1[1])+" " +str(r1[2])+ "\n"        
                        mol_cont_file.write(line)
                    mol_cont_file.write("\n")
                    mol_cont_file.write("\n")
                    count = 0
                    step += 1
                    print(step)

                else:  
                    count +=1
                    
    
    
# file.write(str(int(traj.nres/traj.nstrand))+ "\n")
# file.write("\n")
# for long_array in list_of_quants:
#     for r1 in long_array:
#         line = "C "+ str(r1[0])+" "+str(r1[1])+" " +str(r1[2])+ "\n"        
#         file.write(line)
#     file.write("\n")
#     file.write("\n")

