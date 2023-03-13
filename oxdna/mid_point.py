#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 09:35:39 2023

@author: michaelselby
"""
import sys
import numpy as np


conf = sys.argv[1]
top = sys.argv[2]
xyz = sys.argv[3]


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
                    r_long = (strand1pos + strand2pos) / 2
                    for r1 in r_long:
                        line = "C "+ str(r1[0])+" "+str(r1[1])+" " +str(r1[2])+ "\n"        
                        mol_cont_file.write(line)
                    mol_cont_file.write("\n")
                    mol_cont_file.write("\n")
                    count = 0
                    step += 1
                    if step % 1000 == 0:
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

