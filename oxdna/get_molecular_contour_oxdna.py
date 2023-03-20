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
method = sys.argv[4]

# conf = "/Users/michaelselby/Documents/DPhil/miektils/oxdna/minicircles/circle_336/deltaLk0/MDNA/relax_trajectory_336_0.dat"
# top = "/Users/michaelselby/Documents/DPhil/miektils/oxdna/minicircles/circle_336/deltaLk0/MDNA/circ_336_0.top"
# xyz = "/Users/michaelselby/Documents/DPhil/miektils/oxdna/minicircles/circle_336/deltaLk0/MDNA/circ_336_0.xyz"
# circ = 1

alpha = 0.06

if int(method) == 1:
    wrline = True
else:
    wrline = False

def get_mid_line(bp,r1,r2,n1,n2,b1,b2,alpha):
    mid = [(r1[i,:]+r2[i,:])/2 + alpha * ( np.cross(b1[i,:],n1[i,:]) + np.cross(b2[i,:],n2[i,:]))/2 for i in range(1,bp-1)]
    return mid

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
n1 = np.zeros((bp,3))
b1 = np.zeros((bp,3))
n2 = np.zeros((bp,3))
b2 = np.zeros((bp,3))
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
                bx,by,bz = float(split_line[3]),float(split_line[4]),float(split_line[5])
                nx,ny,nz = float(split_line[6]),float(split_line[7]),float(split_line[8])
                if count < bp:
                    strand1pos[count,:] = np.array([x,y,z])
                    n1[count,:] = np.array([nx,ny,nz])
                    b1[count,:] = np.array([bx,by,bz])
                else:
                    strand2pos[count-bp-1,:] = np.array([x,y,z])
                    n2[count-bp-1,:] = np.array([nx,ny,nz])
                    b2[count-bp-1,:] = np.array([bx,by,bz])
                    
                if count == 2*bp -1:
                    strand2pos = strand2pos[::-1,:]
                    if wrline:
                        r_long = get_molecular_contour(bp, strand1pos, strand2pos, circular=circular)
                    else:
                        r_long = get_mid_line(bp, strand1pos, strand2pos, n1, n2, b1, b2, alpha)
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
                    
    
    

