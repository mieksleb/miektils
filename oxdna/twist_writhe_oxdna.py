#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:44:16 2022

@author: michaelselby
"""
import sys
import numpy as np
from tools import get_spline
from tools_fast_math import get_twist_writhe

top_name = sys.argv[1]
conf_name = sys.argv[2]
file_name = sys.argv[3]
npoints = sys.argv[4]

toplines = open(top_name,"r").readlines()

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
with open(file_name,"w") as twist_writhe_file:
    with open(conf_name,"r") as conf:
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
                    spline1 = get_spline(strand1pos)
                    spline2 = get_spline(strand2pos)
    
                    twist, writhe = get_twist_writhe(spline1,spline2,npoints,circular)
                    twist_writhe_file.write(str(twist) + " " + str(writhe) + "\n")
                    count = 0
                    step += 1
                    # print(step)
                else:  
                    count +=1
    
    
    

