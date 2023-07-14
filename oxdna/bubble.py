#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:44:16 2022

@author: michaelselby
"""
import sys
import numpy as np

top_name = sys.argv[1]
conf_name = sys.argv[2]
file_name = sys.argv[3]

toplines = open(top_name,"r").readlines()

line0 = toplines[0].split()
n = int(line0[0])
strands = int(line0[1])

bp = int(n/strands)

if int(toplines[1].split()[2]) == -1:
    circular = False
else:
    circular = True


vec1 = np.zeros((bp,3))
vec2 = np.zeros((bp,3))
step = 0
with open(file_name,"w") as bubble_file:
    with open(conf_name,"r") as conf:
        count = 0
        for line in conf:
            split_line = line.split()
            if split_line[0]=="t" or split_line[0]=="E" or split_line[0]=="b":
                continue
            else:
                x,y,z = float(split_line[0]),float(split_line[1]),float(split_line[2])
                bx, by, bz = float(split_line[3]),float(split_line[4]),float(split_line[5])
                v = np.array([bx,by,bz])
                v /= np.linalg.norm (v)
                if count < bp:
                    vec1[count,:] = v
                else:
                    vec2[count-bp,:] = v
                    
                if count == 2 * bp -1:
                    vec2 = vec2[::-1,:]
                    step += 1
                    dot = [ np.dot(vec1[i,:], vec2[i,:]) for i in range(bp) ]
                    # ang = [ 180 * np.arccos( - x) / np.pi for x in dot ]
                    bubble = [ 0 if x < -0.5 else 1 for x in dot ]
                    bubble_pos = [ i for i in range(bp) if bubble[i] == 1 ]
                    n_bubs = len(bubble_pos)
                    bubble_file.write(str(n_bubs) + "\n")
                    count = 0
                else:  
                    count +=1

    

