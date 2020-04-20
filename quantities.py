#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:13:53 2019

python script with many different calculable quantities


@author: michaelselby
"""

import math
import numpy as np
import scrape
import tools
import base
import readers


#reads in trajectory file
init_file = open("circle_300-Lk7.dat", "r").readlines()
input_file = open("CPU_input", "r").readlines()
traj_file = open("trajectory.dat", "r").readlines()
top_file = open("generated.top", "r").readlines()

init_file_path = "circle_120-Lk4.dat"
top_file_path = "generated.top"


bp = int((len(init_file)-3)/2) #number of base pairs
N = 2*bp
steps = float(scrape.INPUT(input_file).steps)                       #number of simulation steps
conf_steps = float(scrape.INPUT(input_file).print_conf_interval)    #number of print_conf_interval steps
printed_steps = int(steps/conf_steps)                               #number of steps where we have printed


traj_snippets = tools.snippets(traj_file, N, printed_steps)
#print(traj_snippets[0])




def twist(file, bp, printed_steps): 
    for t in range(printed_steps):
        strand1 = traj_snippets[t][:bp]                #coords of boths strands
        strand2 = traj_snippets[t][bp:]
        pairs = []
        for i in range(len(strand1)):        #coords paired by base pairs
            pairs.append((strand1[i], strand2[-i-1]))
        vecs = []
        for i in range(len(pairs)):          #vectors of base pairs from strand 1 to strand 2
            vecs.append(pairs[i][0] - pairs[i][1])
        indexlist = range(len(vecs))         #indexes which we need to itterate the dot product over 
        dots = []
        for i in range(len(vecs)-1):
            dots.append((np.dot(vecs[indexlist[i]], vecs[indexlist[i+1]])/(np.linalg.norm(vecs[indexlist[i]]) * np.linalg.norm(vecs[indexlist[i+1]]))))
        tw = 0
        for j in range(len(dots)):  
            tw += np.arccos(dots[j])/(2 * math.pi)
    
        return tw
    
twist = twist(init_file, bp, printed_steps)









#twost = twist(traj_file)
#print(twost)