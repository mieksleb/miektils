#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:13:05 2020

This file splits a trajectory file into a list of configuration files, 
saving these files in a temporary directory called 'bin'. We then loop over these config files
and calcuate properties for each. The temporary directory is deleted at the end.

@author: michaelselby
"""

import os
import tools
import readers
#import twist_writhe
import neme_hunter
import shutil
import numpy as np



parent_dir = os.path.dirname(os.path.realpath(__file__))
bin_dir = "bin/"
path = os.path.join(parent_dir, bin_dir) 

try:
    os.mkdir(path)
except Exception:
    pass

output_file = open("spline_analysis.dat", "w")
output_file.write("twist"+" "+"writhe"+" "+"neme1"+" "+"neme2"+"\n")

traj_file = "miekbin/trajectory624.dat"
conf = "miekbin/last_conf624_2_(plect).dat"
init_conf = "miekbin/last_conf624_2.dat"
last_conf = "miekbin/last_conf624_2.dat"
top = "miekbin/sim624.top"
plect_file = open("miekbin/plectoneme624.dat","r").readlines()


"""
circles
"""
#top = "miekbin/circle_400-Lk2.top"
#conf = "miekbin/last_conf_TD.dat"
##conf = "miekbin/circle_400-Lk2.dat"

#traj_file = "miekbin/trajectory_300-Lk7.dat"
##init_conf = "miekbin/last_conf.dat"
#init_conf = "miekbin/circle_300-Lk7.dat"
#top = "miekbin/circle_300-Lk7.top"

num_lines = sum(1 for line in plect_file)

l = readers.LorenzoReader(init_conf, top)
s = l.get_system()
strand1 = s._strands[0]
strand2 = s._strands[1]
bp = len(strand1._nucleotides)
steps = 251000000          #number of simulation steps
printed_steps = num_lines  #results printed every steps
lines_per_print = int(steps/printed_steps)  #number of steps wanted for this calculation
print(printed_steps)
calc_steps = 10
npoints = 1000


'''
Converts a trajectory file into mutiple configuration files so that they can 
be passed through various routines
'''
def traj_2_confs(trajectory, bp, printed_steps, energy_out=True):
    list_of_confs = []
    if energy_out == True:
        lines_per_file = 2*bp+3
    else:
        lines_per_file = 2*bp
    smallfile = None
    with open(trajectory) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % int(lines_per_file*printed_steps*calc_steps) == 0:
                if smallfile:
                    smallfile.close()
                small_filename = 'conf'+str(int(lineno/lines_per_file))
                print(int(lineno/lines_per_file))
                smallfile = open(os.path.join(path, small_filename), "w")
            smallfile.write(line)
        list_of_confs.append(smallfile)
        if smallfile:
            smallfile.close()
        np.insert(list_of_confs,0,init_conf)
        np.insert(list_of_confs,-1,last_conf)   
        return list_of_confs

traj_list = traj_2_confs(traj_file, bp, printed_steps, energy_out=True)



for i in range(0, calc_steps):
    print("step", i*calc_steps)
    conf = os.path.join(path, "conf"+str(i))
    step =  i*int(steps/calc_steps)
    top = top  
    l = readers.LorenzoReader(conf, top)
    s = l.get_system()  
    strand1 = s._strands[0]
    strand2 = s._strands[1]
    line = plect_file[i*int(steps/calc_steps)].split()
    neme_guess = line[3]
    spline1 = tools.get_bb_spline(strand1)
    spline2 = tools.get_bb_spline(strand2, reverse = True)   
#    twist, writhe = twist_writhe.get_twist_writhe(spline1, spline2)
    nemes = neme_hunter.neme_hunter(spline1, spline2, npoints, bp, neme_guess, circular = False)
    if len(nemes) == 1:
        #output_file.write(str(twist)+" "+str(writhe)+" "+str(nemes[0])+"\n")
        output_file.write(str(step)+" "+str(nemes[0])+ " " + str(neme_guess)+ "\n")
    else:   
        output_file.write(str(twist)+" "+str(writhe)+" "+str(nemes[0])+" "+str(nemes[1])+"\n")
    
    
output_file.close()
shutil.rmtree(path)



