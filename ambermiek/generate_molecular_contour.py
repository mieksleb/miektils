#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:13:56 2022

@author: michaelselby
"""
import sys
import miek_tools
import pdb_miek



crd = sys.argv[1]
top = sys.argv[2]
xyz = sys.argv[3]
circ = sys.argv[4]

if int(circ) == 1:
    circular = True
else:
    circular = False


traj = pdb_miek.trajectory(crd,circular=circular)
res = traj.load_topology(top)
traj.process_configurations([("molecular contour", miek_tools.mol_cont)])

with open(xyz,"w") as file:
    list_of_quants = traj.quant_dict["molecular contour"]
    file.write(str(int(traj.nres/traj.nstrand))+ "\n")
    file.write("\n")
    for long_array in list_of_quants:
        for r1 in long_array:
            line = "C "+ str(r1[0])+" "+str(r1[1])+" " +str(r1[2])+ "\n"        
            file.write(line)
        file.write("\n")
        file.write("\n")