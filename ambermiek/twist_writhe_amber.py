#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 09:42:19 2022

@author: michaelselby
"""

import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skspatial.objects import Line
from skspatial.objects import Points

sys.path.append('/Users/michaelselby/Documents/DPhil/miektils/ambermiek')
sys.path.append('/Users/michaelselby/Documents/DPhil/miektils')


crd = "minicircle/Hussain/MDNA/dv160t0_md9_nw.trj"


import pdb_miek
import tools



pdb = "/Users/michaelselby/Documents/DPhil/figures/trefoil/conf.dat.pdb"
traj = pdb_miek.trajectory(pdb,circular=True)
# steps = traj.steps
# traj.add_configurations()


# traj = pdb_miek.trajectory(crd,circular=True)
# seq = traj.load_sequence("40AT.seq")
# res = traj.load_topology("minicircle/Hussain/MDNA/dv160t0.top")
# traj.add_configurations()
# nsteps = traj.nsteps

circular = True

# twist_writhe_file = 'minicircle/Hussain/MDNA/twist_writhe_MDNA.dat'
# file = open(twist_writhe_file,"w")

twist_list = []
writhe_list = []

# main loop is over timesteps

traj.process_configurations(twistwrithe = True)
writhe_list = traj.writhe_list
twist_list = traj.twist_list

# step = 1
# for twist, writhe in zip(twist_list,writhe_list):
#     file.write(str(step) +' '+ str(twist) +' '+str(writhe)+ '\n')
#     step += 1

# step = 1
# for conf in traj.config_list:
#     strand1 = conf.strand_list[0]
#     strand1.get_bb_list()
#     strand1pos = strand1.bb_list
    
#     strand2 = conf.strand_list[1]
#     strand2.get_bb_list()
#     strand2pos = strand2.bb_list
#     spline1 = tools.get_spline(strand1pos, per=circular,reverse=True)
#     spline2 = tools.get_spline(strand2pos,per=circular,reverse=False)
     
#     twist, writhe = tools.get_twist_writhe(spline1, spline2, npoints=2000)
#     print(twist)
#     print(writhe)
#     print(step)
#     # print(twist)
#     # print(writhe)
#     print("\n")
                 
#     file.write(str(step) +' '+ str(twist) +' '+str(writhe)+ '\n')
    

#     step += 1

    
    
    
# file.close()
        
        

