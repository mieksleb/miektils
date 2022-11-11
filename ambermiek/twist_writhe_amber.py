#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 09:42:19 2022

@author: michaelselby
"""

import sys


crd = sys.argv[1]
top = sys.argv[2]
circ = sys.argv[3]
twist_writhe_file = sys.argv[4]


import pdb_miek
from miek_tools import get_twist_writhe_conf


if int(circ) == 1:
    circular = True
else:
    circular = False


traj = pdb_miek.trajectory(crd,circular=circular)
res = traj.load_topology(top)
traj.process_configurations([("twist writhe", get_twist_writhe_conf)])


list_of_quants = traj.quant_dict["twist writhe"]

with open(twist_writhe_file, "w") as file:
    for twist, writhe in list_of_quants:
        file.write(str(twist) +" " +str(writhe)+ '\n')

    
        
        

