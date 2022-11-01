#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:13:56 2022

@author: michaelselby
"""
import molecular_contour
import pdb_miek

# crd = "minicircle/Hussain/MDNA/dv160t0_md9_nw.trj"
# # crd = "minicircle/Hussain/dv160t0tt5_long_nw.trj"
# # crd = "TT_40/MDNA/dv40t_longrun_nw.trj"
# top = "minicircle/Hussain/MDNA/dv160t0.top"


# xyzfile = "minicircle/Hussain/molecular_contour_MTTD5.xyz"
# xyzfile = "minicircle/Hussain/molecular_contour_MDNA.xyz"
# xyzfile = "TT_40/MTTD/molecular_contour_MTTD.xyz"
# xyzfile = "TT_40/MDNA/molecular_contour.xyz"


crd = "minicircle/Hussain/MDNA/dv160t0_md9_nw.trj"
top = "minicircle/Hussain/MDNA/dv160t0.top"
xyzfile = "minicircle/Hussain/molecular_contour_MDNA.xyz"
traj = pdb_miek.trajectory(crd,circular=True)
res = traj.load_topology(top)
traj.process_configurations(molecular_contour=True, mol_cont_file=xyzfile)
