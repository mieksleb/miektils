#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 12:24:06 2022

Collection of tools for calculating quantities using pdb_miek classes such as conf, strand and base

@author: michaelselby
"""
import numpy as np
from tools import get_angle
from molecular_contour import HAXIS, CAXIS, get_twist, setZ
from numba import njit


def bending_angle(conf):
    strand1, strand2 = conf.strand_list[0:2]
    centres = conf.get_duplex_centres()
    DT = strand1.base_list
    DA = strand2.base_list
    tt1  = DT[19]
    tt2  = DT[20]
    a1 = DA[19]
    a2 = DA[20]
    
    r1 = centres[:18,:]
    r2 = centres[22:,:]
    theta = get_angle(r1, r2) 
    if theta > 90:
        theta = 180-theta

        
    return theta
        

def mol_cont(file, conf, buffer=12):
    """
    Function is called during pdb_miek.trajectory.process_configurations()
    to write additonal lines to molecular contour file
    
    """
    strandA = conf.strand_list[0]
    strandApos = np.array(strandA.get_atom_list("C1'"))
    strandB = conf.strand_list[1]
    strandBpos = strandB.get_atom_list("C1'")
    strandBpos.reverse()
    strandBpos = np.array(strandBpos)

    bp = len(strandApos)

    centres = (strandApos + strandBpos)/2
    diff = - strandApos + strandBpos
    diff /= np.sqrt((diff ** 2).sum(-1))[..., np.newaxis]

    new_centres = np.array([(centres[(i+1)%bp,:]+centres[i,:])/2 for i in range(bp)])

    rC = HAXIS(bp, new_centres, strandApos)

    zvals = np.array([rC[(j+1)%bp,:] - rC[(j-1)%bp,:] for j in range(bp)])
    zvals /= np.sqrt((zvals ** 2).sum(-1))[..., np.newaxis]
    
    
    # newdiff = np.array([np.dot( setZ(zvals),diff ) for z, diff in zip(zvals,diff)])

    
    theta_vals = np.array([get_twist(diff[i,:],diff[(i+1)%bp,:],zvals[i]) for i in range(bp) ])
    
    r1 = CAXIS(bp,centres,theta_vals)
    
    
    if strandA.circular==False: 
        r1 = r1[buffer:bp-buffer,:]


    for i in range(np.shape(r1)[0]):
        line = "C "+ str(r1[i,0])+" "+str(r1[i,1])+" " +str(r1[i,2])+ "\n"        
        file.write(line)
        
        