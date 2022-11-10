#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 12:24:06 2022

Collection of tools for calculating quantities using pdb_miek classes such as conf, strand and base

@author: michaelselby
"""
import numpy as np
from tools import get_angle
from molecular_contour import get_molecular_contour
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
        

def mol_cont(conf, buffer=5):
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

    r1 = get_molecular_contour(bp, strandApos, strandBpos, circular=conf.circular)
    
    
    if conf.circular==False: 
        r1 = r1[buffer:bp-buffer,:]
        
    return r1


        
        