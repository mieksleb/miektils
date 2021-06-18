#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 14:15:51 2021

All tools which require base or readers from oxDNA untils

@author: michaelselby
"""

import base
import readers


def get_base_pos(strand, reverse = False):
    base_pos = []
    for nuc in strand._nucleotides:
        base_pos.append(base.Nucleotide.get_pos_base(nuc))
    if reverse:
        base_pos.reverse()
    if strand._circular:
        if reverse:
            base_pos.append(strand._nucleotides[-1].get_pos_base())
        else:
            base_pos.append(strand._nucleotides[0].get_pos_base())
    return base_pos


def get_bb_pos(strand, reverse = False):
    bb_pos = []
    for nuc in strand._nucleotides:
        bb_pos.append(base.Nucleotide.get_pos_back(nuc))
    if reverse:
        bb_pos.reverse()
    if strand._circular:
        if reverse:
            bb_pos.append(strand._nucleotides[-1].get_pos_back())
        else:
            bb_pos.append(strand._nucleotides[0].get_pos_back())
    return bb_pos

  



"""
return a cartesian spline that represents a fit through the bases for the strand 'strand'

args:
strand: base.Strand object
"""

def get_base_spline(strand, reverse = False):

    import scipy

    base_pos = []
    for nuc in strand._nucleotides:
        base_pos.append(nuc.get_pos_base())

    if reverse:
        base_pos.reverse()

    if strand._circular:
        if reverse:
            base_pos.append(strand._nucleotides[-1].get_pos_base())
        else:
            base_pos.append(strand._nucleotides[0].get_pos_base())

    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in base_pos]
    yy = [vec[1] for vec in base_pos]
    zz = [vec[2] for vec in base_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = strand._circular)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = strand._circular)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = strand._circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)




def get_bb_spline(strand, reverse = False, per = True):
    """
    return a cartesian spline that represents a fit through the backbone of a duplex

    args:
    strand: base.Strand object
    """

    from scipy.interpolate import splev, splrep

 
    bb_pos = []
    for nuc in strand._nucleotides:
        bb_pos.append(nuc.get_pos_back())
    if reverse:
        bb_pos.reverse()
    if strand._circular:
        if reverse:
            bb_pos.append(strand._nucleotides[-1].get_pos_back())
        else:
            bb_pos.append(strand._nucleotides[0].get_pos_back())
 
    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in bb_pos]
    yy = [vec[1] for vec in bb_pos]
    zz = [vec[2] for vec in bb_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = splrep(range(len(xx)), xx, k = 3, s = 0, per = strand._circular)
    spline_yy = splrep(range(len(yy)), yy, k = 3, s = 0, per = strand._circular)
    spline_zz = splrep(range(len(zz)), zz, k = 3, s = 0, per = strand._circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)
