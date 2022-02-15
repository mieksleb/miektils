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

pdb = "minicircle/Hussain/dv160t0tt5_md9_nw.pdb"

import pdb_miek
import tools

traj = pdb_miek.trajectory(pdb,circular=True)
# steps = traj.steps
traj.add_configurations()

npoints = 1000
circular = True

twist_writhe_file = 'minicircle/Hussain/twist_writhe2.dat'
file = open(twist_writhe_file,"w")

twist_list = []
writhe_list = []

# main loop is over timesteps

double = True
old_writhe = 0

step = 1
for conf in traj.config_list:
    strand1 = conf.strand_list[0]
    strand1.get_bb_list()
    strand1pos = strand1.bb_list
    
    strand2 = conf.strand_list[1]
    strand2.get_bb_list()
    strand2pos = strand2.bb_list
    spline1 = tools.get_spline(strand1pos, per=circular,reverse=True)
    spline2 = tools.get_spline(strand2pos,per=circular,reverse=False)
     
    if double==True:
        twist, writhe = tools.get_twist_writhe(spline1, spline2, npoints=1000)
        print(twist)
        print(writhe)
        print(twist + writhe)
        Lk = round(twist+writhe)
        print("\n")
        # double=False
        
        
    else:
            
        twist, writhe = tools.get_twist_writhe_fuller(spline1, spline2, old_writhe, t, n, npoints = npoints, circular=circular)
        
        print("Fuller time")
        print(twist)
        print(writhe)
        print(twist + writhe)
        print("\n")
        
        # twist, writhe = tools.get_twist_writhe(spline1, spline2, npoints=1000)
        # print(twist)
        # print(writhe)
        # print(twist + writhe)
        # print("\n")
        
        # if abs(twist+writhe-Lk) > 0.02:
        #     print("Fuller insufficient")
        #     twist, writhe = tools.get_twist_writhe(spline1, spline2, npoints=1000)
        #     print(twist)
        #     print(writhe)
        #     print(twist + writhe)
        #     print("\n")
            
    file.write(str(step) +' '+ str(twist) +' '+str(writhe)+ '\n')
    # file.write(step,twist,writhe)
    step += 1
    twist_list.append(twist)
    writhe_list.append(writhe)
    old_writhe = writhe
    s1xx, s1yy, s1zz = spline1[0]
    s2xx, s2yy, s2zz = spline2[0]
    
    smin = spline1[1][0]
    smax = spline1[1][-1]

    # bpi is the base pair index parameter that is common to both splines
    bpi = np.linspace(smin, smax, npoints)

    # find the midpoint between the input splines, as a function of base pair index
    m1xx = (scipy.interpolate.splev(bpi, s1xx) + scipy.interpolate.splev(bpi, s2xx)) / 2
    m1yy = (scipy.interpolate.splev(bpi, s1yy) + scipy.interpolate.splev(bpi, s2yy)) / 2
    m1zz = (scipy.interpolate.splev(bpi, s1zz) + scipy.interpolate.splev(bpi, s2zz)) / 2
    

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((m1xx[ii+1]-m1xx[ii])**2+(m1yy[ii+1]-m1yy[ii])**2+(m1zz[ii+1]-m1zz[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    # get the splines as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, m1xx, k = 3, s = 0, per = circular)
    msyy = scipy.interpolate.splrep(contour_len, m1yy, k = 3, s = 0, per = circular)
    mszz = scipy.interpolate.splrep(contour_len, m1zz, k = 3, s = 0, per = circular)

    dmxx = scipy.interpolate.splev(ss, msxx, 1)
    dmyy = scipy.interpolate.splev(ss, msyy, 1)
    dmzz = scipy.interpolate.splev(ss, mszz, 1)
    tt = list(range(len(ss)))
    for ii in range(len(ss)):
        tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])
        

    def t(s):
        ddmxx = scipy.interpolate.splev(s, msxx, 1)
        ddmyy = scipy.interpolate.splev(s, msyy, 1)
        ddmzz = scipy.interpolate.splev(s, mszz, 1)
        tan = np.array([ddmxx, ddmyy, ddmzz])
        tan /= np.linalg.norm(tan)
        return tan
  
    
    # get the normal vector via n(s) = dt(s)/ds = d^2/ds^2[r(s)]    
    ddmxx = scipy.interpolate.splev(ss, msxx, 2)
    ddmyy = scipy.interpolate.splev(ss, msyy, 2)
    ddmzz = scipy.interpolate.splev(ss, mszz, 2)
    nn = list(range(len(ss)))
    for ii in range(len(ss)):
        nn[ii] = np.array([ddmxx[ii], ddmyy[ii], ddmzz[ii]])
        
    def n(s):
        ddmxx = scipy.interpolate.splev(s, msxx, 2)
        ddmyy = scipy.interpolate.splev(s, msyy, 2)
        ddmzz = scipy.interpolate.splev(s, mszz, 2)
        return np.array([ddmxx, ddmyy, ddmzz])
    
    
    


    
file.close()
        
        

