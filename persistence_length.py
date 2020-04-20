#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 12:03:09 2020

This script computes the persistence length of a spline

@author: michaelselby
"""
import os
import base
import readers
import tools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy
from scipy import optimize
from itertools import chain
import shutil



parent_dir = os.path.dirname(os.path.realpath(__file__))
bin_dir = "bin/"
path = os.path.join(parent_dir, bin_dir) 

try:
    os.mkdir(path)
except Exception:
    pass


traj_file = "miekbin/pl_trajectory.dat"
conf = "miekbin/init.conf"
top = "miekbin/init.top"
npoints = 200



l = readers.LorenzoReader(conf, top)
s = l.get_system()
strand1 = s._strands[0]
strand2 = s._strands[1]
bp = len(strand1._nucleotides)

spline1 = tools.get_bb_spline(strand1, per = False)
spline2 = tools.get_bb_spline(strand2, reverse = True, per = False) 
    
s1xx, s1yy, s1zz = spline1[0]
s2xx, s2yy, s2zz = spline2[0]

# bpi is the base pair index parameter that is common to both splines
bpi = np.linspace(spline1[1][0], spline1[1][1], npoints)

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

tools.plot_splines([spline1, spline2], [], bpi, marker_size = 1)
plt.show()

def tangent_correlation(spline1, spline2, npoints, length):
          
    s1xx, s1yy, s1zz = spline1[0]
    s2xx, s2yy, s2zz = spline2[0]
    
    # bpi is the base pair index parameter that is common to both splines
    bpi = np.linspace(spline1[1][0], spline1[1][1], npoints)
    
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
    
    
    # get smooth splines as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, m1xx, k = 3, s = 0, per = False)
    msyy = scipy.interpolate.splrep(contour_len, m1yy, k = 3, s = 0, per = False)
    mszz = scipy.interpolate.splrep(contour_len, m1zz, k = 3, s = 0, per = False) 
    
    
    final_spline = [[msxx, msyy, mszz], range(len(ss))]
    
    
    tan0 = [(scipy.interpolate.splev(0, final_spline[0][0], 1)), scipy.interpolate.splev(0, final_spline[0][1], 1), scipy.interpolate.splev(0, final_spline[0][2], 1)]
    tan0 = tan0/np.linalg.norm(tan0)
    
    tan = [(scipy.interpolate.splev(length, final_spline[0][0], 1)), scipy.interpolate.splev(length, final_spline[0][1], 1), scipy.interpolate.splev(length, final_spline[0][2], 1)]
    tan = tan/np.linalg.norm(tan)
    corr =  np.dot(tan0, tan)
    return corr


steps = 10000
printed_steps = 10


'''
Converts a trajectory file into mutiple configuration files so that they can 
be passed through various routines
'''
def traj_2_confs(trajectory, bp, energy_out=True):
    list_of_confs = []
    if energy_out == True:
        lines_per_file = 2*bp+3
    else:
        lines_per_file = 2*bp
    smallfile = None
    with open(trajectory) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                if smallfile:
                    smallfile.close()
                small_filename = 'conf'+str(int(lineno/lines_per_file))
                smallfile = open(os.path.join(path, small_filename), "w")
            smallfile.write(line)
        list_of_confs.append(smallfile)
        if smallfile:
            smallfile.close()
            
        return list_of_confs

traj_list = traj_2_confs(traj_file, bp, energy_out=True)

corr_vals = []
for length in ss:
    corr = 0
    for i in range(0, printed_steps):
        conf = os.path.join(path, "conf"+str(i*int(steps/printed_steps)))
        top = top  
        l = readers.LorenzoReader(conf, top)
        s = l.get_system()  
        strand1 = s._strands[0]
        strand2 = s._strands[1]
        spline1 = tools.get_bb_spline(strand1, per = False)
        spline2 = tools.get_bb_spline(strand2, reverse = True, per = False) 
        corr += tangent_correlation(spline1, spline2, npoints, length)
        
    corr = corr/printed_steps
    corr_vals.append(corr)
    

a = scipy.optimize.curve_fit(lambda t,a: np.exp(-a*t),  ss,  corr_vals,  p0=(1/100))
P = 1/a[0]
print(1/a[0])

plt.plot(ss, corr_vals)
plt.plot(ss, np.exp(-ss/P))
plt.show()

def s_to_pos(spline, s, ss):
    x = scipy.interpolate.splev(ss, spline[0][0])
    y = scipy.interpolate.splev(ss, spline[0][1])
    z = scipy.interpolate.splev(ss, spline[0][2])
    s_index = min(range(len(ss)), key=lambda i: abs(ss[i]-s))
    pos = [x[s_index], y[s_index], z[s_index]]
    return pos, s_index



    
    
shutil.rmtree(path)