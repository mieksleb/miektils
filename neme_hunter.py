#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:13:12 2020

'neme hunter 2k20, this bad boi hunts down 'nemes with an ungodly thirst for 'neme blood

@author: michaelselby
"""
import base
import readers
import origami_utils
import tools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interpolate
import scipy
from itertools import chain

#traj_file = "miekbin/trajectory_300-Lk7.dat"
#conf = "miekbin/last_conf.dat"
#top = "miekbin/circle_300-Lk7.top"

plect_file = open("miekbin/plectoneme624.dat","r").readlines()

conf = "miekbin/last_conf624_2_(plect).dat"
#conf = "miekbin/last_conf624_2.dat"
top = "miekbin/sim624.top"
#
l = readers.LorenzoReader(conf, top)
s = l.get_system()
#
line = plect_file[0].split()
neme_guess = line[3]

strand1 = s._strands[0]
strand2 = s._strands[1]
#
bp = len(strand1._nucleotides)
#
spline1 = tools.get_bb_spline(strand1)
spline2 = tools.get_bb_spline(strand2, reverse = True)  


'''
Given the two splines which run through the base positions in a duplex, neme hunter will return
a list of plectonemes called 'nemes'. For circular = True, it assumes there are two nemes located
at approximately the other side of the minicircle
'''

def neme_hunter(spline1, spline2, npoints, bp, guess, circular = True, k=0.5):
          
    s1xx, s1yy, s1zz = spline1[0]
    s2xx, s2yy, s2zz = spline2[0]
    
    # bpi is the base pair index parameter that is common to both splines
    bpi = np.linspace(spline1[1][0], spline1[1][1], npoints)
    
    # find the midpoint between the input splines, as a function of base pair index
    m1xx = (scipy.interpolate.splev(bpi, s1xx) + scipy.interpolate.splev(bpi, s2xx)) / 2
    m1yy = (scipy.interpolate.splev(bpi, s1yy) + scipy.interpolate.splev(bpi, s2yy)) / 2
    m1zz = (scipy.interpolate.splev(bpi, s1zz) + scipy.interpolate.splev(bpi, s2zz)) / 2
    
   
    dogs = [[m1xx[i], m1yy[i], m1zz[i]] for i in range(len(m1xx))]   #list of coordinates in bpi
    
    
    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((m1xx[ii+1]-m1xx[ii])**2+(m1yy[ii+1]-m1yy[ii])**2+(m1zz[ii+1]-m1zz[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)
    
    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)
    
    
    # get smooth splines as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, m1xx, k = 3, s = 10, per = True)
    msyy = scipy.interpolate.splrep(contour_len, m1yy, k = 3, s = 10, per = True)
    mszz = scipy.interpolate.splrep(contour_len, m1zz, k = 3, s = 10, per = True) 
    
    
    final_spline = [[msxx, msyy, mszz], range(len(ss))]
    
    
    xx = scipy.interpolate.splev(ss, msxx)
    yy = scipy.interpolate.splev(ss, msyy)
    zz = scipy.interpolate.splev(ss, mszz)
    
    cats = [[xx[i], yy[i], zz[i]] for i in range(len(xx))]
    
    
    nxx = (scipy.interpolate.splev(ss, final_spline[0][0], 2)) 
    nyy = (scipy.interpolate.splev(ss, final_spline[0][1], 2)) 
    nzz = (scipy.interpolate.splev(ss, final_spline[0][2], 2)) 
    

    curvs = []       
    for i in range(len(ss)):
        curvs.append(nxx[i]*nxx[i]+nyy[i]*nyy[i]+nzz[i]*nzz[i])
        

    
    
    '''
    This converts an arclength to a postion (x,y,z)
    '''
    def s_to_pos(spline, s, ss):
        x = scipy.interpolate.splev(ss, spline[0][0])
        y = scipy.interpolate.splev(ss, spline[0][1])
        z = scipy.interpolate.splev(ss, spline[0][2])
        s_index = min(range(len(ss)), key=lambda i: abs(ss[i]-s))
        pos = [x[s_index], y[s_index], z[s_index]]
        return pos, s_index
    
    def pos_to_s(spline, pos, ss):
        x = scipy.interpolate.splev(ss, spline[0][0])
        y = scipy.interpolate.splev(ss, spline[0][1])
        z = scipy.interpolate.splev(ss, spline[0][2])
        
        s_pos = [[msxx[i], msyy[i], mszz[i]] for i in range(len(msxx))] 
        
        s_index = min(range(len(ss)), key=lambda i: np.linalg.norm(np.array[pos]-s_pos[i]))
        s = ss[s_index]
        return s
    
    
    
    if circular:
        jsp1 = max(range(len(ss)), key=lambda i: curvs[i])
        
        midzone = s_to_pos(final_spline, (ss[jsp1]+ss[-1]/2) % ss[-1], ss)
        
        if (midzone[1]+10 < npoints):
            region = range(midzone[1]-10, (midzone[1]+10))
        else:
            region = chain(range(0,(midzone[1]+10) % npoints),  range(midzone[1]-10, npoints))
            
        jsp2 = max(region, key=lambda i: curvs[i])
        
        pos1 = s_to_pos(final_spline, ss[jsp1], ss)[0]
        pos2 = s_to_pos(final_spline, ss[jsp2], ss)[0]
          
        boyos = [pos1, pos2]
        
        neme1 = int(bpi[tools.closest_point(pos1, dogs)])
        neme2 = int(bpi[tools.closest_point(pos2, dogs)])
        
        nemes = [neme1, neme2]
        
    else:
        
        a = 5      # cutoff paramter to prevent large unphysical end curvatures from contributing
        ss = ss[a:]
        curvs = curvs[:-a]
        curvs = curvs[a:]
        
        guess_pos = dogs[int(float(guess))]  #have to convert the guess which is a bp parameter to an arclength parameter
        
        
        j = min(range(len(ss)), key=lambda i: np.linalg.norm(np.array(guess_pos)-np.array(cats[i])))

        curvs = [curvs[i]*np.exp(-k*(i-j)**2) for i in range(len(curvs))]
        
        jsp1 = max(range(len(ss)-2*a), key=lambda i: curvs[i])
        pos1 = s_to_pos(final_spline, ss[jsp1], ss)[0]             
        neme1 = int(bpi[tools.closest_point(pos1, dogs)])
        
        nemes = [neme1]
        

    return nemes





ohno = neme_hunter(spline1, spline2, 500, 624, circular = False, guess = neme_guess)
