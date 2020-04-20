#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:13:05 2020

@author: michaelselby
"""
import sys
import base
import numpy as np
import tools
import readers

import matplotlib.pyplot as plt
from numba import njit
from numba import jit
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt
import scipy
import cmath
import numba


# conf = "hannah/round2/LH_0_10.conf"
# #conf = "miekbin/last_conf624_2.dat"
# top = "hannah/round2/LH_0_10.top"


# conf = "miekbin/last_conf.dat"
# # conf = "miekbin/circle_300-Lk7.dat"
# top = "miekbin/circle_300-Lk7.top"

conf = "paktep/conf.conf"
top = "paktep/top.top"


l = readers.LorenzoReader(conf, top)
s = l.get_system()


strand1 = s._strands[0]
strand2 = s._strands[1]

bp = len(strand1._nucleotides)

spline1 = tools.get_bb_spline(strand1)
spline2 = tools.get_bb_spline(strand2, reverse = True)   

def get_twist_writhe(spline1, spline2, npoints = 1000, circular = False, integral_type = "simple"):
    
    """
    return the twist and writhe for a given configuration and 
    Using integral_type = 'simple' 

    args:
    spline1: list of 3 splines corresponding to x, y and z spline through strand 1's backbone
    spline2: list of 3 splines corresponding to x, y and z spline through strand 2's backbone -- NB the splines should run in the same direction, i.e. one must reverse one of the splines if they come from get_base_spline (e.g. use get_base_spline(reverse = True))

    npoints: number of points for the discrete integration
    """

    import scipy.interpolate
    
         

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
    

    xx = scipy.interpolate.splev(ss, msxx)
    yy = scipy.interpolate.splev(ss, msyy)
    zz = scipy.interpolate.splev(ss, mszz)
    
    

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of normalised tangent vectors; __call__(xxx, 1) returns the first derivative
    # the tangent vector is a unit vector
    dmxx = scipy.interpolate.splev(ss, msxx, 1)
    dmyy = scipy.interpolate.splev(ss, msyy, 1)
    dmzz = scipy.interpolate.splev(ss, mszz, 1)

    
    
    tt = list(range(len(ss)))
    for ii in range(len(ss)):
        tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])
  
    
    # get the normal vector via n(s) = dt(s)/ds = d^2/ds^2[r(s)]    
    ddmxx = scipy.interpolate.splev(ss, msxx, 2)
    ddmyy = scipy.interpolate.splev(ss, msyy, 2)
    ddmzz = scipy.interpolate.splev(ss, mszz, 2)
    nn = list(range(len(ss)))
    for ii in range(len(ss)):
        nn[ii] = np.array([ddmxx[ii], ddmyy[ii], ddmzz[ii]])
 

    # we also need the 'normal' vector u(s) which points between the base pairs. (or between the spline fits through the backbones in this case)
    # n.b. these uxx, uyy, uzz are not normalised
    uxx_bpi = scipy.interpolate.splev(bpi, s2xx) - scipy.interpolate.splev(bpi, s1xx)
    uyy_bpi = scipy.interpolate.splev(bpi, s2yy) - scipy.interpolate.splev(bpi, s1yy)
    uzz_bpi = scipy.interpolate.splev(bpi, s2zz) - scipy.interpolate.splev(bpi, s1zz)

    # get the normal vector spline as a function of contour length
    suxx = scipy.interpolate.splrep(contour_len, uxx_bpi, k = 3, s = 0, per = circular)
    suyy = scipy.interpolate.splrep(contour_len, uyy_bpi, k = 3, s = 0, per = circular)
    suzz = scipy.interpolate.splrep(contour_len, uzz_bpi, k = 3, s = 0, per = circular)
    
    dogs = [[xx[i], yy[i], zz[i]] for i in range(len(xx))]

    tools.plot_splines([spline1], [dogs], range(bp), marker_size = 1)
    plt.show()
    
    
    tools.plot_splines([spline1, spline2], [], range(bp), marker_size = 1)
    plt.show()

    # evaluate the normal vector spline as a function of contour length
    uxx = scipy.interpolate.splev(ss, suxx)
    uyy = scipy.interpolate.splev(ss, suyy)
    uzz = scipy.interpolate.splev(ss, suzz)
    


    uu = list(range(len(ss)))
    for ii in list(range(len(ss))):
        uu[ii] = np.array([uxx[ii], uyy[ii], uzz[ii]])
        uu[ii] = uu[ii] - np.dot(tt[ii], uu[ii]) * tt[ii]
        # the normal vector should be normalised
        uu[ii] = tools.norm(uu[ii])
        
    print(xx[1],yy[1],zz[1])
    print(dmxx[1],dmyy[1],dmzz[1])    
        
    print(uxx_bpi[100])
    print(uzz[777])
    print(uu[27])
    print(uu[100])
    # and finally we need the derivatives of that vector u(s). It takes a bit of work to get a spline of the normalised version of u from the unnormalised one
    nuxx = [vec[0] for vec in uu]
    nuyy = [vec[1] for vec in uu]
    nuzz = [vec[2] for vec in uu]
    nusxx = scipy.interpolate.splrep(ss, nuxx, k = 3, s = 0, per = circular)
    nusyy = scipy.interpolate.splrep(ss, nuyy, k = 3, s = 0, per = circular)
    nuszz = scipy.interpolate.splrep(ss, nuzz, k = 3, s = 0, per = circular)
    duxx = scipy.interpolate.splev(ss, nusxx, 1)
    duyy = scipy.interpolate.splev(ss, nusyy, 1)
    duzz = scipy.interpolate.splev(ss, nuszz, 1)
    duu = list(range(len(ss)))
    for ii in list(range(len(ss))):
        duu[ii] = np.array([duxx[ii], duyy[ii], duzz[ii]])

    
    ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
    # do the integration w.r.t. s
 
    
    # triple vector product using numba fastmast
    @njit(fastmath=True)
    def multidet(u,v,w):
        d=np.empty(3)
        d=\
        u[0]*(v[1]*w[2]-v[2]*w[1])+\
        u[1]*(v[2]*w[0]-v[0]*w[2])+\
        u[2]*(v[0]*w[1]-v[1]*w[0])  # 14 operations / det
        return d

    
    
    # now for the integration which we use numba fastmath capabilities for speed
    @njit(fastmath=True)
    def discrete_dbl_int(tt, uu, duu, xx, yy, zz, dx, dy, ss, circular = circular):
        if circular:
            srange = range(len(ss)-1)
        else:
            srange = range(len(ss))
        twist_integral = 0
        writhe_integral = 0
        for ii in srange:
            triple_scalar_product = multidet(tt[ii], uu[ii], duu[ii])
            twist_integral += triple_scalar_product 
            for jj in srange:
                # skip ii=jj and use symmetry in {ii, jj}
                if ii > jj:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = multidet((tt[ii]), tt[jj], diff_frac)
                    writhe_integral += triple_scalar_product
            
        twist_integral *= dx / (2 * np.pi)
        writhe_integral *= dx * dy / (2 * np.pi) # multiplied by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
        return twist_integral, writhe_integral
    
    
    twist, writhe = discrete_dbl_int(tt, uu, duu, xx, yy, zz, ds, ds, ss, circular = True)
    
    return twist, writhe
    
    
 
        
#        
#l = readers.LorenzoReader(conf, top)
#s = l.get_system()
#
#
#strand1 = s._strands[0]
#strand2 = s._strands[1]
#
#bp = len(strand1._nucleotides)
#
#spline1 = tools.get_bb_spline(strand1)
#spline2 = tools.get_bb_spline(strand2, reverse = True)        
#        
#        
#snoot = hoot = spline1[0]
#hoot = spline1[0][0]      
#
#
#h = 10.67


twist, writhe = get_twist_writhe(spline1, spline2)
boyo = twist + writhe 
#
#twisty = twist - 656/h
#
#boyo = twisty + writhe
#
#print(twisty, "twist")
#print(writhe, "writhe")
#print(boyo)