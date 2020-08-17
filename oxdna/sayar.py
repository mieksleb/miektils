#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:43:38 2019

@author: michaelselby
"""

import base
try:
    import numpy as np
except:
    import mynumpy as np
import sys
import subprocess
import pickle
import os
import tempfile

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

def get_pos_midpoint(r1, r2, box):
    """
    Returns the midpoint of two vectors r1 and r2.
    
    Uses the minimum image: i.e. make sure we get a sensible answer if the                                                                     
    positions r1 and r2 correspond to nucleotides which were put at opposite                                                                  
    ends of the box due to periodic boundary conditions.
    Args:
        r1: first vector
        r2: second vector
        box: box dimensions in a numpy array                                                                                            
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    return r1 - min_distance(r1, r2, box)/2

def array_index(mylist, myarray):
    """
    Returns the index in a list mylist of a numpy array myarray.
    Args:
        mylist: list
        myarray: numpy array
    """
    return map(lambda x: (myarray == x).all(), mylist).index(True)

def min_distance (r1, r2, box):
    """
    Returns the minimum image distance in going from r1 to r2, in a box of size box.
    Same as base.Nucleotide.distance().
    Args:
        r1: first vector
        r2: second vector
        box: box dimensions in a numpy array
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    dr = r2 - r1
    dr -= box * np.rint (dr / box)
    return dr

def vecs2spline(vecs, per):
    """
    This function was not documented when it was written.
    The author of this documentation believes that it returns a spline interpolating the given vectors.
    Args:
        vecs: list of vectors (arrays of length 3)
        per: boolean. If true, forces spline to be periodic.
    """
    import scipy.interpolate
    # interpolate vecs by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in vecs]
    yy = [vec[1] for vec in vecs]
    zz = [vec[2] for vec in vecs]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = per)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = per)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = per)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def get_base_spline(strand, reverse = False):
    """
    Returns a cartesian spline that represents a fit through the bases for the strand 'strand'.
    Args:
        strand: base.Strand object
        reverse: boolean. If true, strand is reversed
    """
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

def get_sayar_twist(s1, s2, smin, smax, npoints = 1000, circular = False, integral_type = "simple"):
    """
    Returns the twist for a given pair of spline fits, one through the bases of each strand. 
    From Sayar et al. Phys. Rev. E, 81, 041916 (2010)
    Just need to integrate along the contour parameter s that is common to both splines. 
    We need the normalised tangent vector to the spline formed by 
    the midpoint of the two input splines t(s), 
    the normalised normal vector formed by the vectors between the splines u(s), 
    and the derivative of the normalised normal vector between the splines d/ds (u(s)).
    NB, the normal u(s) vector should be orthogonal to the tangent vector t(s); we ensure this by using only the component orthogonal to t(s).
    Using integral_type = 'simple' and npoints = 200, it will give a correct twist, or at least one that gives a conserved linking number when combined with get_sayar_writhe.
    Args:
        s1: list of 3 splines corresponding to 3-D spline through strand 1's bases (e.g. use get_base_spline())
        s2: list of 3 splines corresponding to 3-D spline through strand 2's bases -- NB the splines should run in the same direction, i.e. one must reverse one of the splines if they come from get_base_spline (e.g. use get_base_spline(reverse = True))
        smin: minimum value for s, which parameterises the splines
        smax: maximum value for s, which parameterises the splines
        npoints: number of points for the discrete integration
        circular: boolean. If true, spline is circular.
        integral_type: "simple"
    """

    import scipy.interpolate
    import scipy.integrate
    
    s1xx, s1yy, s1zz = s1
    s2xx, s2yy, s2zz = s2

    # bpi is the base pair index parameter that common to both splines
    bpi = np.linspace(smin, smax, npoints)

    # find the midpoint between the input splines, as a function of base pair index
    mxx = (scipy.interpolate.splev(bpi, s1xx) + scipy.interpolate.splev(bpi, s2xx)) / 2
    myy = (scipy.interpolate.splev(bpi, s1yy) + scipy.interpolate.splev(bpi, s2yy)) / 2
    mzz = (scipy.interpolate.splev(bpi, s1zz) + scipy.interpolate.splev(bpi, s2zz)) / 2

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((mxx[ii+1]-mxx[ii])**2+(myy[ii+1]-myy[ii])**2+(mzz[ii+1]-mzz[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    # get the midpoint spline as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, mxx, k = 3, s = 0, per = circular)
    msyy = scipy.interpolate.splrep(contour_len, myy, k = 3, s = 0, per = circular)
    mszz = scipy.interpolate.splrep(contour_len, mzz, k = 3, s = 0, per = circular)

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of normalised tangent vectors; __call__(xxx, 1) returns the first derivative
    # the tangent vector is a unit vector
    dmxx = scipy.interpolate.splev(ss, msxx, 1)
    dmyy = scipy.interpolate.splev(ss, msyy, 1)
    dmzz = scipy.interpolate.splev(ss, mszz, 1)
    tt = range(len(ss))
    for ii in range(len(ss)):
        tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])

    # we also need the 'normal' vector u(s) which points between the base pairs. (or between the spline fits through the bases in this case)
    # n.b. these uxx, uyy, uzz are not normalised
    uxx_bpi = scipy.interpolate.splev(bpi, s2xx) - scipy.interpolate.splev(bpi, s1xx)
    uyy_bpi = scipy.interpolate.splev(bpi, s2yy) - scipy.interpolate.splev(bpi, s1yy)
    uzz_bpi = scipy.interpolate.splev(bpi, s2zz) - scipy.interpolate.splev(bpi, s1zz)

    # get the normal vector spline as a function of contour length
    suxx = scipy.interpolate.splrep(contour_len, uxx_bpi, k = 3, s = 0, per = circular)
    suyy = scipy.interpolate.splrep(contour_len, uyy_bpi, k = 3, s = 0, per = circular)
    suzz = scipy.interpolate.splrep(contour_len, uzz_bpi, k = 3, s = 0, per = circular)

    # evaluate the normal vector spline as a function of contour length
    uxx = scipy.interpolate.splev(ss, suxx)
    uyy = scipy.interpolate.splev(ss, suyy)
    uzz = scipy.interpolate.splev(ss, suzz)

    uu = range(len(ss))
    for ii in range(len(ss)):
        uu[ii] = np.array([uxx[ii], uyy[ii], uzz[ii]])
        uu[ii] = uu[ii] - np.dot(tt[ii], uu[ii]) * tt[ii]
        # the normal vector should be normalised
        uu[ii] = norm(uu[ii])

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
    duu = range(len(ss))
    for ii in range(len(ss)):
        duu[ii] = np.array([duxx[ii], duyy[ii], duzz[ii]])

    ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
    # do the integration w.r.t. s
    if circular:
        srange = range(len(ss)-1)
    else:
        srange = range(len(ss))
    if integral_type == "simple":
        integral = 0
        for ii in srange:
            #print np.dot(uu[ii], tt[ii])
            triple_scalar_product = np.dot(tt[ii], np.cross(uu[ii], duu[ii]))
            integral += triple_scalar_product * ds
    elif integral_type == "quad":
        assert False, "not currently supported; shouldn't be difficult to implement if wanted"
        integral, err = scipy.integrate.quad(twist_integrand, ss[0], ss[-1], args = (msxx, msyy, mszz, nusxx, nusyy, nuszz), limit = 500)
        print >> sys.stderr, "error estimate:", err
        
    twist = integral/(2 * np.pi)

    return twist

def get_sayar_writhe(splines1, smin, smax, splines2 = False, npoints = 1000, debug = False, circular = False, integral_type = "simple"):
    """
    Returns the writhe for a 3D spline fit through a set of duplex midpoints.
    From Sayar et al. Phys. Rev. E, 81, 041916 (2010).
    Using integral_type = 'simple' and npoints = 200, it will give a correct writhe, or at least one that gives a conserved linking number when combined with get_sayar_twist.
    Args:
        splines1: list of 3 splines corresponding to either (if not splines2) a 3D spline through the duplex or (if splines2) strand 1's bases
        smin: minimum value for s, which parameterises the splines
        smax: maximum value for s, which parameterises the splines
        splines2: optionally, (see splines1) list of 3 splines corresponding to a 3D spline through strand2's bases
        npoints: number of points for the discrete integration
        debug: print a load of debugging information
        circular: boolean. If true, spline is circular.
        integral_type: "simple"
    """
    import scipy.integrate

    # bpi is the base pair index parameter that common to both strands' splines
    bpi = np.linspace(smin, smax, npoints)

    ## get midpoint splines sxx, syy, szz
    if not splines2:
        # splines1 is the midpoint 3D spline as a function of base pair index
        sxx_bpi, syy_bpi, szz_bpi = splines1
        xx_bpi = scipy.interpolate.splev(bpi, sxx_bpi)
        yy_bpi = scipy.interpolate.splev(bpi, syy_bpi)
        zz_bpi = scipy.interpolate.splev(bpi, szz_bpi)
    else:
        # take splines1 and splines2 to be the splines through the bases of each strand; in that case we need to find the midpoint here first
        s1xx_bpi, s1yy_bpi, s1zz_bpi = splines1
        s2xx_bpi, s2yy_bpi, s2zz_bpi = splines2

        # find the midpoint as a function of base pair index between the input splines 
        xx_bpi = (scipy.interpolate.splev(bpi, s1xx_bpi) + scipy.interpolate.splev(bpi, s2xx_bpi)) / 2
        yy_bpi = (scipy.interpolate.splev(bpi, s1yy_bpi) + scipy.interpolate.splev(bpi, s2yy_bpi)) / 2
        zz_bpi = (scipy.interpolate.splev(bpi, s1zz_bpi) + scipy.interpolate.splev(bpi, s2zz_bpi)) / 2

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((xx_bpi[ii+1]-xx_bpi[ii])**2+(yy_bpi[ii+1]-yy_bpi[ii])**2+(zz_bpi[ii+1]-zz_bpi[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    sxx = scipy.interpolate.splrep(contour_len, xx_bpi, k = 3, s = 0, per = circular)
    syy = scipy.interpolate.splrep(contour_len, yy_bpi, k = 3, s = 0, per = circular)
    szz = scipy.interpolate.splrep(contour_len, zz_bpi, k = 3, s = 0, per = circular)
    xx = scipy.interpolate.splev(ss, sxx)
    yy = scipy.interpolate.splev(ss, syy)
    zz = scipy.interpolate.splev(ss, szz)

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of tangent vectors; __call__(xxx, 1) returns the first derivative
    dxx = scipy.interpolate.splev(ss, sxx, 1)
    dyy = scipy.interpolate.splev(ss, syy, 1)
    dzz = scipy.interpolate.splev(ss, szz, 1)
    tt = range(len(ss))
    for ii in range(len(ss)):
        tt[ii] = np.array([dxx[ii], dyy[ii], dzz[ii]])
    
    # do the double integration w.r.t. s and s'
    if integral_type == "simple":
        integral = 0
        if circular:
            srange = range(len(ss)-1)
            ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
        else:
            srange = range(len(ss))
            ds = float(contour_len[-1] - contour_len[0]) / npoints
        for ii in srange:
            for jj in srange:
                # skip ii=jj and use symmetry in {ii, jj}
                if ii > jj:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = np.dot(np.cross(tt[ii], tt[jj]), diff_frac)
                    integral += triple_scalar_product * ds * ds
        # multiply by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
        integral *= 2
    elif  integral_type == "dblquad":
        # contour_len[0] to contour[-1] SHOULD be a complete integral on the closed curve; i.e. sxx(contour_len[0]) = sxx(contour_len[-1]) etc.
        val, err = scipy.integrate.dblquad(writhe_integrand, ss[0], ss[-1], lambda x: ss[0], lambda x: ss[-1], args = (sxx, syy, szz, ss[-1]))
        print(sys.stderr, err)
        integral = val
    elif integral_type == "chopped dblquad":
        integral = 0
        for ss_coarse in np.linspace(ss[0], ss[-1], 10):
            for ss_coarse_prime in np.linspace(ss[0], ss[-1], 10):
                val, err = scipy.integrate.dblquad(writhe_integrand, ss_coarse, ss_coarse + float(ss[-1]-ss[0])/9, lambda x: ss_coarse_prime, lambda x: ss_coarse_prime + float(ss[-1]-ss[0])/9, args = (sxx, syy, szz, contour_len[-1]))
                print(err)
                integral += val
    elif integral_type == "quad":
        integral, err = scipy.integrate.quad(writhe_integrand2, ss[0], ss[-1], args = (sxx, syy, szz, ss[0], ss[-1]), limit = 100, epsabs = 1e-5, epsrel = 0)
    elif integral_type == "simps":
        srange = range(len(ss))
        integrand = [[] for ii in srange]
        for ii in srange:
            for jj in srange:
                # skip ii=jj
                if ii == jj:
                    triple_scalar_product = 0
                else:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = np.dot(np.cross(tt[ii], tt[jj]), diff_frac)
                integrand[ii].append(triple_scalar_product)
        integral = scipy.integrate.simps(scipy.integrate.simps(integrand, ss), ss)
    else:
        assert False
        
    writhe = float(integral) / (4*np.pi)

    return writhe



strand = base.Strand()
get_base_spline(strand)
