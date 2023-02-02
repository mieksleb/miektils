#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:04:17 2021

Tools_fast_math

All functions in tools.py which can be sped up using njit are here. 
Functions are generally wrappable in this way if they use numpy arrays.

@author: michaelselby
"""
import numpy as np
import numba
from numba import njit
from numba import jit
# from numba.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

# warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
# warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


@njit(fastmath=True)
def multidet(u, v, w):
    #d = np.empty(3)
    d =\
        u[0]*(v[1]*w[2]-v[2]*w[1]) +\
        u[1]*(v[2]*w[0]-v[0]*w[2]) +\
        u[2]*(v[0]*w[1]-v[1]*w[0])  # 14 operations / det
    return d

@njit(fastmath=True)
def multidet2(u, v, w):
    return np.linalg.det(np.dstack([u,v,w]))


@njit(fastmath=True)
def rot_fast(v, k, theta):  # rotates vector v about axis k using Euler-Rodrigues formula
    return v*np.cos(theta)+np.cross(k, v)*np.sin(theta)+k*np.dot(k, v)*(1-np.cos(theta))


@njit(fastmath=True)
def trip(a, b, c):
    return np.dot(a, np.cross(b, c))


# @njit(fastmath=True)
def normalize(vec):
    return vec / np.sqrt(np.dot(vec, vec))


@njit(fastmath=True)
def norm(vec):
    return np.sqrt(np.dot(vec, vec))
    # return np.linalg.norm(vec)


@njit(fastmath=True)
def get_tangent(r, circular):
    n = np.shape(r)[0]  # number of particles
    t = np.array([r[(i+1) % n, :]-r[i, :] for i in range(n)])
    t /= np.sqrt((t ** 2).sum(-1))[..., np.newaxis]
    if not circular:
        np.delete(t[-1, :])
    return t


# now for the integration which we use numba fastmath capabilities for speed
@njit(fastmath=True)
def discrete_dbl_int(tt, uu, duu, xx, yy, zz, dx, dy, ss, circular=True):
    if circular:
        srange = range(len(ss)-1)
    else:
        srange = range(len(ss))
    twist_integral = 0
    writhe_integral = 0
    for ii in srange:
        triple_scalar_product = multidet(tt[ii], uu[ii], duu[ii])
        #triple_scalar_product = np.dot(tt[ii], np.cross(uu[ii], duu[ii]))
        twist_integral += triple_scalar_product
        for jj in srange:
            # skip ii=jj and use symmetry in {ii, jj}
            if ii > jj:
                diff = np.array(
                    [xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                diff_mag = np.sqrt(np.dot(diff, diff))
                diff_frac = diff / (diff_mag ** 3)
                triple_scalar_product = multidet((tt[ii]), tt[jj], diff_frac)
                writhe_integral += triple_scalar_product
        # multiply by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
    twist_integral *= dx / (2 * np.pi)
    writhe_integral *= 1 * dx * dy / (2 * np.pi)
    return twist_integral, writhe_integral


def get_twist_writhe(spline1, spline2, npoints=1000, circular=False, integral_type="simple"):
    """
    return the twist and writhe for a given configuration and 
    Using integral_type = 'simple' 

    args:
    spline1: list of 3 splines corresponding to x, y and z spline through strand 1's backbone
    spline2: list of 3 splines corresponding to x, y and z spline through strand 2's backbone -- NB the splines should run in the same direction, i.e. one must reverse one of the splines if they c
    from get_base_spline (e.g. use get_base_spline(reverse = True))

    npoints: number of points for the discrete integration
    """
    import scipy
    from scipy.interpolate import splev, splrep

    s1xx, s1yy, s1zz = spline1[0]
    s2xx, s2yy, s2zz = spline2[0]

    smin = spline1[1][0]
    smax = spline1[1][-1]

    # bpi is the base pair index parameter that is common to both splines
    bpi = np.linspace(smin, smax, npoints)

    # find the midpoint between the input splines, as a function of base pair index
    m1xx = (splev(bpi, s1xx) + splev(bpi, s2xx)) / 2
    m1yy = (splev(bpi, s1yy) + splev(bpi, s2yy)) / 2
    m1zz = (splev(bpi, s1zz) + splev(bpi, s2zz)) / 2

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = np.array([np.sqrt((m1xx[ii+1]-m1xx[ii])**2+(m1yy[ii+1]-m1yy[ii])
                       ** 2+(m1zz[ii+1]-m1zz[ii])**2) for ii in range(len(bpi)-1)])

    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    # get the splines as a function of contour length
    msxx = splrep(contour_len, m1xx, k=3, s=0, per=circular)
    msyy = splrep(contour_len, m1yy, k=3, s=0, per=circular)
    mszz = splrep(contour_len, m1zz, k=3, s=0, per=circular)

    xx, yy, zz = splev(ss, msxx), splev(ss, msyy), splev(ss, mszz)

    # find the tangent of the midpoint spline.
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of normalised tangent vectors; __call__(xxx, 1) returns the first derivative
    # the tangent vector is a unit vector
    dmxx, dmyy, dmzz = splev(ss, msxx, 1), splev(
        ss, msyy, 1), splev(ss, mszz, 1)

    tt = np.array((dmxx, dmyy, dmzz)).T

    # we also need the 'normal' vector u(s) which points between the base pairs. (or between the spline fits through the backbones in this case)
    # n.b. these uxx, uyy, uzz are not normalised
    uxx_bpi = splev(bpi, s2xx) - splev(bpi, s1xx)
    uyy_bpi = splev(bpi, s2yy) - splev(bpi, s1yy)
    uzz_bpi = splev(bpi, s2zz) - splev(bpi, s1zz)

    # get the normal vector spline as a function of contour length
    suxx = splrep(contour_len, uxx_bpi, k=3, s=0, per=circular)
    suyy = splrep(contour_len, uyy_bpi, k=3, s=0, per=circular)
    suzz = splrep(contour_len, uzz_bpi, k=3, s=0, per=circular)

    # evaluate the normal vector spline as a function of contour length
    uxx, uyy, uzz = splev(ss, suxx), splev(ss, suyy), splev(ss, suzz)

    uu = np.array((uxx, uyy, uzz)).T

    uu -= np.array([np.dot(tt, uu) * tt for tt, uu in zip(tt, uu)])
    uu /= np.sqrt((uu ** 2).sum(-1))[..., np.newaxis]

    # and finally we need the derivatives of that vector u(s). It takes a bit of work to get a spline of the normalised version of u from the unnormalised one
    nuxx = np.array(np.array([vec[0] for vec in uu]))
    nuyy = np.array(np.array([vec[1] for vec in uu]))
    nuzz = np.array(np.array([vec[2] for vec in uu]))
    nusxx = splrep(ss, nuxx, k=3, s=0, per=circular)
    nusyy = splrep(ss, nuyy, k=3, s=0, per=circular)
    nuszz = splrep(ss, nuzz, k=3, s=0, per=circular)

    duxx, duyy, duzz = splev(ss, nusxx, 1), splev(
        ss, nusyy, 1), splev(ss, nuszz, 1)

    duu = np.array((duxx, duyy, duzz)).T

    ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
    # do the integration w.r.t. s

    twist, writhe = discrete_dbl_int(
        tt, uu, duu, xx, yy, zz, ds, ds, ss, circular=True)

    return twist, writhe


def get_twist_writhe_TEP(centre_line, tan_vals, normals, npoints=1000, circular=False, integral_type="simple"):
    """
    return the twist and writhe for a given configuration and 
    Using integral_type = 'simple' 

    args:


    npoints: number of points for the discrete integration
    """

    import scipy.interpolate

    m1xx = [vec[0] for vec in centre_line]
    m1yy = [vec[1] for vec in centre_line]
    m1zz = [vec[2] for vec in centre_line]

 # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((m1xx[ii+1]-m1xx[ii])**2+(m1yy[ii+1]-m1yy[ii])
                       ** 2+(m1zz[ii+1]-m1zz[ii])**2) for ii in range(len(m1xx)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((m1xx[ii+1]-m1xx[ii])**2+(m1yy[ii+1]-m1yy[ii])
                       ** 2+(m1zz[ii+1]-m1zz[ii])**2) for ii in range(len(m1xx)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # get the splines as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, m1xx, k=3, s=0, per=circular)
    msyy = scipy.interpolate.splrep(contour_len, m1yy, k=3, s=0, per=circular)
    mszz = scipy.interpolate.splrep(contour_len, m1zz, k=3, s=0, per=circular)

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
    uxx_bpi = [vec[0] for vec in normals]
    uyy_bpi = [vec[1] for vec in normals]
    uzz_bpi = [vec[2] for vec in normals]

    # get the normal vector spline as a function of contour length
    suxx = scipy.interpolate.splrep(
        contour_len, uxx_bpi, k=3, s=0, per=circular)
    suyy = scipy.interpolate.splrep(
        contour_len, uyy_bpi, k=3, s=0, per=circular)
    suzz = scipy.interpolate.splrep(
        contour_len, uzz_bpi, k=3, s=0, per=circular)

    # evaluate the normal vector spline as a function of contour length
    uxx = scipy.interpolate.splev(ss, suxx)
    uyy = scipy.interpolate.splev(ss, suyy)
    uzz = scipy.interpolate.splev(ss, suzz)

    uu = list(range(len(ss)))
    for ii in list(range(len(ss))):
        uu[ii] = np.array([uxx[ii], uyy[ii], uzz[ii]])
        uu[ii] = uu[ii] - np.dot(tt[ii], uu[ii]) * tt[ii]
        # the normal vector should be normalised
        uu[ii] = normalize(uu[ii])

    # and finally we need the derivatives of that vector u(s). It takes a bit of work to get a spline of the normalised version of u from the unnormalised one
    nuxx = [vec[0] for vec in uu]
    nuyy = [vec[1] for vec in uu]
    nuzz = [vec[2] for vec in uu]

    nusxx = scipy.interpolate.splrep(ss, nuxx, k=3, s=0, per=circular)
    nusyy = scipy.interpolate.splrep(ss, nuyy, k=3, s=0, per=circular)
    nuszz = scipy.interpolate.splrep(ss, nuzz, k=3, s=0, per=circular)
    duxx = scipy.interpolate.splev(ss, nusxx, 1)
    duyy = scipy.interpolate.splev(ss, nusyy, 1)
    duzz = scipy.interpolate.splev(ss, nuszz, 1)
    duu = list(range(len(ss)))
    for ii in list(range(len(ss))):
        duu[ii] = np.array([duxx[ii], duyy[ii], duzz[ii]])

    ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
    # do the integration w.r.t. s

    twist, writhe = discrete_dbl_int(
        tt, uu, duu, xx, yy, zz, ds, ds, ss, circular=True)

    return twist, writhe


@njit(fastmath=True)
def radius_of_gyration(r):
    n = np.shape(r)[0]
    rad2 = 0
    for i in range(n):
        for j in range(n):
            rad2 += np.dot(r[i, :]-r[j, :], r[i, :]-r[j, :])
    rad2 /= (2*n**2)
    rad = np.sqrt(rad2)
    return rad



# @njit(fastmath=True)
def difference(r, circular):
    bp = np.shape(r)[0]
    diff = [r[(j+1), :] - r[j, :] for j in range(bp-1)]
    diff.append(r[0, :] - r[-1, :])
    diff = np.array(diff)
    return diff


# @njit(fastmath=True)
def disc_curvature(r, circular):
    """
    Parameters
    ----------
    r : numpy array (bp,3)
        array of positions, r(i,j) where i is base-pair index and j is coordinate index
    circular : boolean
        true if DNA is closed
        false of DNA is linear

    Returns
    -------
    curv : numpy array (bp)
        curvature at each base-pair index

    """
    bp = np.shape(r)[0]
    curv = np.zeros(bp)
    diff = difference(r, circular)
    delta_s = np.sqrt((diff ** 2).sum(-1))[..., np.newaxis][:,0]
    # diff = np.diff(r, axis=0, append=0)
    # diff /= np.sqrt((diff ** 2).sum(-1))[..., np.newaxis]
    diff = np.array([ diff / norm(diff) for diff in diff])
    
    
    length = bp - 1
    if circular:
        length += 1
        
    diff2 = difference(diff,circular)
    curv = np.array([norm(diff2[i,:])/delta_s[i] for i in range(length)])
    # curv = np.sqrt( (dtanz*tany - dtany*tanz) **2 + (dtanx*tanz - dtanz*tanx) **2 + (dtany*tanx - dtanx*tany) **2  ) / (tanx **2 + tany **2 + tanz **2)**(1.5)
    
    return curv



@njit(fastmath=True)
def get_angles(r, extent, circular):
    """
    Parameters
    ----------
    r : numpy array (bp,3)
        array of positions, r(i,j) where i is base-pair index and j is coordinate index
    circular : boolean
        true if DNA is closed
        false of DNA is linear

    Returns
    -------
    angles : numpy array (bp)
        angle formed at each base-pair index given extent in each direction

    """
    bp = np.shape(r)[0]
    curv = np.zeros(bp)
    diff1 = np.array([r[(j+extent) % bp, :] - r[j, :] for j in range(bp)])
    diff1 /= np.sqrt((diff1 ** 2).sum(-1))[..., np.newaxis]

    diff2 = np.array([- r[(j-extent) % bp, :] + r[j, :] for j in range(bp)])
    diff2 /= np.sqrt((diff2 ** 2).sum(-1))[..., np.newaxis]
    length = bp-1
    if circular:
        length += 1
    angles = np.array([np.arccos(np.dot(diff1[i], diff2[i]))
                      for i in range(length)])
    angles *= 180/np.pi
    return angles






