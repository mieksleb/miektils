#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:41:58 2019

A python script which defines all the tools required when using miektils

@author: michaelselby
"""
import numpy as np
import quaternion as quat
import base
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
import readers
from scipy import interpolate
import numba
from numba import njit
from numba import jit


pi2 = 2*np.pi

def get_array(line):                #turns a line in a trajectpory file into a numpy array
    return np.array([line[0][:-1], line[1][:-1], line[2]])


def get_centre(A, B):               #find the centre point of two quantities, for minimum image convention see get_pos_midpoint
    return 0.5*(A+B)


def vector_string(array):           #turns a numpy array into a readable string
    return str(array[0])+", "+str(array[1])+", "+str(array[2])


def get_pos(file, line):            #gets a numpy array position from line in file 
    traj_line = file[int(line)].split()     
    return get_array(traj_line)


def quat_rot(rate, t, axis, vector): #rotates vector about axis using quaternion rotation
    vec =  np.array([0.] + vector)
    rot_axis = np.array([0.] + axis)
    axis_angle = (rate*t*0.5) * rot_axis/np.linalg.norm(rot_axis)
    vec_quat = quat.quaternion(*vec)
    qlog = quat.quaternion(*axis_angle)
    q = np.exp(qlog)
    new_vec = q * vec_quat * np.conjugate(q)
    return new_vec.imag

def rot(omega, t, A, B):  #rotates vector B about axis A using Euler-Rodrigues formula
    return A*(np.dot(A,B))+np.cos(omega*t)*np.cross((np.cross(A, B)), A)+np.sin(omega*t)*(np.cross(A,B))


#takes a snippet of a trajectory file at every (printed) timestep and returns a list of these snippets 
#hence returns a list of lists of arrays

def snippets(file, N, steps): 
    snippets = []           
    for i in range(steps):
        snippet = []
        for j in range(N):
            snippet.append(get_pos(file, i*(N+3)+j+3))
        snippets.append(snippet)
    return snippets


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

  

    
def get_bb_midpoint(file, bp):
    bbs = []
    bases = []
    comps = []
    for i in range(2*bp):
        bases.append(get_pos(file, i*(2*bp+3)+i+3))
        comps.append(get_pos(file, 2*bp - (i*(2*bp+3)+i+3)))
    bbs = get_centre(bases, comps)
        
    return bbs


def get_spline(base_pos, k = 3, s = 0, per = True, reverse = False):
    """
    return a cartesian spline that represents a fit through the backbone of a duplex

    args:
    base_pos - numpy array of base positions, a strand based version can be found in origami_utils.py :-)
    k - order (default 3, i.e. cubic)
    per - Boolean, True = assume circular
    """
    import scipy
 
    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in base_pos]
    yy = [vec[1] for vec in base_pos]
    zz = [vec[2] for vec in base_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 10, per = per)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 10, per = per)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3 , s = 10, per = per)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)




def spline2frenet(splinex, spliney, splinez, s):
    import scipy.interpolate
    
    def r(s):
        return np.array([scipy.interpolate.splev(s, msxx, 0), scipy.interpolate.splev(s, msyy, 0), scipy.interpolate.splev(s, mszz, 0)])
    
    def T(s):
        return np.array([scipy.interpolate.splev(s, msxx, 1), scipy.interpolate.splev(s, msyy, 1), scipy.interpolate.splev(s, mszz, 1)])

    def N(s):
        return np.array([scipy.interpolate.splev(s, msxx, 2), scipy.interpolate.splev(s, msyy, 2), scipy.interpolate.splev(s, mszz, 2)])
    
    def B(s):
        return np.cross(T(s), N(s))
    
    return np.array[(T(s), N(s), B(s))]
      
'''

This function plots a spline fit and the nucleotide positions.
Function takes two b-spline objects and 2 np.arrays defining the bp positions and the number of bp
returns an axis object (i.e. requires plt.show() after when called)

'''

def plot_spline(spline1, spline2, base_pos1, base_pos2, bp):

    fig2 = plt.figure()
    
    ax = fig2.add_subplot(111, projection='3d')
    
    for i in range(len(base_pos1)):
#        ax.scatter(bbpos1[i][0], bbpos1[i][1], bbpos1[i][2], s=1, marker=",")
        ax.scatter(base_pos1[i][0], base_pos1[i][1], base_pos1[i][2], s=1, marker=",")
    
    for i in range(len(base_pos2)):
#        ax.scatter(bbpos2[i][0], bbpos2[i][1], bbpos2[i][2], s=1, marker=",")
        ax.scatter(base_pos2[i][0], base_pos2[i][1], base_pos2[i][2], s=1, marker=",")
    
    tck, u = spline1[0], spline1[1]
    x_knots, y_knots, z_knots = interpolate.splev(range(0,bp+1), tck[0]), interpolate.splev(range(0,bp+1), tck[1]), interpolate.splev(range(0,bp+1),tck[2])
    #u_fine = np.linspace(0,u,100)
    #x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck[0])[1],interpolate.splev(u_fine, tck[1])[1],interpolate.splev(u_fine, tck[2])[1]  
    ax.plot(x_knots, y_knots, z_knots)
    
    tck2, u2 = spline2[0], spline2[1]
    x_knots2, y_knots2, z_knots2 = interpolate.splev(range(0,bp+1), tck2[0]), interpolate.splev(range(0,bp+1), tck2[1]), interpolate.splev(range(0,bp+1),tck2[2])
    #u_fine2 = np.linspace(0,u,100)
    #x_fine2, y_fine2, z_fine2 = interpolate.splev(u_fine2, tck2[0])[1],interpolate.splev(u_fine2, tck2[1])[1],interpolate.splev(u_fine2, tck2[2])[1]  
    ax.plot(x_knots2, y_knots2, z_knots2)
    
    return ax.plot(x_knots, y_knots, z_knots), ax.plot(x_knots2, y_knots2, z_knots2) 



'''
Plots a projection of a list of splines and a list of positions

xper, yper and zper are booleans which if set to TRUE project out that particular dimension

if best_proj = True, then xper, yper and zper are found using the best_projection_boolean function
''' 
     
def plot_2D_spline(splines, pos, spline_domain, best_proj = True, xper=True, yper=False, zper=False, marker_size = 1):

    fig, ax = plt.subplots()
    
    x = []
    y = []
    z = []
    plots = []
    tcks = []
    us = []
    
    npoints = 1000

    N = 20
    # points for annotation, for a point in points, point[0] is the spline_domain position,
    # and point[1] is the index that this point occurs in the list spline_domain
    points  = [[int(spline_domain[int(i*(npoints-1)/N)]), int(i*(npoints-1)/N)] for i in range(0,N)]
    

    
    if best_proj == True:
        xper = best_projection_boolean(pos[0])[0]
        yper = best_projection_boolean(pos[0])[1]
        zper = best_projection_boolean(pos[0])[2]
    
    
    for j in range(len(splines)):   
        tcks.append(splines[j][0])
        us.append(splines[j][1])
        x.append(interpolate.splev(us[j], tcks[j][0]))
        y.append(interpolate.splev(us[j], tcks[j][1]))
        z.append(interpolate.splev(us[j], tcks[j][2]))
 


    if xper == True:
        for i in range(len(pos)):
            for j in range(len(pos[i])):
                plots.append(ax.scatter(pos[i][j][1], pos[i][j][2], s=marker_size, marker="o"))
        for i in range(len(splines)):      
            x, y, z = interpolate.splev(spline_domain, tcks[i][0]), interpolate.splev(spline_domain, tcks[i][1]), interpolate.splev(spline_domain,tcks[i][2])
            plots.append(ax.plot(y, z))
            for boi in points:
                plots.append(ax.annotate(boi[0], (y[boi[1]], z[boi[1]])))
            
    elif yper == True:
        for i in range(len(pos)):
            for j in range(len(pos[i])):
                plots.append(ax.scatter(pos[i][j][0], pos[i][j][2], s=marker_size, marker="o"))
        for i in range(len(splines)):      
            x, y, z = interpolate.splev(spline_domain, tcks[i][0]), interpolate.splev(spline_domain, tcks[i][1]), interpolate.splev(spline_domain,tcks[i][2])
            plots.append(ax.plot(x, z))
            for boi in points:
                plots.append(ax.annotate(boi[0], (x[boi[1]], z[boi[1]])))
            
            
    elif zper == True:
        for i in range(len(pos)):
            for j in range(len(pos[i])):
                plots.append(ax.scatter(pos[i][j][0], pos[i][j][1], s=marker_size, marker="o"))
        for i in range(len(splines)):      
            x, y, z = interpolate.splev(spline_domain, tcks[i][0]), interpolate.splev(spline_domain, tcks[i][1]), interpolate.splev(spline_domain,tcks[i][2])
            plots.append(ax.plot(x, y))
            for boi in points:
                plots.append(ax.annotate(boi[0], (x[boi[1]], y[boi[1]])))
 

    return plots




'''
Takes an array of n splines, and an array of m bp positions, and returns an ax.plot object
'''
def plot_splines(splines, base_pos, spline_domain, npoints = 1000, annotate = False, marker_size=1):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig2 = plt.figure()
    
    ax = fig2.add_subplot(111, projection='3d')
    
    N = 1
    points  = [[int(spline_domain[int(i*(npoints-1)/N)]), int(i*(npoints-1)/N)] for i in range(0,N)]
     

    plots = []
    for j in range(len(splines)):
        tck, u = splines[j][0], splines[j][1]
        x, y, z = interpolate.splev(spline_domain, tck[0]), interpolate.splev(spline_domain, tck[1]), interpolate.splev(spline_domain,tck[2])
        plots.append(ax.plot(x, y, z))
        
    for i in range(len(base_pos)):
        for j in range(len(base_pos[i])):
            plots.append(ax.scatter(base_pos[i][j][0], base_pos[i][j][1], base_pos[i][j][2], s=marker_size, marker="o"))
    
    
    if annotate:
        for point in points:
            ax.text(x[point[1]], y[point[1]], z[point[1]],  '%s' % (str(point[0])), size=20)
#            plots.append(ax.annotate(point[0], [x[point[1]], y[point[1]], z[point[1]]] ))
                   
    
    return plots 





def plot_surface(f, x_array, y_array, xlabel, ylabel, zlabel):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    def fun(x, y):
        return f(x,y)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x,y = x_array, y_array
    X, Y = np.meshgrid(x, y)
    zs = np.array(fun(np.ravel(X), np.ravel(Y)))
    Z = zs.reshape(X.shape)
    
    ax.plot_surface(X, Y, Z)
    
    ax.set_xlabel(str(xlabel))
    ax.set_ylabel(str(ylabel))
    ax.set_zlabel(str(zlabel))
    
    plt.show()

    
    
    
'''
All these functions calculate the triple product of three vectors a, b and c 

using the quaternion 
triple product who's real part is equal to the negative of the normal triple product

Note: multidet proven to be the fastest

'''

#using quaternions
def triple_product(a, b, c): 
    veca = np.array([0.] + a)
    vecb =  np.array([0.] + b)
    vecc =  np.array([0.] + c)

    vec_quata = quat.quaternion(*veca)
    vec_quatb = quat.quaternion(*vecb)
    vec_quatc = quat.quaternion(*vecc)

    new_vec =  -vec_quata * vec_quatb * vec_quatc
    return new_vec.real


a = np.array([11,0,1])
b = np.array([1,5,1])
c = np.array([2,2,2])


#algebraic
def multidet_normal(u,v,w):
    n=a.shape[0]
    d=np.empty(n)
    d=\
    u[0]*(v[1]*w[2]-v[2]*w[1])+\
    u[1]*(v[2]*w[0]-v[0]*w[2])+\
    u[2]*(v[0]*w[1]-v[1]*w[0])  # 14 operations / det
    return d


#using numba C wrapper on algebraic   FASTEST
@njit(fastmath=True)
def multidet(u,v,w):
    d=np.empty(3)
    d=\
    u[0]*(v[1]*w[2]-v[2]*w[1])+\
    u[1]*(v[2]*w[0]-v[0]*w[2])+\
    u[2]*(v[0]*w[1]-v[1]*w[0])  # 14 operations / det
    return d


#numpy
def trip(a,b,c):
    return np.dot(a,np.cross(b,c))


@njit(fastmath=True)
def norm(vec):
    return vec / np.sqrt(np.dot(vec,vec))


'''
another spline function, possibly redundant
'''
def vecs2spline(vecs, per, k):
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


def get_pos_midpoint(r1, r2, box):
    """
    return the midpoint of two vectors r1 and r2
    
    use the minimum image: i.e. make sure we get a sensible answer if the                                                                     
    positions r1 and r2 correspond to nucleotides which were put at opposite                                                                  
    ends of the box due to pbc's                                                                                                            
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    return r1 - min_distance(r1, r2, box)/2





@njit(fastmath=True)
def discrete_dbl_int(tt, uu, duu, xx, yy, zz, dx, dy, ss, circular = True):
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
                diff = np.array([float(xx[ii]-xx[jj]), float(yy[ii] - yy[jj]), float(zz[ii] - zz[jj])])
                diff_mag = np.sqrt(np.dot(diff, diff))
                diff_frac = diff / (diff_mag ** 3)
                triple_scalar_product = multidet((tt[ii]), tt[jj], diff_frac)
                writhe_integral += triple_scalar_product
        # multiply by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
    twist_integral *= dx / (2 * np.pi)
    writhe_integral *= 1 * dx * dy / (2 * np.pi)
    return twist_integral, writhe_integral




"""
Inputting a list of arrays, this function will return 3 booleans: xper, yper, zper 
One value being TRUE if it is the best dimension to project out of, i.e. if output reads
False, True, False, then you should project out the y component
"""
def best_projection_boolean(pos1):
    
    xmax = max(pos1, key=lambda item: item[0])[0]
    xmin = min(pos1, key=lambda item: item[0])[0]
    ymax = max(pos1, key=lambda item: item[1])[1]
    ymin = min(pos1, key=lambda item: item[1])[1]
    zmax = max(pos1, key=lambda item: item[2])[2]
    zmin = min(pos1, key=lambda item: item[2])[2]
    
    xrange = xmax - xmin 
    yrange = ymax - ymin
    zrange = zmax - zmin
    
    best_projection = np.argmin([xrange, yrange, zrange])
        
        
    if best_projection == 0:
        xper = True
        yper = False
        zper = False
    elif best_projection == 1:
        xper = False
        yper = True
        zper = False
    elif best_projection == 2:
        xper = False
        yper = False
        zper = True
    return xper, yper, zper




'''
Finds the member of the list 'points' which is closest to 'point'
'''
def closest_point(point, points):
    points = np.asarray(points)
    deltas = points - point
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)


#tt = np.array(inp_list, dtype=np.float32)
tt = np.array([[1.0,1.0,1.0], [1.0,0.0,0.0], [2.0,0.0,1.0]])
xx = [1,2,2]
yy = [2,2,1]
zz = [2,3,5]
ds = 1
xrange = yrange = range(0,3)
ss = [1,2,3,4]

#yes = len(xx)
#hootenanny = discrete_dbl_int(tt, xx, yy, zz, ds, ds, ss)

