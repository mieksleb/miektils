#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 17:07:55 2021

Miek_generator

A self-contained configuration and topology generator for people that hate oxDNA utils

The idea is to create some centre-line structure which you can then form a duplex around


@author: michaelselby
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from random import choice
from scipy.interpolate import splev, splrep



speedup = True

if speedup:
    from tools_fast_math import get_twist_writhe,multidet,norm
else:
    from tools import multidet,norm,get_twist_writhe




'''
List of all the different structures we can make

'''

def circle(t):
        x = np.cos(t)
        y = np.sin(t)
        z = 0
        return x,y,z

def figure_eight(t):
        x = (2+np.cos(2*t))*np.cos(3*t)
        y = (2+np.cos(2*t))*np.sin(3*t)
        z = np.sin(4*t)
        return x,y,z

def right_trefoil(t):
        x = (np.sin(t)+2*np.sin(2*t))
        y = (np.cos(t)-2*np.cos(2*t))
        z = -np.sin(3*t)
        return x,y,z
    
    
def bernoulli_lemniscate(t):
        x = np.cos(t)/(1+(np.sin(t))**2)
        y = np.sin(t)*np.cos(t)/(1+(np.sin(t))**2)
        z = 0.1*np.sin(t)
        return x,y,z
    
    


'''
Parameters and checks

'''
func = lambda t: right_trefoil(t)
bp = 300
circular = True         # circular boolean is true when structure is closed
periodic = True         # periodic boolean is true when x(n)=x(1), i.e. the final position is equal to the first
sequence_input = False  # sequence = True means a sequence file will be given
strands = 2
oxdna_unit_len = 0.8518 # nm
lpbp = 0.33 / oxdna_unit_len        # length per base pair in simulation units
dup_rad = (1.15 / oxdna_unit_len)/2 # duplex radius in simualtion units
L = bp*lpbp                         # contour length 
base_length = 0.7                   # distance from backbone to base
centre_line = []
pitch = 10.36                       # helical pitch of dsDNA
npoints = bp

Lk0 = bp/pitch # this is the relaxed linking number, deltaLk=0
deltaLk = 0
    
theta0 = Lk0*2*np.pi/bp

if circular == False:
    if periodic == False:
        print('Configration cannot be periodic but not circular!')

sequence = []
if sequence_input == False:
    print('No sequence given, generating a random sequence')
    for count in range(bp):
       sequence.append(choice("CGTA"))
       
def generate(bp,centre_line,tan_vals,theta,v1):
    """

    Parameters
    ----------
    bp : integer
        number of base pairs
    centre_line : list of numpy arrays
        ith member of list is numpy array of centreline position
    tan_vals : list
        list of tangent values
    theta : real
        angle per base pair of rotation
    v1 : numpy array (3d vector)
        DESCRIPTION.

    Returns
    -------
    strand1pos : list
        list of positions of strand1
    strand2pos : list
        list of positions of strand1
    normals1 : TYPE
        DESCRIPTION.
    normals2 : TYPE
        DESCRIPTION.

    """
    # skip first base pairs
    for i in range(1,bp):
        
        v2 = v1
        # remove component of v2 in direction of tan
        v2 -= v2*(np.dot(v2,tan_vals[i]))
        v2 /= np.linalg.norm(v2)
        
        v2 = rot(v2,tan_vals[i],theta)
        v2 /= np.linalg.norm(v2)    
        
        '''
        # rotate the previous v1 by theta, then project (other option)
        v2 = rot(v1,tan_vals[i],theta)
        v2 /= np.linalg.norm(v2)
        
        # remove component of v2 in direction of v1
        v2 -= v2*(np.dot(v2,v1))
        v2 /= np.linalg.norm(v2)
        '''
             
        rb1 = centre_line[i] + dup_rad*v2
        rb2 = centre_line[i] - dup_rad*v2
        
        normals1.append(-v2)
        normals2.append(v2)
        strand1pos.append(rb1)
        strand2pos.append(rb2)
        v1 = v2
    return strand1pos,strand2pos,normals1,normals2,v2

 
    
def get_spline(bb_pos,reverse):
    if reverse:
        bb_pos.reverse()
    if circular:
        if reverse:
            bb_pos.append(bb_pos[-1])
        else:
            bb_pos.append(bb_pos[0])
    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in bb_pos]
    yy = [vec[1] for vec in bb_pos]
    zz = [vec[2] for vec in bb_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = splrep(range(len(xx)), xx, k = 3, s = 0, per = circular)
    spline_yy = splrep(range(len(yy)), yy, k = 3, s = 0, per = circular)
    spline_zz = splrep(range(len(zz)), zz, k = 3, s = 0, per = circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

comp_sequence1 = []
for base in sequence:
    if base=="A":
        comp_sequence1.append("T")
    elif base=="T":
         comp_sequence1.append("A")
    elif base=="C":
         comp_sequence1.append("G")
    else:
         comp_sequence1.append("C")
         
comp_sequence = []
for i in range(bp):
    comp_sequence.append(comp_sequence1[bp-1-i])
    
def rot(vec,axis,theta):
    # rotates a vector vec about axis by counter-clockwise angle theta
    return axis*(np.dot(vec,axis))*(1-np.cos(theta))+(np.cross(axis,vec))*np.sin(theta)+vec*np.cos(theta)

if periodic == True:
    length = bp+1
else:
    length = bp


# create the arclength range
svals = np.linspace(0,L,length)
ds = L/bp



'''
Constructing you centre-line

We begin with a general configuration describing the centre line of the dsDNA

We then rescale such that the contour length matches the expected length for a sequence of the length
'''

n=1
io = 1
while io!=3:
    centre_line = []
    for s in svals:
        t = s*2*np.pi/L
        x,y,z = func(t)
        x=n*x
        y=n*y
        z=n*z
        point = np.array([x,y,z])
        centre_line.append(point) 
    # We must now calculate the contour length of this shape
    xx = [vec[0] for vec in centre_line]
    yy = [vec[1] for vec in centre_line]
    zz = [vec[2] for vec in centre_line]
    xrange = max(xx)-min(xx)
    yrange = max(yy)-min(yy)
    zrange = max(zz)-min(zz)
    delta_s = [np.sqrt((xx[ii+1]-xx[ii])**2+(yy[ii+1]-yy[ii])**2+(zz[ii+1]-zz[ii])**2) for ii in range(bp-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)
    total_len = contour_len[-1]
    n = L/total_len
    io+=1



box = max(xrange,yrange,zrange)


# if the structure is circular and periodic, then we remove the last element of the centre_lines
if periodic == True:
    del centre_line[-1]
    


# calculate the unit tangents at every base pair
tan_vals = []
old_point = centre_line[-1]
for point in centre_line:
    diff = np.array([point[0]-old_point[0],point[1]-old_point[1],point[2]-old_point[2]])
    diff /= np.linalg.norm(diff)
    old_point = point
    tan_vals.append(diff)
    
# need a vector perpendicular to the initial tangent tan_vals[0]
v1 = np.random.random_sample(3)
v1 -= tan_vals[0] * (np.dot(tan_vals[0], v1))
v1 /= np.sqrt(sum(v1*v1))
v10 = v1    # save the very first vector

phi = 0


'''
This part of the code ensures that the configuration has the correct linking number
'''

for j in range(0,2):
    strand1pos = []
    strand2pos = []
    normals1 = []
    normals2 = []
    
    rb1 = centre_line[0] + dup_rad*v1
    rb2 = centre_line[0] - dup_rad*v1

    normals1.append(-v1)
    normals2.append(v1)
    strand1pos.append(rb1)
    strand2pos.append(rb2)

    theta = theta0 + phi
    
    # skip first base pairs
    strand1pos,strand2pos,normals1,normals2,v2 = generate(bp, centre_line, tan_vals, theta, v1)
    
    spline1 = get_spline(strand1pos,reverse=False)
    spline2 = get_spline(strand2pos,reverse=False)
    twist, writhe = get_twist_writhe(spline1, spline2, npoints, circular = True, integral_type = "simple")
    Lk = twist + writhe
    print('twist is '+ str(twist))
    print('writhe is '+ str(writhe))
    print('linking number is  '+ str(Lk))
    
    deltaLk_actual = Lk - Lk0
    difference = deltaLk_actual-deltaLk
    phi = -difference*2*np.pi/bp
    
    
# The structure now has the correct linking number! We now need to check the
# angle between the last base-pair and the first
    
v2 -= v2*(np.dot(v2,tan_vals[-1]))
v2 /= np.linalg.norm(v2)
last_ang = np.arccos(np.dot(v2,v10))
print(last_ang)
print(theta)
delta_ang = last_ang - theta
delta_ang_pbp = delta_ang/bp # chnage in angle per base pair
theta += delta_ang_pbp
    


strand1pos = []
strand2pos = []
normals1 = []
normals2 = []

rb1 = centre_line[0] + dup_rad*v1
rb2 = centre_line[0] - dup_rad*v1


normals1.append(-v1)
normals2.append(v1)
strand1pos.append(rb1)
strand2pos.append(rb2)

strand1pos,strand2pos,normals1,normals2,v2 = generate(bp, centre_line, tan_vals, theta, v10)



spline1 = get_spline(strand1pos,reverse=False)
spline2 = get_spline(strand2pos,reverse=False)
twist, writhe = get_twist_writhe(spline1, spline2, npoints = 1000, circular = True, integral_type = "simple")
Lk = twist + writhe
print('twist is '+ str(twist))
print('writhe is '+ str(writhe))
print('linking number is  '+ str(Lk))

deltaLk = Lk - Lk0

    
def plot(list_vecs,scatter=False):
    vecs = np.array(list_vecs)
    x = vecs[:,0]
    y = vecs[:,1]
    z = vecs[:,2]
    if scatter==False:
        return ax.plot(x,y,z)
    else:
        return ax.scatter(x,y,z)

    

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-L/4, L/4)
ax.set_ylim3d(-L/4, L/4)
ax.set_zlim3d(-L/4, L/4)

plot(strand1pos,scatter=False)
plot(strand2pos,scatter=False)
plot(centre_line)

plt.show()







'''
Now we write the topology and configuration files
'''


top_file_name = "top.top"
conf_file_name = "conf.dat"

top = open(top_file_name,"w")
conf = open(conf_file_name,"w")

top.write(str(2*bp)+ " " + str(strands) + "\n")
for s in range(strands):
    s = s + 1
    if s==1:
        if circular == True:
            top.write(str(s) + " " + sequence[0] + " " + str(bp-1) + " " + str(1) + "\n")
        else:
            top.write(str(s) + " " + sequence[0] + " " + str(-1) + " " + str(1) + "\n")
        for i in range(1,bp-1):
            top.write(str(s) + " " + sequence[i] + " " + str(i-1) + " " + str(i+1) + "\n")
        if circular == True:
            top.write(str(s) + " " + sequence[-1] + " " + str(bp-2) + " " + str(0) + "\n")
        else:
            top.write(str(s) + " " + sequence[-1] + " " + str(-1) + " " + str(1+bp) + "\n")
    else:
        if circular == True:
            top.write(str(s) + " " + comp_sequence[0] + " " + str(2*bp-1) + " " + str(bp+1) + "\n")
        else:
            top.write(str(s) + " " + comp_sequence[0] + " " + str(-1) + " " + str(1+bp) + "\n")
        for i in range(1,bp-1):
            top.write(str(s) + " " + comp_sequence[i] + " " + str(i-1+bp) + " " + str(i+1+bp) + "\n")
        if circular == True:
            top.write(str(s) + " " + comp_sequence[-1] + " " + str(2*bp-2) + " " + str(bp) + "\n")
        else:
            top.write(str(s) + " " + comp_sequence[-1] + " " + str(-1) + " " + str(1+bp) + "\n")

    
top.close()



conf.write("t = 0" + "\n")
conf.write("b = "+str(box)+" " +str(box)+" "+str(box) + "\n")
conf.write("E = 0" + "\n")
for s in range(strands):
    s = s + 1
    if s==1:
        for i in range(0,bp):
            conf.write(str(strand1pos[i][0]) + " " + str(strand1pos[i][1]) + " " + str(strand1pos[i][2])  \
                        +" " + str(normals1[i][0]) + " " + str(normals1[i][1]) + " " + str(normals1[i][2]) \
                             + " " +str(-tan_vals[i][0]) + " " + str(-tan_vals[i][1]) + " " + str(-tan_vals[i][2]) \
                                 + " 0 0 0 0 0 0" +"\n")
    else:
        for i in range(0,bp):
            conf.write(str(strand2pos[bp-i-1][0]) + " " + str(strand2pos[bp-i-1][1]) + " " + str(strand2pos[bp-i-1][2]) \
                       +" " + str(normals2[bp-i-1][0]) + " " + str(normals2[bp-i-1][1]) + " " + str(normals2[bp-i-1][2]) \
                           + " " + str(-tan_vals[bp-i-1][0]) + " " + str(-tan_vals[bp-i-1][1]) + " " + str(-tan_vals[bp-i-1][2]) \
                               + " 0 0 0 0 0 0" +"\n")

conf.close()



















