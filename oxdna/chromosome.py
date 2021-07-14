#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 10:42:47 2021

@author: michaelselby
"""
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
from miek_generator import write_conf,write_top,plot,generate,gen_seq,pos_2_spline,evaluate_spline
from hilbertcurve.hilbertcurve import HilbertCurve

speedup = True

if speedup:
    from tools_fast_math import get_twist_writhe,multidet,normalize,rot_fast
else:
    from tools import multidet,normalize,get_twist_writhe,rot


circular = True
periodic = True
bp = 10000
npoints = int(bp/300) # each segment is the length of 300 base pairs
lpbp = 0.3897628551303122
L = bp*lpbp
box = 5
pitch = 10.36



def get_possible_directions(point):
    """Point is in form (x, y, z)"""
    directions = [
        (point[0]+1, point[1], point[2]),  # right
        (point[0]-1, point[1], point[2]),  # left
        (point[0], point[1]+1, point[2]),  # forward
        (point[0], point[1]-1, point[2]),  # backward
        (point[0], point[1], point[2]+1),  # up
        (point[0], point[1], point[2]-1)   # down
    ]
    return directions


def random_walk_3D(N,box):
    Nsteps = range(N)
    current_position = (0, 0, 0)
    visited_points = []
    forbidden_points = []
    for i in range(-box,box):
        for j in range(-box,box):
            forbidden_points.append((i,j,box))
            forbidden_points.append((i,box,j))
            forbidden_points.append((box,i,j))
            forbidden_points.append((i,j,-box))
            forbidden_points.append((i,-box,j))
            forbidden_points.append((-box,i,j))
    for n in Nsteps:
        visited_points.append(current_position)
        all_directions = get_possible_directions(current_position)
        not_visited_directions = [direction for direction in all_directions if direction not in visited_points and direction not in forbidden_points]
        current_position = random.choice(not_visited_directions)

    xp, yp, zp = zip(*visited_points)
    return xp, yp, zp  # returns tuples. If you want lists, just do list(xp), ...


def space_filling3D(N,box):
    Nsteps = range(N)
    current_position = (0, 0, 0)
    visited_points = []
    forbidden_points = []
    for i in range(-box,box):
        for j in range(-box,box):
            forbidden_points.append((i,j,box))
            forbidden_points.append((i,box,j))
            forbidden_points.append((box,i,j))
            forbidden_points.append((i,j,-box))
            forbidden_points.append((i,-box,j))
            forbidden_points.append((-box,i,j))
    for n in Nsteps:
        visited_points.append(current_position)
        all_directions = get_possible_directions(current_position)
        not_visited_directions = [direction for direction in all_directions if direction not in visited_points and direction not in forbidden_points]
        current_position = random.choice(not_visited_directions)

    xp, yp, zp = zip(*visited_points)
    return xp, yp, zp  # returns tuples. If you want lists, just do list(xp), ...


def closed_space_fill(width,height):
    points = []
    start = np.array([0,0,0])
    points.append(start)
    point = start
    for i in range(width):
        for j in range(height):
            points.append(point + np.array([-1,0,0]))
            points.append(point + np.array([-2,0,0]))
            points.append(point + np.array([-2,1,0]))
            points.append(point + np.array([-1,1,0]))
            points.append(point + np.array([0,1,0]))
            points.append(point + np.array([0,2,0]))
            point = point + np.array([0,2,0])
        points.pop()
        point = point + np.array([0,-1,0])
        point = point + np.array([0,0,1])
        points.append(point)
        for j in range(height):
            points.append(point + np.array([-1,0,0]))
            points.append(point + np.array([-2,0,0]))
            points.append(point + np.array([-2,-1,0]))
            points.append(point + np.array([-1,-1,0]))
            points.append(point + np.array([0,-1,0]))
            points.append(point + np.array([0,-2,0]))
            point = point + np.array([0,-2,0])
        points.pop()
        point = point + np.array([0,1,0])
        point = point + np.array([0,0,1])
        points.append(point)
    points.pop()
    point = point + np.array([0,0,-1])
    point = point + np.array([1,0,0])
    points.append(point)
    for i in range(width):
        for j in range(height):
            points.append(point + np.array([1,0,0]))
            points.append(point + np.array([2,0,0]))
            points.append(point + np.array([2,1,0]))
            points.append(point + np.array([1,1,0]))
            points.append(point + np.array([0,1,0]))
            points.append(point + np.array([0,2,0]))
            point = point + np.array([0,2,0])
        points.pop()
        point = point + np.array([0,-1,0])
        point = point + np.array([0,0,-1])
        points.append(point)
        for j in range(height):
            points.append(point + np.array([1,0,0]))
            points.append(point + np.array([2,0,0]))
            points.append(point + np.array([2,-1,0]))
            points.append(point + np.array([1,-1,0]))
            points.append(point + np.array([0,-1,0]))
            points.append(point + np.array([0,-2,0]))
            point = point + np.array([0,-2,0])
        points.pop()
        point = point + np.array([0,1,0])
        point = point + np.array([0,0,-1])
        points.append(point)
    points.pop()
    points.append(start)
    
    return points
        
        
        
    

points = closed_space_fill(2,1)
# if the structure is circular and periodic, then we remove the last element of the centre_lines
if periodic == True:
    del points[-1]

bp = 4900
circular = True         # circular boolean is true when structure is closed
periodic = True         # periodic boolean is true when x(n)=x(1), i.e. the final position is equal to the first
sequence_input = False  # sequence = True means a sequence file will be given
strands = 2
oxdna_unit_len = 0.8518 # nm
lpbp = 0.3897628551303122        # length per base pair in simulation units
dup_rad = 0.6 # duplex radius in simualtion units
L = bp*lpbp                         # contour length 
base_length = 0.7                   # distance from backbone to base
pitch = 10.34                       # helical pitch of dsDNA
npoints = bp


Lk0 = bp/pitch # this is the relaxed linking number, deltaLk=0
deltaLk = 0
    
theta0 = Lk0*2*np.pi/bp

if circular == False:
    if periodic == False:
        print('Configration cannot be periodic but not circular!')
        
sequence = gen_seq(bp)


if periodic == True:
    length = bp+1
else:
    length = bp


# create the arclength range
svals = np.linspace(0,49,length)
ds = L/bp



'''
Constructing you centre-line

We begin with a general configuration describing the centre line of the dsDNA

We then rescale such that the contour length matches the expected length for a sequence of the length
'''



# We must now calculate the contour length of this shape
xx = [vec[0] for vec in points]
yy = [vec[1] for vec in points]
zz = [vec[2] for vec in points]
xrange = max(xx)-min(xx)
yrange = max(yy)-min(yy)
zrange = max(zz)-min(zz)
delta_s = [np.sqrt((xx[ii+1]-xx[ii])**2+(yy[ii+1]-yy[ii])**2+(zz[ii+1]-zz[ii])**2) for ii in range(len(xx)-1)]
contour_len = np.cumsum(delta_s)
contour_len = np.insert(contour_len, 0, 0)
total_len = contour_len[-1]
n = L/total_len

centre_line = []
for point in points:
    point = n*point
    centre_line.append(point)
    

box = 5*max(xrange,yrange,zrange)


# if the structure is circular and periodic, then we remove the last element of the centre_lines
# if periodic == True:
#     del centre_line[-1]


spline = pos_2_spline(centre_line,circular,reverse=False)

    
centre_line2 = evaluate_spline(spline,svals)


strand1pos,strand2pos,normals1,normals2,tan_vals = generate(bp,centre_line2,theta0,circular)

 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


plot(strand1pos,ax,scatter=False)


    

# phi = 0

# '''
# This part of the code ensures that the configuration has the correct linking number
# '''
# print("Correcting the Linking number")
# for j in range(0,1):

#     theta = theta0 + phi
    
#     strand1pos,strand2pos,normals1,normals2,tan_vals = generate(bp,centre_line,theta,circular)
#     spline1 = pos_2_spline(strand1pos,circular,reverse=False)
#     spline2 = pos_2_spline(strand2pos,circular,reverse=False)
#     twist, writhe = get_twist_writhe(spline1, spline2, npoints, circular, integral_type = "simple")
#     Lk = twist + writhe
#     print('twist is '+ str(twist))
#     print('writhe is '+ str(writhe))
#     print('linking number is  '+ str(Lk))
    
#     deltaLk_actual = Lk - Lk0
#     difference = deltaLk_actual-deltaLk
#     phi = difference*2*np.pi/bp
    




    
# # The structure now has the correct linking number and major and minor grooving! We now need to check the
# # angle between the last base-pair and the first
# v10 = normals1[0]
# v2 = normals2[-1]
# v2 -= v2*(np.dot(v2,tan_vals[-1]))
# v2 /= np.linalg.norm(v2)
# last_ang = np.arccos(np.dot(v2,v10))

# delta_ang = last_ang - theta
# delta_ang_pbp = delta_ang/bp # change in angle per base pair
# theta += delta_ang_pbp
# # strand1pos,strand2pos,normals1,normals2,tan_vals = generate(bp,centre_line,theta,circular,periodic)

# expected_twist = theta*bp/(2*np.pi)
# print("expected twist is "+str(expected_twist))
    


# spline1 = pos_2_spline(strand1pos,circular,reverse=False)
# spline2 = pos_2_spline(strand2pos,circular,reverse=False)
# twist, writhe = get_twist_writhe(spline1, spline2, npoints, circular = True, integral_type = "simple")
# Lk = twist + writhe
# print('twist is '+ str(twist))
# print('writhe is '+ str(writhe))
# print('linking number is  '+ str(Lk))

# deltaLk = Lk - Lk0

    
# def plot(list_vecs,scatter=False):
#     vecs = np.array(list_vecs)
#     x = vecs[:,0]
#     y = vecs[:,1]
#     z = vecs[:,2]
#     if scatter==False:
#         return ax.plot(x,y,z)
#     else:
#         return ax.scatter(x,y,z)

    

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlim3d(-L/5, L/5)
# ax.set_ylim3d(-L/5, L/5)
# ax.set_zlim3d(-L/5, L/5)

# plot(strand1pos,scatter=False)
# plot(strand2pos,scatter=False)


# plt.show()

# box = 100


write_conf("nucleoid.dat",strand1pos,normals1,strand2pos,normals2,tan_vals,bp-1,box)
write_top("nucleoid.top",bp-1,sequence,circular)













