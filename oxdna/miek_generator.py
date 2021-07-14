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


def rot(B,A,theta):  #rotates vector B about axis A using Euler-Rodrigues formula
    return A*(np.dot(A,B))+np.cos(theta)*np.cross((np.cross(A, B)), A)+np.sin(theta)*(np.cross(A,B))


speedup = True

if speedup:
    from tools_fast_math import get_twist_writhe,multidet,normalize,rot_fast
else:
    from tools import multidet,normalize,get_twist_writhe,rot



'''
Parameters and checks

'''

bp = 500
circular = True         # circular boolean is true when structure is closed
periodic = True         # periodic boolean is true when x(n)=x(1), i.e. the final position is equal to the first
sequence_input = False  # sequence = True means a sequence file will be given
strands = 2
oxdna_unit_len = 0.8518 # nm
lpbp = 0.3897628551303122        # length per base pair in simulation units
dup_rad = 0.6 # duplex radius in simualtion units
L = bp*lpbp                         # contour length 
base_length = 0.7                   # distance from backbone to base
centre_line = []
pitch = 10.34                       # helical pitch of dsDNA
npoints = bp


Lk0 = bp/pitch # this is the relaxed linking number, deltaLk=0
deltaLk = 0
    
theta0 = Lk0*2*np.pi/bp

def gen_seq(bp):
    sequence = []
    if sequence_input == False:
        print('No sequence given, generating a random sequence')
        for count in range(bp):
           sequence.append(choice("CGTA"))
    return sequence

       
def generate(bp,centre_line,theta,circular):
    """
    Main function for generating the strand positions and normals for a duplex
    
    As we move around the centre-line, we rotate the bases counter clockwise
    by some angle theta (to form a right handed helix)
    
    v1 is a vector perpendicular to the tangent at the first base
    
    """
    
    strand1pos = []
    strand2pos = []
    normals1 = []
    normals2 = []
    
    
    # calculate the unit tangents at every base pair
    # for circular dna len(tan_vals)=bp, otherwise len(tan_vals)=bp-1
    tan_vals = []
    old_point = centre_line[-1]
    point = centre_line[0]
    diff = np.array([point[0]-old_point[0],point[1]-old_point[1],point[2]-old_point[2]])
    diff /=np.linalg.norm(diff)
    tan_vals.append(diff)
    old_point = centre_line[0]
    for i in range(1,bp):
        point = centre_line[i]
        diff = np.array([point[0]-old_point[0],point[1]-old_point[1],point[2]-old_point[2]])
        diff /=np.linalg.norm(diff)
        tan_vals.append(diff)
        old_point = point
        
    rad = dup_rad
    
    # we need a vector perpendicular to the inital tangent
    v1 = np.random.random_sample(3)
    v1 -= tan_vals[0] * (np.dot(tan_vals[0], v1))
    v1 /= np.linalg.norm(v1)

    rb1 = centre_line[0] + dup_rad*v1
    rb2 = centre_line[0] - dup_rad*v1
    normals1.append(-v1)
    normals2.append(v1)
      
    for i in range(1,len(tan_vals)):   
        # remove component of v2 in the direction of tangent
        v1 -= v1*(np.dot(v1,tan_vals[i]))
        v1 /= np.linalg.norm(v1)
        # rotate by theta and normalise
        v1 = rot(v1,tan_vals[i],theta)
        v1 /= np.linalg.norm(v1)
        
             
        rb1 = centre_line[i] + rad*v1
        rb2 = centre_line[i] - rad*v1
        
        normals1.append(-v1)
        normals2.append(v1)
        strand1pos.append(rb1)
        strand2pos.append(rb2)

    return strand1pos,strand2pos,normals1,normals2,tan_vals

def plot(list_vecs,ax,scatter=False):
    vecs = np.array(list_vecs)
    x = vecs[:,0]
    y = vecs[:,1]
    z = vecs[:,2]
    if scatter==False:
        return ax.plot(x,y,z)
    else:
        return ax.scatter(x,y,z)
    
def pos_2_spline(pos,circular,reverse = False):
    if reverse:
        pos.reverse()
    if circular:
        if reverse:
            pos.append(pos[-1])
        else:
            pos.append(pos[0])
    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in pos]
    yy = [vec[1] for vec in pos]
    zz = [vec[2] for vec in pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = splrep(range(len(xx)), xx, k = 3, s = 0, per = circular)
    spline_yy = splrep(range(len(yy)), yy, k = 3, s = 0, per = circular)
    spline_zz = splrep(range(len(zz)), zz, k = 3, s = 0, per = circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def evaluate_spline(spline,domain):
    vec_list = []
    splinex,spliney,splinez = spline[0]
    xx = splev(domain, splinex)
    yy = splev(domain, spliney)
    zz = splev(domain, splinez)
    for i in range(len(xx)):
        vec_list.append(np.array([xx[i],yy[i],zz[i]]))
    return vec_list
    


def write_top(top_file_name,bp,sequence,circular):
    
    top = open(top_file_name,"w")
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



def write_conf(conf_file_name,strand1pos,normals1,strand2pos,normals2,tan_vals,bp,box):
    conf = open(conf_file_name,"w")
    conf.write("t = 0" + "\n")
    conf.write("b = "+str(box)+" " +str(box)+" "+str(box) + "\n")
    conf.write("E = 0" + "\n")
    for s in range(strands):
        s = s + 1
        if s==1:
            for i in range(0,bp):
                conf.write(str(strand1pos[i][0]) + " " + str(strand1pos[i][1]) + " " + str(strand1pos[i][2])  \
                            +" " + str(normals1[i][0]) + " " + str(normals1[i][1]) + " " + str(normals1[i][2]) \
                                 + " " +str(tan_vals[i][0]) + " " + str(tan_vals[i][1]) + " " + str(tan_vals[i][2]) \
                                     + " 0 0 0 0 0 0" +"\n")
        else:
            for i in range(0,bp):
                conf.write(str(strand2pos[bp-i-1][0]) + " " + str(strand2pos[bp-i-1][1]) + " " + str(strand2pos[bp-i-1][2]) \
                           +" " + str(normals2[bp-i-1][0]) + " " + str(normals2[bp-i-1][1]) + " " + str(normals2[bp-i-1][2]) \
                               + " " + str(-tan_vals[bp-i-1][0]) + " " + str(-tan_vals[bp-i-1][1]) + " " + str(-tan_vals[bp-i-1][2]) \
                                   + " 0 0 0 0 0 0" +"\n")
    
    conf.close()

















