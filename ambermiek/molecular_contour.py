#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:41:05 2022

This script creates the molecular contour with removed intrisnic helicity, intakes a pdb and exports mol_contour.xyz
using WrLine's algorithm, fucntions Arctan360, setZ, setX, CAXIS and HAXIS are modified from caxislib.py

@author: michaelselby
"""


import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from numba import njit


sys.path.append('/Users/michaelselby/Documents/DPhil/miektils/ambermiek')
sys.path.append('/Users/michaelselby/Documents/DPhil/miektils')
sys.path.append('/Users/michaelselby/WrLINE')



@njit(fastmath=True)
def setZ(V):
    x = V[0]
    y = V[1]
    z = V[2]
    # defining rotation matrix about x,y axis and euler angles as: 
    # Rxy(A,B) = dot( Ry(B),Rx(A) )
    # Given r = (x,y,z) and uz = (0,0,1)
    # solve the equation uz = dot( Rxy, r )
    A = np.arctan2(y,z)
    B = np.arctan( -x / np.sqrt(y**2+z**2) )
    # calculating Rxy
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(A), -np.sin(A)], [0.0, np.sin(A), np.cos(A)]])
    Ry = np.array([[np.cos(B), 0.0, np.sin(B)], [0.0, 1.0, 0.0], [-np.sin(B), 0.0, np.cos(B)]])
    Rxy = np.dot( Ry, Rx )
    return Rxy

@njit(fastmath=True)
def setX(V):
    x = V[0]
    y = V[1]
    # z = V[2]
    # defining rotation matrix about z axis and euler angles as: 
    # Rz(C)
    C = - np.arctan2(y,x)
    Rz = np.array([[np.cos(C), -np.sin(C), 0.0], [np.sin(C), np.cos(C), 0.0], [0.0, 0.0, 1.0]])
    return Rz

@njit(fastmath=True)
def CAXIS(bp, r, theta_vals, linear=False): 
    r1 = np.zeros(np.shape(r)) 
    for j in range(bp):
        Tw = 0
        Sum = np.zeros((3))          # summation of coordinates
        Tw += theta_vals[j]        # Twist of the central bp step
        Sum[:] += r[j,:]
        k = 0
        while (Tw < 360.0):
            k += 1
            prev = Tw             # store previous total twist
            if linear and (j-k < 0 or j+k >= bp):
                break
            else:
                Tw += theta_vals[(j-k)%bp] + theta_vals[(j+k)%bp] # adding two more flanking steps then ttw would exceed 360.0
                Sum[:] += r[(j-k)%bp,:] + r[(j+k)%bp,:]  # sum up the single helix position
        
        if linear and (j-k < 0 or j+k >= bp):
            w = 0
        else:
            w = (360.0 - prev) / (Tw - prev)       # weighting
            Sum[:] -= (1-w)*(r[(j-k)%bp,:] + r[(j+k)%bp,:]) # now adding the flanks with the weight w<1
        
        W = np.array( [w, w, w] )      #  make it works for 3D
        r1[j,:] = Sum[:]/(2*(k+W)-1) 
        
    return r1

# average position of 2x5 neighbours of C1'-midpoints C1' 
# @njit(fastmath=True)
def HAXIS(bp, r, linear = False):
    haxis = np.zeros((bp,3))
    for j in range(bp):
        Sum = np.zeros((3))
        Sum[:] += r[j,:]
        k = 0
        while k < 5:
            k += 1
            if linear:
                try:
                    Sum[:] += (r[j-k, :] + r[j+k,:])
                except IndexError:
                    k -= 1
                    break
            else:       
                Sum[:] += r[(j-k)%bp,:] + r[(j+k)%bp,:]  # sum up the single helix position
        haxis[j,:] = Sum[:]/(2*k+1) # average helix (almost full helical turn)
    return haxis

@njit(fastmath=True)
def get_twist(diff1, diff2, z):
    # Z is the vector defining Z axis
    # the two vectors
    # rotate the vector Z to z axis
    Rxy = setZ(z)
    # rotate r1 and r2 along with Z using rotation matrix Rxy 
    diff1 = np.dot( Rxy, diff1 )
    diff2 = np.dot( Rxy, diff2 )
    # rotate r1 about z axis so that y becomes 0 
    Rz = setX(diff1)
    diff2 = np.dot(Rz, diff2 )
    # now calculating the twist
    return np.arctan2(diff2[1],diff2[0])*180.0/np.pi


# @njit(fastmath=True)
def full_twist(bp, haxis, diff, linear = False):
    """
    Calculates twist
    """
    twist = np.zeros(bp)
    for j in range(bp):
            # Linear special cases
        if linear:
            # Does haxis[:, :, j+1] exist?
            if np.shape(haxis)[0] > j + 1:
                # If haxis[j-1,:] doesn't exist,
                # approximate using half the range
                if j > 0:
                    z = haxis[j+1, :] - haxis[j-1, :]
                else:
                    z = 2 * (haxis[j+1, :] - haxis[j, :])
                diff1 = diff[j,:]
                diff2 = diff[j+1, :]
                twist[j] = get_twist(diff1, diff2, z)
            # If not, just return zero
            # This may not be the best way to handle this
            else:
                twist[j] = 0
        else:
            z = haxis[(j+1) % bp, :] - haxis[(j-1) % bp, :]
            diff1 = diff[j%bp, :]
            diff2 = diff[ (j+1)%bp, :]
            twist[j] = twist(diff1, diff2, z)
            
    return twist



def get_molecular_contour(bp, strandApos, strandBpos, linear=False):
    centres = (strandApos + strandBpos)/2
    diff = - strandApos + strandBpos
    diff /= np.sqrt((diff ** 2).sum(-1))[..., np.newaxis]
    
    if linear:
        new_centres = centres
    else:
        new_centres = np.array([ (centres[(i+1)%bp,:] + centres[i,:]) / 2 for i in range(bp)])

    haxis = HAXIS(bp, new_centres, linear=linear)

    twist_vals = full_twist(bp, haxis, diff, linear=linear)
    
    mol_cont = CAXIS(bp, new_centres, twist_vals, linear=linear)
    
    return mol_cont

     

     
 