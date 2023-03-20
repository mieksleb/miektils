#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 14:58:15 2023

@author: michaelselby
"""
import sys
import numpy as np


conf = sys.argv[1]
top = sys.argv[2]

# conf = "linear/linear_128/MDNA/relax_trajectory_dv128t.dat"
# top = "linear/linear_128/MDNA/dv128t.top"

# conf = "linear/linear_128/MTTD/relax_trajectory_dv128tt.dat"
# top = "linear/linear_128/MTTD/dv128tt.top"


toplines = open(top,"r").readlines()

line0 = toplines[0].split()
n = int(line0[0])
strands = int(line0[1])

bp = int(n/strands)

if int(toplines[1].split()[2]) == -1:
    circular = False
else:
    circular = True
    
    
# def rot(theta, axis, vec):  #rotates vector B about axis A using Euler-Rodrigues formula
#     return vec*np.cos(theta) + np.sin(theta)*np.cross(axis, vec) + axis*(np.dot(axis,vec))*(1-np.cos(theta))
    
def rot(omega, A, B):  #rotates vector B about axis A using Euler-Rodrigues formula
    return A*(np.dot(A,B))+np.cos(omega)*np.cross((np.cross(A, B)), A)+np.sin(omega)*(np.cross(A,B))

def bsp(triad1,triad2):
    """
    
    Calculates the 6 base-step parameters: roll,tilt,twist,slide,shift,rise
    

    Parameters
    ----------
    triad1 : list of length 4, triad[i] is a numpy array of length 3
        triad for base 1
        triad1[0] is the position of the origin, triad1[1],triad1[2],triad1[2] are the x,y,z axes for the ring respectively
    triad2 : list of length 4, triad[i] is a numpy array of length 3
        triad for base 2
        triad2[0] is the position of the origin, triad2[1],triad2[2],triad2[2] are the x,y,z axes for the ring respectively

    Returns
    -------
    roll,tilt,twist,slide,shift,rise

    """
    r01,x1,y1,z1 = triad1
    r02,x2,y2,z2 = triad2
    
    gamma = np.arccos(np.dot(z1,z2))
    q = np.cross(z1,z2)
    q /= np.linalg.norm(q)
    
    r0mst = (r01+r02)/2
    x1prime,y1prime,z1prime = rot(gamma/2,q,x1),rot(gamma/2,q,y1),rot(gamma/2,q,z1)
    x2prime,y2prime,z2prime = rot(-gamma/2,q,x2),rot(-gamma/2,q,y2),rot(-gamma/2,q,z2)
    x1prime /= np.linalg.norm(x1prime)
    y1prime /= np.linalg.norm(y1prime)
    z1prime /= np.linalg.norm(z1prime)
    x2prime /= np.linalg.norm(x2prime) 
    y2prime /= np.linalg.norm(y2prime)
    z2prime /= np.linalg.norm(z2prime)
    xmst,ymst,zmst = (x1prime+x2prime)/2,(y1prime+y2prime)/2,(z1prime+z2prime)/2
    xmst /= np.linalg.norm(xmst)
    ymst /= np.linalg.norm(ymst)
    zmst /= np.linalg.norm(zmst)
    twist = np.arccos(np.dot(y1prime,y2prime))
    sign = np.dot(np.cross(y1prime,y2prime),zmst)
    if sign > 0:
        twist = abs(twist)
    else:
        twist = -abs(twist)
        
    twist *= 180/np.pi
    gamma *= 180/np.pi
        
    phi = np.arccos(np.dot(q,ymst))
    sign = np.dot(np.cross(q,ymst),zmst)
    if sign>0:
        phi = abs(phi)
    else:
        phi = -abs(phi)
        
    # phi must be bwteen -180 and 180
        
    shift = 0
    rise = 0
    slide = 0
    roll = gamma*np.cos(phi)
    tilt = gamma*np.sin(phi)
    
    diff = r02 - r01
    shift,slide,rise = np.dot(diff,xmst),np.dot(diff,ymst),np.dot(diff,zmst)
    return roll,tilt,twist,slide,shift,rise

def get_triad(r1,r2,n1,n2,b1,b2):
    n2 *= -1
    r0 = (r1 + r2) / 2
    z = (r2-r1)/np.linalg.norm((r2-r1))
    y = (n1 + n2)/2
    y /= np.linalg.norm(y)
    x = np.cross(y,z)
    x /= np.linalg.norm(x)
    return r0,x,y,z
    

r1 = np.zeros((bp,3))
r2 = np.zeros((bp,3))
n1 = np.zeros((bp,3))
b1 = np.zeros((bp,3))
n2 = np.zeros((bp,3))
b2 = np.zeros((bp,3))
step = 0

ang = np.zeros((bp))

with open(conf,"r") as conf:
    count = 0
    for line in conf:
        split_line = line.split()
        if split_line[0]=="t" or split_line[0]=="E" or split_line[0]=="b":
            continue
        else:
            x,y,z = float(split_line[0]),float(split_line[1]),float(split_line[2])
            bx,by,bz = float(split_line[3]),float(split_line[4]),float(split_line[5])
            nx,ny,nz = float(split_line[6]),float(split_line[7]),float(split_line[8])
            if count < bp:
                r1[count,:] = np.array([x,y,z])
                n1[count,:] = np.array([nx,ny,nz])
                b1[count,:] = np.array([bx,by,bz])
            else:
                # print("cheg")
                r2[-(count-bp)-1,:] = np.array([x,y,z])
                n2[-(count-bp)-1,:] = np.array([nx,ny,nz])
                b2[-(count-bp)-1,:] = np.array([bx,by,bz])
                
            if count == 2*bp -1:
                for i in range(bp-1):
                    # print(np.dot(n1[i,:],n2[i,:]))
                    triad1 = get_triad(r1[i,:], r2[i,:], n1[i,:], n2[i,:], b1[i,:], b2[i,:])
                    triad2 = get_triad(r1[i+1,:], r2[i+1,:], n1[i+1,:], n2[i+1,:], b1[i+1,:], b2[i+1,:])
                    # triad2 = get_triad(r1[bp-i-1,:], r2[bp-i-1,:], n1[bp-i-1,:], n2[bp-i-1,:], b1[bp-i-1,:], b2[bp-i-1,:])
                    roll,tilt,twist,slide,shift,rise = bsp(triad1,triad2)
                    # print(twist)
                    # print(np.linalg.norm(triad1[0]-triad2[0]))
                    ang[i] += 180*np.arccos(np.dot(triad1[-1],triad2[-1]))/np.pi
                    

                step += 1
                count = 0
                # if step > 1:
                #     break
                if step % 1000 == 0:
                    print(step)

            else:  
                count +=1
                    
    
ang /= step

with open("twist_bps.dat","w") as file:
    for angle in ang:
        file.write(str(angle) + "\n")

