"""
Created on Mon Jun 14 17:07:55 2021

Miek_generate

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
from miek_generator import generate,write_conf,write_top,gen_seq,pos_2_spline


def rot(B,A,theta):  #rotates vector B about axis A using Euler-Rodrigues formula
    return A*(np.dot(A,B))+np.cos(theta)*np.cross((np.cross(A, B)), A)+np.sin(theta)*(np.cross(A,B))

speedup = True

if speedup:
    from tools_fast_math import get_twist_writhe,multidet,normalize,rot_fast
else:
    from tools import multidet,normalize,get_twist_writhe,rot




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
bp = 1000
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

if circular == False:
    if periodic == False:
        print('Configration cannot be periodic but not circular!')
        
sequence = gen_seq(bp)


if periodic == True:
    length = bp+1
else:
    length = bp


# create the arclength range
svals = np.linspace(0,L,length)
ds = L/bp



'''
Constructing your centre-line

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

box = 5*max(xrange,yrange,zrange)

dist = np.linalg.norm(centre_line[0]-centre_line[1])

# we now centre the centre_line in the box

com = []
xcom=0
ycom=0
zcom=0
for i in range(len(xx)):
    xcom += xx[i]
    ycom += yy[i]
    zcom += zz[i]
xcom /= len(xx)
ycom /= len(yy)
zcom /= len(zz)

centre_line = []
for i in range(len(xx)):
    centre_line.append(np.array([xx[i]-xcom,yy[i]-ycom,zz[i]-zcom]))



# if the structure is circular and periodic, then we remove the last element of the centre_lines
if periodic == True:
    del centre_line[-1]
    

phi = 0

'''
This part of the code ensures that the configuration has the correct linking number
'''
print("Correcting the Linking number")
for j in range(0,2):

    theta = theta0 + phi
    
    strand1pos,strand2pos,normals1,normals2,tan_vals = generate(bp,centre_line,theta,circular)
    spline1 = pos_2_spline(strand1pos,circular,reverse=False)
    spline2 = pos_2_spline(strand2pos,circular,reverse=False)
    twist, writhe = get_twist_writhe(spline1, spline2, npoints, circular, integral_type = "simple")
    Lk = twist + writhe
    print('twist is '+ str(twist))
    print('writhe is '+ str(writhe))
    print('linking number is  '+ str(Lk))
    
    deltaLk_actual = Lk - Lk0
    difference = deltaLk_actual-deltaLk
    phi = difference*2*np.pi/bp
    


# The structure now has the correct linking number/ We now need to check the
# angle between the last base-pair and the first
v10 = normals1[0]
v2 = normals2[-1]
v2 -= v2*(np.dot(v2,tan_vals[-1]))
v2 /= np.linalg.norm(v2)
last_ang = np.arccos(np.dot(v2,v10))

delta_ang = last_ang - theta
delta_ang_pbp = delta_ang/bp # change in angle per base pair
# theta += delta_ang_pbp

expected_twist = theta*bp/(2*np.pi)
print("expected twist is "+str(expected_twist))
    


spline1 = pos_2_spline(strand1pos,circular,reverse=False)
spline2 = pos_2_spline(strand2pos,circular,reverse=False)
twist, writhe = get_twist_writhe(spline1, spline2, npoints, circular = True, integral_type = "simple")
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
ax.set_xlim3d(-L/5, L/5)
ax.set_ylim3d(-L/5, L/5)
ax.set_zlim3d(-L/5, L/5)

plot(strand1pos,scatter=False)
plot(strand2pos,scatter=False)


plt.show()

box = 100


write_conf("conf.dat",strand1pos,normals1,strand2pos,normals2,tan_vals,bp,box)
write_top("top.top",bp,sequence,circular)













