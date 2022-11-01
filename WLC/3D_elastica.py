#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 13:41:44 2021

@author: michaelselby
"""

from numba import njit
from scipy.optimize import fsolve, newton, minimize
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tikzplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from tools import set_axes_equal, sn, E, am, K, Pi_inc, Pi
from tools_fast_math import multidet, norm, multidet2
import warnings
import copy
warnings.filterwarnings("ignore")


def rot(theta, A, B):  # rotates vector B about axis A using Euler-Rodrigues formula
    return A*(np.dot(A, B))+np.cos(theta)*np.cross((np.cross(A, B)), A)+np.sin(theta)*(np.cross(A, B))


# parameters
# KT is kBT in nmpN
# a1 is the bending modulus of the floppy segment with units
# L is the length of the floppy segment in nm
# a is the bending modulus of the rest of the curve
# f is the applied force in pN
kT = 4.114  # kT in pN nm at 298K
A = 45*kT  # bend persistence length in units of nm kbT
C = 92*kT
f = 1
d0 = 3  # radius of DNduplex for self-excluded volume interactions
m_spec = fsolve(lambda m: K(m)-2*E(np.pi/2,m),0.8261)[0]
u_spec = -1 + 2 * m_spec

source_pos = np.array([-20, -20, -20])


def plt_sphere(ax, fig, list_center, list_radius):
    for c, r in zip(list_center, list_radius):
        ax = fig.gca(projection='3d')

        # draw sphere
        u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
        x = r*np.cos(u)*np.sin(v)
        y = r*np.sin(u)*np.sin(v)
        z = r*np.cos(v)

        ax.plot_surface(x-c[0], y-c[1], z-c[2], color='yellow')


# define the rotation matrix R1
def rot1(theta):
    return np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
# defines the rotation matrix R3


def rot3(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
# defines the rotation R3(psi)*R1(theta)*R3(phi)


def euler_rot(psi, theta, phi):
    return np.dot(np.dot(rot3(psi), rot1(theta)), rot3(phi))


def get_vectors(vector0, psi, theta, phi):
    '''
    given an initial vector0, determines the evolution in vector determined
    by the euler angles psi,theta,phi

    '''
    val = np.array([np.dot(euler_rot(psi1, theta1, phi1), vector0)
                   for psi1, theta1, phi1 in zip(psi, theta, phi)])
    return val


def evolve_frame(triad0, psi, theta, phi):
    '''
    given an initial vector0, determines the evolution in vector determined
    by the euler angles psi,theta,phi

    '''
    normal0, binormal0, tangent0 = triad0
    mat = np.array([euler_rot(psi1, theta1, phi1).T for psi1,
                   theta1, phi1 in zip(psi, theta, phi)])
    npoints = np.shape(psi)[0]
    normal = np.array([mat[i, 0, 0]*normal0 + mat[i, 1, 0] *
                      binormal0 + mat[i, 2, 0]*tangent0 for i in range(npoints)])
    binormal = np.array([mat[i, 0, 1]*normal0 + mat[i, 1, 1] *
                        binormal0 + mat[i, 2, 1]*tangent0 for i in range(npoints)])
    tangent = np.array([mat[i, 0, 2]*normal0 + mat[i, 1, 2] *
                       binormal0 + mat[i, 2, 2]*tangent0 for i in range(npoints)])

    # bigtens = np.einsum("ijk,kl->ijl",mat,triad0)
    # normal = bigtens[:,0,:]
    # binormal = bigtens[:,1,:]
    # tangent = bigtens[:,2,:]

    return normal, binormal, tangent


def init_new_euler(new_basis, triad):
    """
    Given a triad of vectors, transform the basis to new_basis 
    
    Called when switching basis while wanting to preserve triad continuity

    """

    mat2 = np.dot(np.linalg.inv(new_basis), triad)
    
    phi = np.arctan(mat2[2,0]/mat2[2,1])
    theta = np.arccos(mat2[2,2])
    
    # getting psi is more diffuclt, must solve for cos(psi) using a quadratic
    a = np.cos(phi)
    b = - np.cos(theta)* np.sin(phi)
    c = triad[0,0]
    cospsi0plus = ( a * c + b * ( a**2 + b**2 - c**2 )**0.5 ) / ( a**2 + b**2 )
    # cospsi0minus = ( a * c - b * ( a**2 + b**2 - c**2 )**0.5 ) / ( a**2 + b**2 )
    psi = np.arccos(cospsi0plus)
    # psi = np.arccos(cospsi0minus)
    
    return theta, phi, psi
    
    


def concatenate_elastica(elas):

    for att_str in ["r", "tangent", "normal", "binormal"]:
        long = np.concatenate((getattr(elas.tail, att_str)[
                              ::-1, :], getattr(elas.braid, att_str)[::-1, :]), axis=0)
        long = np.concatenate(
            (long, getattr(elas.loop, att_str)[::-1, :]), axis=0)
        long = np.concatenate((long, getattr(elas.loop_ref, att_str)), axis=0)
        long = np.concatenate((long, getattr(elas.braid_ref, att_str)), axis=0)
        long = np.concatenate((long, getattr(elas.tail_ref, att_str)), axis=0)
        setattr(elas, att_str, long)


def invert_xz(vec):
    vec[:, 0] = -vec[:, 0]
    vec[:, 2] = -vec[:, 2]
    return vec


def invert_z(vec):
    vec[:, 2] = - vec[:, 2]

    return vec

def invert_xy(vec):
    vec[:, 0] = -vec[:, 0]
    vec[:, 1] = -vec[:, 1]
    return vec

def invert_y(vec):
    vec[:, 1] = -vec[:, 1]
    return vec

def y2z(vec):
    y = copy.deepcopy(vec[:, 1])
    vec[:, 1] = vec[:, 2]
    vec[:, 2] = y

def z2minusy(elastica):
    for vec in [elastica.r, elastica.tangent, elastica.normal, elastica.binormal]:
        y = copy.deepcopy(vec[:, 1])
        vec[:, 1] = vec[:, 2]
        vec[:, 2] = -y


def reflectxy(elastica):
    """
    Flips an elastica object in the y-axis

    """
    invert_xy(elastica.r)
    invert_z(elastica.tangent)
    invert_z(elastica.normal)
    invert_z(elastica.binormal)


def reflecty(elastica):
    invert_y(elastica.r)
    invert_y(elastica.tangent)
    invert_y(elastica.normal)
    invert_y(elastica.binormal)


def plot_surface(elastica, ax, source_pos, divisions=100):
    surfaces = []
    color_list = []
    npoints = elastica.npoints
    r = elastica.r
    tangent = elastica.tangent
    normal = elastica.normal
    binormal = elastica.binormal
    ratio = int(npoints/divisions)
    length = elastica.L
    for j in range(0, npoints):
        if j % ratio == 0:
            t, b, n = tangent[j], normal[j], binormal[j]
            x = r[j, 0]
            y = r[j, 1]
            z = r[j, 2]
            l = 2*length/(npoints/ratio)
            w = d0

            vec = np.array([x, y, z]) - source_pos
            vec /= np.linalg.norm(vec)
            dot = -np.dot(vec, n)

            X = [x-w*b[0]/2+l*t[0]/2, x-w*b[0]/2-l*t[0]/2,
                 x+w*b[0]/2-l*t[0]/2, x+w*b[0]/2+l*t[0]/2]
            Y = [y-w*b[1]/2+l*t[1]/2, y-w*b[1]/2-l*t[1]/2,
                 y+w*b[1]/2-l*t[1]/2, y+w*b[1]/2+l*t[1]/2]
            Z = [z-w*b[2]/2+l*t[2]/2, z-w*b[2]/2-l*t[2]/2,
                 z+w*b[2]/2-l*t[2]/2, z+w*b[2]/2+l*t[2]/2]
            verts = [list(zip(X, Y, Z))]
            surfaces.append(verts)
            # color = [color[i]*0.25*]
            if dot == 0:
                color_list.append([0, 0, 0])
            else:
                # color = [dot**gamma,0,0]
                color = [0, 1, 1]
                color_list.append(color)

    for i, surface in enumerate(surfaces):
        ax.add_collection3d(Poly3DCollection(
            surface, facecolors=color_list[i], edgecolors='k'))


def plot_tubes(elastica, ax, r0, divisions=100):
    """
    Parameters
    ----------
    elastica : Elastica3D or Solenoid3D
        elastica object created by 
    ax : matplotlib.axes._subplots.Axes3DSubplot
        3D axis supplied by elastica.plot() method
    divisions : int, optional
        The number of cylinders which comprise the tube. The default is 100.

    Returns
    -------
    None.

    """
    r = elastica.r
    ratio = int(elastica.npoints/divisions)
    length = elastica.L
    for j in range(0, elastica.npoints):
        if j % ratio == 0:
            t = elastica.tangent[j, :]
            n = elastica.normal[j, :]
            b = elastica.binormal[j, :]
            x = r[j, 0]
            y = r[j, 1]
            z = r[j, 2]
            l = 2*length/(elastica.npoints/ratio)
            l *= 1.01

            theta = np.linspace(0, 2*np.pi, 10)
            radius = r0
            s = np.linspace(-l/2, l/2, 10)

            thetas, ss = np.meshgrid(theta, s)
            x += radius*n[0]*np.cos(thetas)+b[0]*np.sin(thetas) + t[0]*ss
            y += radius*n[1]*np.cos(thetas)+b[1]*np.sin(thetas) + t[1]*ss
            z += radius*n[2]*np.cos(thetas)+b[2]*np.sin(thetas) + t[2]*ss

            ax.plot_surface(x, y, z, color='red')



def writhe_fuller(tan, tan_ref, writhe_ref):

    
    n = np.shape(tan)[0]
    writhe = 0
    for i in range(n-2):
        t1 = tan[i,:]
        t01 = tan_ref[i,:]
        t2 = tan[i+1,:]
        t02 = tan_ref[i,:]
        
        val = multidet(t01, t1, (t02+t2)-(t01+t1))/(1+np.dot(t1,t01))
        if np.isnan(val) or np.isinf(val):
            print("chumbo")
            print(i)
            print(t1,t01)
            pass
        else:
            writhe += val
        
    writhe /= (2*np.pi)
    writhe += writhe_ref
    return writhe

class Plectoneme3D:

    def __init__(self, length, A, C, f, r0, n, L_loop, pitch, p_twist):
        self.length = length
        self.L = self.length/2
        self.A = A
        self.C = C
        self.f = f
        self.r0 = r0
        self.n = n
        self.L_loop = L_loop
        self.pitch = pitch
        self.p_twist = p_twist
        self.neme_bool = True
        self.lam = (self.A/self.f)**0.5

    def generate_plectoneme(self, npoints):

        self.righthanded = True
    


        c_loop = -1
        a_loop, b_loop = loop_solver2(c_loop, self.L_loop, self.A, self.f, self.r0, niter=1000)
        m_loop = (b_loop-c_loop)/(a_loop-c_loop)
        self.loop = Elastica3D(a_loop, b_loop, c_loop, self.L_loop, self.A, self.C, self.f)
        self.npoints = npoints
        n_points_loop = round(self.npoints/2 * self.loop.L/self.L)
        
        origin = np.array([0, 0, 0])
        normal0, binormal0, tangent0 = np.array([1, 0, 0]), np.array([0, 0, -1]) ,np.array([0, 1, 0])
        self.loop.triad0 = np.array([normal0, binormal0, tangent0])

        self.loop.generate_euler_angles(n_points_loop)
        self.loop.euler2cartesian(origin)
        
        costhetaL = self.loop.uvals[-1] 
        thetaL = np.arccos(costhetaL)
        phiL, psiL = self.loop.phi[-1], self.loop.psi[-1]
        
        tangent0, binormal0 = np.array([0, 0, -1]), np.array([0, -1, 0])
        triad0 = np.array([normal0, binormal0, tangent0])
        triadloopL = np.array([self.loop.normal[-1,:],self.loop.binormal[-1,:],self.loop.tangent[-1,:]])

        
        theta0, phi0, psi0 = init_new_euler(triad0, triadloopL)
        b_braid = np.cos(theta0)
        
        self.braid = Solenoid3D(self.A, self.C, self.n,
                                b_braid, self.pitch, self.p_twist, self.righthanded)
        self.braid.triad0 = triad0
        self.phi_braid0 = phi0
        self.psi_braid0 = psi0
        n_points_braid = round(self.npoints/2 * self.braid.L/self.L)
        self.braid.generate_euler_angles(n_points_braid)
        
        self.braid.euler2cartesian(origin, phi0, psi0)
        self.loop.r -= np.array([0, 0, self.loop.r[-1, 2]])
        self.braid.r -= np.array([0, 0, self.braid.r[0, 2]])
        self.braid.r += np.array([self.loop.r[-1, 0], 0, 0])





        binormal0 = np.array([0, 0, -1])
        alpha = -0.1
        # for alpha in np.linspace(0,np.pi,100):
        tangent0 = np.array([np.sin(alpha), np.cos(alpha), 0])
        normal0 = np.array([np.cos(alpha), - np.sin(alpha), 0])
        
        tail_triad0 = np.array([normal0, binormal0, tangent0])
        
        triadbraidL = np.array([self.braid.normal[-1,:],self.braid.binormal[-1,:],self.braid.tangent[-1,:]])

        alphavals = [0]
        for alpha in alphavals:
            print(alpha)
            tangent0 = np.array([np.sin(alpha), np.cos(alpha), 0])
            normal0 = np.array([np.cos(alpha), - np.sin(alpha), 0])
            tail_triad0 = np.array([normal0, binormal0, tangent0])
            theta0, phi0, psi0 = init_new_euler(tail_triad0, triadbraidL)
            print(np.cos(theta0),phi0,psi0)
            
        b_tail = 1
        u0_tail = costhetaL
        a_tail = 1.015

        
        
        L_tail = self.L - self.loop.L - self.braid.L
        # if the tail is too small, plectoneme is forbidden
        if L_tail < 2:
            self.neme_bool = False
            
        c_tail = fsolve(lambda c: c+(b_tail-c)*sn(L_tail*((a_tail-c)/2)**0.5/self.lam+K((b_tail-c)/(a_tail-c)),(b_tail-c)/(a_tail-c))**2-u0_tail,-0.9)[0]
        print("ctail = ", c_tail)
        # a_tail = fsolve(lambda a: c_tail+(b_tail-c_tail)*sn(L_tail*((a-c_tail)/2)**0.5/self.lam+K((b_tail-c_tail)/(a-c_tail)),(b_tail-c_tail)/(a-c_tail))**2-costhetaL,1.1)[0]
            
        # a_tail = fsolve(lambda a: L_tail - (self.A/self.f)**0.5 *
        #                 (2/(a-c_tail))**0.5*K((b_tail-c_tail)/(a-c_tail)), 1.5)[0]
        
        # print(a_tail)


        self.tail = Elastica3D(a_tail, b_tail, c_tail, L_tail, self.A, self.C, self.f, increase=False)
        n_points_tail = round(self.npoints/2 * self.tail.L/self.L)
        self.tail.triad0 = tail_triad0
        self.tail.generate_euler_angles(n_points_tail)

        origin = self.braid.r[-1, :]

    
            

        self.tail.euler2cartesian(origin, 0, 1.9)
        
        print(np.dot(self.braid.tangent[-1,:],self.tail.tangent[0,:]))
        print(np.dot(self.braid.normal[-1,:],self.tail.normal[0,:]))
        print(np.dot(self.braid.binormal[-1,:],self.tail.binormal[0,:]))


        print("tangents")
        print(self.loop.tangent[0, :])
        print(self.loop.tangent[-1, :])
        print(self.braid.tangent[0, :])
        print(self.braid.tangent[-1, :])
        print(self.tail.tangent[0, :])
        print(self.tail.tangent[-1, :])


        min_z = self.braid.r[-1, 2]
        self.height = min_z
        self.elastica_list = [self.loop, self.braid, self.tail]
        for elastica in self.elastica_list:
            elastica.r[:, 2] -= min_z

        self.loop_ref = copy.deepcopy(self.loop)
        reflectxy(self.loop_ref)
        self.braid_ref = copy.deepcopy(self.braid)
        reflectxy(self.braid_ref)
        self.tail_ref = copy.deepcopy(self.tail)
        reflectxy(self.tail_ref)
        
        concatenate_elastica(self)
        

    def get_twist(self):
        self.twist = self.loop.get_twist() + self.braid.get_twist() + self.tail.get_twist()
        self.twist *= 2
        return self.twist

    
    @staticmethod
    @njit(fastmath=True)
    def get_writhe(tangent, r, npoints, length):

        writhe = 0
        for i in range(npoints):
            ti = tangent[i, :]
            ri = r[i, :]
            for j in range(npoints):
                if i > j:
                    tj = tangent[j, :]
                    rj = r[j, :]
                    val = multidet(ti, tj, ri-rj)
                    diff = norm(ri-rj)
                    if diff < 0.001:
                        pass
                    else:
                        val /= diff**3
                    if np.isnan(val):
                        pass
                    else:
                        writhe += val
        writhe /= (2*np.pi)
        writhe *= (length/npoints)**2
        
        return writhe

    def get_link(self):
        self.link = self.get_twist() + self.get_writhe(self.tangent,self.r,self.npoints,self.length)
        return self.link

    def plot_plectoneme(self, style= None):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_ylabel("y")
        ax.set_xlabel("x")
        ax.set_zlabel("z")
        self.height = 60

        ax.set_xlim(-self.height/2, self.height/2)
        ax.set_ylim(-self.height/2, self.height/2)
        ax.set_zlim(0, self.height)


        if style == "tubes":
            plot_tubes(self, ax, self.r0, divisions=50)

        elif style == "ribbon":
            plot_surface(self, ax, source_pos, divisions=100)
            
        else:
            ax.plot(*self.r.T, color="red")

        plt.show()
        
        
    def get_energy(self):
        self.energy = self.loop.get_energy() + self.braid.get_energy() + self.tail.get_energy()
        self.energy *= 2
        return self.energy

    def printvals(self):
        print("Overall:")
        print("Total Length = " + str(self.length))
        print(" ")
        print("Loop:")
        self.loop.printvals()
        print("")
        print("Braid:")
        self.braid.printvals()
        print("")
        print("Tail:")
        self.tail.printvals()
        print("")
        print("Total Twist = " + str(self.get_twist()))
        print("Writhe = " + str(self.get_writhe(self.tangent, self.r, self.npoints, self.length)))
        print("Linking Number = " + str(self.get_link()))


class Solenoid3D:
    def __init__(self, A, C, n, u0, pitch, p_twist, righthanded=True):
        self.A = A
        self.C = C
        self.n = n
        self.u0 = u0
        self.pitch = pitch
        self.rad = self.pitch*np.tan(np.arccos(self.u0))
        self.L = 2*np.pi*self.n*(self.pitch**2+self.rad**2)**0.5
        self.p_psi = p_twist
        self.righthanded = righthanded
        if righthanded:
            self.factor = 1
        else:
            self.factor = -1

    def generate_euler_angles(self, npoints):
        self.npoints = npoints

        self.uvals = np.full(self.npoints, self.u0)
        self.theta = np.full(self.npoints, np.arccos(self.u0))
        self.s = np.linspace(0, self.L, self.npoints)
        self.ds = self.L/self.npoints
        self.phi = 2*np.pi*self.n*self.s/self.L
        self.psi = (self.p_psi/self.C - 2*np.pi*self.n*self.u0/self.L)*self.s

    def euler2cartesian(self, origin, phi0, psi0):
        self.phi += phi0
        self.psi += psi0
        self.normal, self.binormal, self.tangent = evolve_frame(
            self.triad0, self.psi, self.theta, self.phi)
        self.r = np.cumsum(self.tangent, axis=0)*self.ds
        if self.righthanded:
            reflecty(self)
        self.r += origin

    def plot_solenoid(self, source_pos, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(*self.r.T, color='red')
        ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
        set_axes_equal(ax)
        plot_surface(self.r, self.tangent, self.normal, self.binormal,
                     ax, source_pos, self.npoints, divisions=divisions)
        plt.show()

    def get_twist(self):
        self.twist = self.p_psi*self.L/(2*np.pi*self.C)
        return self.twist
        
    def get_energy(self):
        self.energy = 2*self.A*np.pi**2*(1-self.u0**2)/self.L + 2*np.pi**2*self.C*(self.get_twist())**2/self.L
        return self.energy

    def printvals(self):
        print("Length = " + str(self.L))
        print("Number of braids = " + str(self.n))
        print("Radius = " + str(self.rad))
        print("Pitch = " + str(self.pitch))
        print("Twist = " + str(self.get_twist()))


class Elastica3D:

    def __init__(self, a, b, c, L, A, C, f, increase = True):
        self.a = a
        self.b = b
        self.c = c
        self.A = A
        self.C = C
        self.f = f
        self.lam = (self.A/self.f)**0.5
        self.m = (self.b-self.c)/(self.a-self.c)
        self.L = L
        self.increase = increase
        p1 = (2*self.A*self.f*(self.a+1)*(self.b+1)*(self.c+1))**0.5
        p2 = (2*self.A*self.f*(self.a-1)*(self.b-1)*(self.c-1))**0.5

        self.p_psi = (p1 + p2)/2
        self.p_phi = (p1 - p2)/2

    def generate_euler_angles(self, npoints):

        self.npoints = npoints
        self.s = np.linspace(0, self.L, self.npoints)
        self.ds = self.L/self.npoints

        self.uvals = costheta_incomp(self.a, self.b, self.c, self.L,
                              self.A, self.f, self.npoints, increase = self.increase)
        self.theta = np.arccos(self.uvals)

        if np.isclose(self.p_psi, self.p_phi):
            phidot = self.p_phi / (self.A*(1 + self.uvals))
            psidot = self.p_psi*(1/self.C-1/self.A) + phidot
        elif np.isclose(self.p_psi, -self.p_phi):
            phidot = -self.p_psi / (self.A*(1 - self.uvals))
            psidot = self.p_psi*(1/self.C-1/self.A) - phidot
            
        else:
            phidot = (self.p_phi - self.p_psi*self.uvals) / \
                (self.A*(1-self.uvals**2))
            psidot = self.p_psi*(1/self.C-1/self.A) + (self.p_psi - self.p_phi *
                                                       self.uvals)/(self.A*(1-self.uvals**2))
        phidot[0] = 0
        psidot[0] = 0
        ds = self.s[1] - self.s[0]
        self.phi = np.cumsum(phidot)*ds
        self.psi = np.cumsum(psidot)*ds
        

    def euler2cartesian(self, origin, phi0=0, psi0=0):
        self.phi += phi0
        self.psi += psi0

        self.normal, self.binormal, self.tangent = evolve_frame(
            self.triad0, self.psi, self.theta, self.phi)
        self.r = np.cumsum(self.tangent, axis=0)*self.ds
        self.r += origin
        self.r = np.insert(self.r, 0, origin, axis=0)

    def plot_elastica(self, source_pos, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(*self.r.T, color='red')

        ax.set_ylabel("y")
        ax.set_xlabel("x")
        ax.set_zlabel("z")
        set_axes_equal(ax)
        plot_surface(self, ax, source_pos, divisions)
        plt.show()

    def get_twist(self):
        self.twist = self.L*self.p_psi/(2*np.pi*self.C)
        return self.twist
    
    def get_energy(self):
        self.energy = self.f * (self.a + self.b + self.c) * self.L
        self.energy += self.p_psi**2 * (1/self.C - 1/self.A) * self.L / 2
        self.energy += 2*(2*self.A*self.f*(self.a-self.c))**0.5*E(np.pi/2,self.m) - 2*self.a*self.f*self.L
        return self.energy

        

    def printvals(self):
        print("L = " + str(self.L))
        print("a = " + str(self.a))
        print("b = " + str(self.b))
        print("c = " + str(self.c))
        print("m = " + str(self.m))
        print("Twist = " + str(self.get_twist()))


def costheta(a, b, c, A, f, npoints):
    """
    Parameters
    ----------
    a,b,c : float
        three roots which describe elastica, ordered c<=b<=a
    A: float
        Bending modulus
    f : float
        force
    npoints : int
        number of evaluation points

    Returns
    -------
    costheta : numpy array of size (npoints)
        values of cosine of theta euler angle along an elastica
    """
    m = (b-c)/(a-c)
    lam = (A/f)**0.5
    L = lam*(2/(a-c))**0.5*K(m)
    if np.isclose(b, c):
        costheta = np.full(npoints, b)
    else:
        s = np.linspace(0, L, npoints)
        costheta = c + (b-c)*sn(s*((a-c)/2)**0.5/lam, m)**2

    costheta[0] = c

    return costheta


def costheta_incomp(a, b, c, L, A, f, npoints, increase=True):

    m = (b-c)/(a-c)
    lam = (A/f)**0.5
    if np.isclose(b, c):
        costheta = np.full(npoints, b)
    elif increase:
        s = np.linspace(0, L, npoints)
        costheta = c + (b-c)*sn(s*((a-c)/2)**0.5/lam, m)**2
    else:
        s = np.linspace(0, L, npoints)
        s = s[::-1]
        costheta = c + (b-c)*sn(s*((a-c)/2)**0.5/lam + K(m), m)**2

    return costheta


def phi_func(a, b, c, A, f, npoints):
    m = (b-c)/(a-c)
    lam = (A/f)**0.5
    L = lam*(2/(a-c))**0.5*K(m)
    p1 = (2*A*f*(a+1)*(b+1)*(c+1))**0.5
    p2 = (2*A*f*(a-1)*(b-1)*(c-1))**0.5
    p_psi = (p1 + p2)/2
    p_phi = (p1 - p2)/2

    uvals = costheta(a, b, c, A, f, npoints)

    if b == 1:
        # np.array([(p_phi)/(A*(1+u)) for u in uvals])
        phidot = (p_phi)/(A*(1+uvals))
    elif c == -1:
        phidot = p_phi/(A*(1-uvals))
    else:
        phidot = (p_phi - p_psi*uvals)/(A*(1-uvals**2))
    ds = L/npoints
    phi = np.cumsum(phidot)*ds

    return phi

def phi_func_incomp(a, b, c, L, A, f, npoints):
    m = (b-c)/(a-c)
    lam = (A/f)**0.5
    p1 = (2*A*f*(a+1)*(b+1)*(c+1))**0.5
    p2 = (2*A*f*(a-1)*(b-1)*(c-1))**0.5
    p_psi = (p1 + p2)/2
    p_phi = (p1 - p2)/2

    uvals = costheta_incomp(a, b, c, L, A, f, npoints)

    if b == 1:
        # np.array([(p_phi)/(A*(1+u)) for u in uvals])
        phidot = (p_phi)/(A*(1+uvals))
    elif c == -1:
        phidot = p_phi/(A*(1-uvals))
    else:
        phidot = (p_phi - p_psi*uvals)/(A*(1-uvals**2))
    ds = L/npoints
    phi = np.cumsum(phidot)*ds

    return phi


def phi_exact(a, b, c, A, f, npoints):
    m = (b-c)/(a-c)
    lam = (A/f)**0.5
    k = m**0.5
    L = lam*k*(2/(b-c))**0.5*K(m)
    p1 = (2*A*f*(a+1)*(b+1)*(c+1))**0.5
    p2 = (2*A*f*(a-1)*(b-1)*(c-1))**0.5
    p_psi = (p1 + p2)/2
    p_phi = (p1 - p2)/2
    s = np.linspace(0, L, npoints)
    angles = am(s*((a-c)/2)**0.5/lam, m)
    phi = np.array([((p_psi+p_phi)/(1+c)*Pi_inc((c-b)/(c+1), ang, m) +
                   (p_phi-p_psi)/(1-c)*Pi_inc((b-c)/(1-c), ang, m)) for ang in angles])
    phi /= ((a-c)**0.5*(2*A*f)**0.5)

    return phi


def zclosure(a, b, c):
    m = (b-c)/(a-c)
    val = a*K(m) - (a-c)*E(np.pi/2, m)
    return val


def rfunc(a, b, c, A, f, triad0, npoints):
    e1, e2, e3 = triad0
    lam = (A/f)**0.5
    m = (b-c)/(a-c)
    L = lam*(2/(a-c))**0.5*K(m)
    ds = L/npoints
    uvals = costheta(a, b, c, A, f, npoints)
    sintheta = (1-uvals**2)**0.5
    phi = phi_func(a, b, c, A, f, npoints)
    x = np.cumsum(sintheta*np.sin(phi))*ds
    y = np.cumsum(sintheta*np.cos(phi))*ds
    z = np.cumsum(uvals)*ds
    r = x*e1 + y*e2 + z*e3
    return r


def rend(a, b, c, A, f, triad0, npoints):
    e1, e2, e3 = triad0
    lam = (A/f)**0.5
    m = (b-c)/(a-c)
    L = lam*(2/(a-c))**0.5*K(m)
    ds = L/npoints
    uvals = costheta(a, b, c, A, f, npoints)
    sintheta = (1-uvals**2)**0.5
    phi = phi_func(a, b, c, A, f, npoints)
    x = np.sum(sintheta*np.sin(phi))*ds
    y = np.sum(sintheta*np.cos(phi))*ds
    z = np.sum(uvals)*ds
    r = x*e1 + y*e2 + z*e3
    return r


def xfunc(a, b, c, A, f, npoints):
    lam = (A/f)**0.5
    m = (b-c)/(a-c)
    L = lam*(2/(a-c))**0.5*K(m)
    ds = L/npoints
    uvals = costheta(a, b, c, A, f, npoints)
    sintheta = (1-uvals**2)**0.5
    phi = phi_func(a, b, c, A, f, npoints)
    val = np.cumsum(sintheta*np.sin(phi))*ds
    return val


def yfunc(a, b, c, A, f, npoints):
    lam = (A/f)**0.5
    m = (b-c)/(a-c)
    k = m**0.5
    L = lam*k*(2/(b-c))**0.5*K(m)
    ds = L/npoints
    uvals = costheta(a, b, c, A, f, npoints)
    sintheta = (1-uvals**2)**0.5
    phi = phi_func(a, b, c, A, f, npoints)
    val = np.cumsum(sintheta*np.cos(phi))*ds
    return val


def zfunc(a, b, c, A, f, npoints):
    lam = (A/f)**0.5
    m = (b-c)/(a-c)
    k = m**0.5
    L = lam*k*(2/(b-c))**0.5*K(m)
    ds = L/npoints
    uvals = costheta(a, b, c, A, f, npoints)
    val = np.cumsum(uvals)*ds
    return val


def xend(a, b, c, A, f, npoints):
    lam = (A/f)**0.5
    m = (b-c)/(a-c)
    L = lam*(2/(a-c))**0.5*K(m)
    ds = L/npoints
    uvals = costheta(a, b, c, A, f, npoints)
    sintheta = (1-uvals**2)**0.5
    phi = phi_func(a, b, c, A, f, npoints)
    val = np.sum(sintheta*np.sin(phi))*ds
    return val

def xend_incomp(a, b, c, L, A, f, npoints):
    ds = L/npoints
    uvals = costheta_incomp(a, b, c, L, A, f, npoints)
    sintheta = (1-uvals**2)**0.5
    phi = phi_func_incomp(a, b, c, L, A, f, npoints)
    xend = np.sum(sintheta*np.sin(phi))*ds
    return xend

def zend_incomp(a, c, m, L, A, f, npoints):
    lam = (A/f)**0.5
    ds = L/npoints
    s = np.linspace(0, L, npoints)
    uvals = c + m*(a-c)*sn(s*((a-c)/2)**0.5/lam, m)**2
    zend = np.sum(uvals)*ds
    return zend


def mfunc(m, c, A, f, npoints):

    rat = K(m)/E(np.pi/2, m)
    a = c/(1-rat)
    b = m*(a-c)+c
    return xend(a, b, c, A, f, npoints)


def tail_solver(b, c, A, f, npoints, r0):
    a = newton(lambda a: xend(a, b, c, A, f, npoints)-r0, 1.2)
    return a


def loop_func(a, c, A, f, r0, npoints):
    m = fsolve(lambda m: a*K(m)-(a-c)*E(np.pi/2, m), 0.8)[0]
    b = c + (a-c)*m
    return xend(a, b, c, A, f, npoints)

def loop_func2(a, c, L, A, f, r0, npoints):
    lam = (A/f)**0.5
    m = fsolve(lambda m: a*L - lam*(2*(a-c))**0.5*E(am(L*((a-c)/2)**0.5/lam, m), m), 0.8)[0]
    # m = fsolve(lambda m: zend_incomp(a, c, m, L, A, f, 1000),0.8)[0]
    b = c + (a-c)*m
    return xend_incomp(a, b, c, L, A, f, npoints) 



def braid_func(a, c, A, f, r0, npoints):
    m = fsolve(lambda m: a*K(m)-(a-c)*E(np.pi/2, m), 0.0001)[0]
    b = c + (a-c)*m
    return xend(a, b, c, A, f, npoints)


# loop solver considerably quicker than the alternatives, depsite being able to elimate variables by hand
def loop_solver(c, A, f, r0, niter):
    a = newton(lambda a: loop_func(a, c, A, f, r0, niter)+r0, 1.0001)
    m = fsolve(lambda m: a*K(m)-(a-c)*E(np.pi/2, m), 0.8)[0]
    b = c + (a-c)*m
    return a, b


def loop_solver2(c, L, A, f, r0, niter):
    lam = (A/f)**0.5
    a = newton(lambda a: loop_func2(a, c, L, A, f, r0, niter)+r0, 1.0001)
    # m = fsolve(lambda m: zend_incomp(a, c, m, L, A, f, niter),0.8)[0]
    m = fsolve(lambda m: a*L - lam*(2*(a-c))**0.5*E(am(L*((a-c)/2)**0.5/lam, m), m), 0.8)[0]
    # print(a*L - lam*(2*(a-c))**0.5*E(am(L*((a-c)/2)**0.5/lam,m),m))
    b = c + (a-c)*m
    return a, b

def braid_solver(c, A, f, r0, npoints):
    a = newton(lambda a: loop_func(a, c, A, f, r0, npoints)-2*r0, 2)
    m = fsolve(lambda m: a*K(m)-(a-c)*E(np.pi/2, m), 0.01)[0]
    b = c + (a-c)*m
    return a, b


nemey_daniel = Plectoneme3D(200, A, C, f, 1.5, 1, 35, 3, 0.01)
nemey_daniel.generate_plectoneme(4000)
# nemey_daniel.generate_coordinates(4000)
nemey_daniel.plot_plectoneme(style="ribbon")
# nemey_daniel.printvals()


# origin = np.array([0, 0, 0])
# tangent0 = np.array([0, 1, 0])
# normal0 = np.array([1, 0, 0])
# binormal0 = np.cross(tangent0, normal0)
# triad0 = np.array([normal0, binormal0, tangent0])
# test_loop = Elastica3D(1.01,-0.221,-1,A,C,f)
# test_loop.generate_euler_angles(1000)
# test_loop.euler2cartesian(triad0, origin)
# test_loop.plot_elastica(source_pos,20)


r0 = 1.5
L = 200
bp = L/0.34
pitch_bdna = 10.5
Lk0 = bp/pitch_bdna
sigma = -0.02
deltaLk = sigma*Lk0
nlist = [1,2,3]
tol = 0.1

print("--------------------------------------------------------")
print("Neme Factory")
print("")

best_neme_list = []

# lowest_energy = 100000
# first = True
# L_loop = 30
# for n in nlist:
#     best_neme = None
#     for p in np.linspace(1, 10, 20):
#         neme = Plectoneme3D(L, A, C, f, r0, n, L_loop, p, 0)
#         neme.generate_plectoneme()
#         if neme.neme_bool:
#             neme.generate_coordinates(1000)
#             twist = neme.get_twist()
#             writhe = neme.get_writhe(neme.tangent,neme.r,neme.npoints,neme.length)
#             extra_twist = deltaLk - writhe 
#             tail_twist = neme.tail.get_twist()
#             loop_twist = neme.loop.get_twist()
#             extra_twist /= 2
#             extra_twist -= (tail_twist+loop_twist)
#             # we guess that the twist required to correct for sigma will be mostly controlled 
#             # by the braid, so p_twist_guess must produce half the desired twist
#             p_twist_guess = 2*extra_twist*np.pi*neme.C/(neme.braid.L)

        
#         neme = Plectoneme3D(L, A, C, f, r0, n, L_loop, p, p_twist_guess)
#         neme.generate_plectoneme()
#         if neme.neme_bool:
#             neme.generate_coordinates(1000)
#             twist = neme.get_twist()
#             twist_braid = neme.braid.get_twist()
#             writhe = neme.get_writhe(neme.tangent,neme.r,neme.npoints,neme.length)
                
                
                
#                 # if first:
#                 #     writhe = neme.get_writhe(neme.tangent,neme.r,neme.npoints,neme.length)
    
#                 #     print(writhe)
#                 #     writhe_ref = writhe
#                 #     tan_ref = neme.tangent
#                 #     first = False
#                 # else:
    
#                 #     writhe = neme.get_writhe(neme.tangent,neme.r,neme.npoints,neme.length)
#                 #     tan = neme.tangent
#                 #     writhe2 = writhe_fuller(tan, tan_ref, writhe_ref)
#                 #     print("writhe:", writhe)
#                 #     print("fuller writhe:", writhe2)
#                 #     writhe_ref = writhe
#                 #     tan_ref = neme.tangent
                    
#             Lk = neme.get_link()
#             if (Lk - deltaLk)**2 < tol:
#                 print("neme wrangled!")
#                 energy = neme.get_energy()
#                 if energy < lowest_energy:
#                     best_neme = neme
#                     lowest_energy = energy
                    

                
#     if best_neme is not None:          
#         best_neme_list.append(best_neme) 
#     else:
#         print("this neme has been forsaken")

# for neme in best_neme_list:     
#     neme.plot_plectoneme(style="ribbon") 
 
        
# lowest_energy = 100000
# for neme in best_neme_list:
#     energy = neme.energy
#     if energy < lowest_energy:
#         bestest_neme = neme
#         lowest_energy = energy
            
# print("The bestest of all the nemes is the neme with " + str(bestest_neme.n) + " braids :-)")
            

