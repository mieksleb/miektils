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
from scipy.special import kn
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from tools import set_axes_equal, sn, cn, dn, E, am, K, Pi_inc, Pi, F
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
d0 = 3  # radius of DNAduplex for self-excluded volume interactions

parms = {}
parms[10] = 3, 1.97, 0.7
parms[50] = 1.34, 4.33, 0.7
parms[100] = 0.95, 6.24, 0.7
parms[500] = 0.42, 26.6, 0.7


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

# defines the rotation matrix R3
def rot2(theta):
    return np.array([[np.cos(theta), 0, -np.sin(theta)], [0, 1, 0], [np.sin(theta), 0, np.cos(theta)]])

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

    return normal, binormal, tangent


def evolve_tangent(triad0, psi, theta, phi):
    # normal0, binormal0, tangent0 = triad0
    # tangent = np.array([np.sin(theta) * np.sin(phi) * normal0 for theta,phi in zip(theta,phi)])
    # tangent += np.array([np.sin(theta) * np.cos(phi) * binormal0 for theta,phi in zip(theta,phi)])
    # tangent += np.array([np.cos(theta) * tangent0 for theta,phi in zip(theta,phi)])
    tangent = np.array([[np.sin(theta) * np.sin(phi), np.cos(theta), - np.sin(theta) * np.cos(phi)] for theta,phi in zip(theta,phi)])
    return tangent

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
        try:
            long = np.concatenate((getattr(elas.tail, att_str)[
                                  ::-1, :], getattr(elas.braid, att_str)[::-1, :]), axis=0)
            long = np.concatenate(
                (long, getattr(elas.loop, att_str)[::-1, :]), axis=0)
            long = np.concatenate((long, getattr(elas.loop_ref, att_str)), axis=0)
            long = np.concatenate((long, getattr(elas.braid_ref, att_str)), axis=0)
            long = np.concatenate((long, getattr(elas.tail_ref, att_str)), axis=0)
            setattr(elas, att_str, long)
        except AttributeError:
            pass


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

def invert_x(vec):
    vec[:, 0] = -vec[:, 0]
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
    try:
        invert_z(elastica.tangent)
        invert_z(elastica.normal)
        invert_z(elastica.binormal)
    except AttributeError:
        pass


def reflecty(elastica):
    invert_y(elastica.r)
    try:
        invert_z(elastica.tangent)
        invert_z(elastica.normal)
        invert_z(elastica.binormal)
    except AttributeError:
        pass
    
def reflectx(elastica):
    invert_x(elastica.r)
    invert_x(elastica.tangent)
    invert_x(elastica.normal)
    invert_x(elastica.binormal)


def plot_surface(elastica, ax, divisions=100):
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

            vec = np.array([x, y, z]) 
            vec /= np.linalg.norm(vec)
            dot = - np.dot(vec, n)

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

def electrostatic_energy(radius, L, delta, salt_conc):
    a = 1/np.tan(delta)
    m1 = 0.828
    m2 = 0.864
    c1 = 0.042
    c2 = 0.312
    lambdad, nu, bjerrum = parms[salt_conc]
    zeta = 2 * L * kT * bjerrum * nu**2 
    zed = 1 + m1 / a**2 + m2 / a**4
    ene = zeta*kn(0,2*radius/lambdad)*zed + zeta*c1*lambdad**2 / (radius**2 * a**2 * (1+c2*a**2))
    return ene


@njit(fastmath=True)
def get_writhe(r, npoints, length):
    writhe = 0
    tangent = np.empty((npoints-1,3))
    for i in range(0,npoints-1):
        tangent[i,:] = (r[i+1,:]-r[i,:])/norm((r[i+1,:]-r[i,:]))
    
    for i in range(0,npoints-1):
        ti = tangent[i, :]
        ri = r[i, :]
        for j in range(0,npoints-1):
            if i > j:
                tj = tangent[j, :]
                rj = r[j, :]
                val = multidet(ti, tj, ri-rj)
                diff = norm(ri-rj)
                # if diff < 0.001 or np.isnan(val) or np.isnan(val):
                #     val = 0
                # else:
                val /= diff**3
                writhe += val
                
    writhe /= (2*np.pi)
    writhe *= (length/npoints)**2
    
    return writhe
    
@njit(fastmath=True)
def writhe_fuller(r, r_ref, writhe_ref):

    n = min(np.shape(r)[0],np.shape(r_ref)[0])
    tan = np.empty((n-1,3))
    tan_ref = np.empty((n-1,3))
    for i in range(0,n-1):
        tan[i,:] = (r[i+1,:]-r[i,:])/norm((r[i+1,:]-r[i,:]))
        tan_ref[i,:] = (r_ref[i+1,:]-r_ref[i,:])/norm((r_ref[i+1,:]-r_ref[i,:]))
    writhe = 0
    for i in range(n-2):
        t1 = tan[i,:]
        t01 = tan_ref[i,:]
        t2 = tan[i+1,:]
        t02 = tan_ref[i+1,:]
        n0 = (t02 - t01)
        n = (t2 - t1)
        val = multidet(t01, t1, n + n0)/(1+np.dot(t1,t01))
        if np.isnan(val) or np.isinf(val) or val > 1:
            val = 0
        else:
            writhe += val
        
    writhe /= (2*np.pi)
    writhe += writhe_ref
    return writhe

class Plectoneme3D:

    def __init__(self, length, A, C, f, r0, salt_conc, n, L_loop, f_loop, p_twist, righthanded):
        self.length = length
        self.L = self.length/2
        self.A = A
        self.C = C
        self.f = f
        self.r0 = r0
        self.salt_conc = salt_conc
        self.n = n
        self.L_loop = L_loop
        self.f_loop = f_loop
        self.p_twist = p_twist
        self.righthanded = righthanded
        self.neme_bool = True
        self.lam = (self.A/self.f)**0.5
            

    def generate_plectoneme(self, npoints):
        self.npoints = npoints
        
    # Construct the loop
        # theta0 = np.pi - 15*np.pi/180
        c_loop = -1
        # c_loop = np.cos(theta0)
        try:
            a_loop, b_loop = loop_solver2(c_loop, self.L_loop, self.A, self.f_loop, self.r0, niter=100)
        except:
            a_loop = 1.1
            b_loop = 0.5
            self.neme_bool = False
            
        self.loop = Elastica3D(a_loop, b_loop, c_loop, self.L_loop, self.A, self.C, self.f_loop, overtwist=self.righthanded)
        
        if self.loop.L_complete < self.loop.L:
            self.neme_bool = False
        
        n_points_loop = round(self.npoints/2 * self.loop.L/(self.L))
        origin = np.array([0, 0, 0])
        normal0, binormal0, tangent0 = np.array([1, 0, 0]), np.array([0, 0, -1]) ,np.array([0, 1, 0])
        self.loop.triad0 = np.array([normal0, binormal0, tangent0])
        self.loop.generate_euler_angles(n_points_loop)
        self.loop.euler2cartesian(origin, phi0 = 0)
        costhetaL = self.loop.uvals[-1]
        if costhetaL < 0:
            self.neme_bool = False
        thetaL = np.arccos(costhetaL)
        
        
        
    # Construct the Braid
        tangent0, binormal0 = np.array([0, 0, -1]), np.array([0, -1, 0])
        triad0 = np.array([normal0, binormal0, tangent0])

        self.pitch = self.r0 * np.tan(thetaL)
        b_braid = (1 - costhetaL**2)**0.5 
        L_braid = np.pi*self.pitch*self.n / b_braid
        self.braid = Solenoid3D(self.A, self.C, self.n,
                                b_braid, self.r0, self.righthanded)
        self.braid.triad0 = triad0
        n_points_braid = round(self.npoints/2 * self.braid.L/self.L)
        phi0 = 0
        self.braid.generate_solenoid(n_points_braid, phi0)
        self.loop.r -= np.array([0, 0, self.loop.r[-1, 2]])
        
    # Construct the Tail

        L_tail = self.L - self.loop.L - L_braid
        
        # if the tail is too small, plectoneme is forbidden
        if L_tail < 2:
            self.neme_bool = False
            
        c_tail = costhetaL 
        self.thetac = np.arccos(c_tail)
        ang_c = (np.pi - self.thetac)/2
        m_tail = fsolve(lambda m: L_tail/(self.lam*m**0.5) - F(ang_c,m) - K(m), 0.999)
        # if diff > 0.001:
        m_tail = 1      
        # a_tail = fsolve(lambda a: L_tail - (self.A/self.f)**0.5 *
        #                 (2/(a-c_tail))**0.5*K((b_tail-c_tail)/(a-c_tail)), 1.001)[0]
        
        
        normal0, binormal0, tangent0 = np.array([1, 0, 0]), np.array([0, 0, -1]), np.array([0, 1, 0])
        self.tail = TailElastica3D(0,self.thetac, m_tail, L_tail, self.A, self.C, self.f, p_twist=self.p_twist, overtwist = False )
        self.tail.triad0 = np.array([normal0, binormal0, tangent0])
        n_points_tail = round(self.npoints/2 * abs(self.tail.L/self.L))
        if n_points_tail < 2:
            n_points_tail = 10
            self.neme_bool = False
        self.tail.generate_tail(n_points_tail)

        self.tail.r -= self.tail.r[-1,:]

        self.tail.r += np.array([self.r0, 0, self.braid.r[-1, 2]])
        self.tail.euler2cartesian(origin, 0, 0) 

    # Construct whole Neme from pieces

        min_z = self.braid.r[-1, 2]
        self.height = min_z
        self.elastica_list = [self.loop, self.braid, self.tail]
        for elastica in self.elastica_list:
            elastica.r[:, 2] -= min_z
            
        if self.n % 2 == 0:
            self.tail.r = self.tail.r[::-1,:]
        else:   
            self.tail.r = self.tail.r[::-1,:]
            reflectxy(self.tail)
        
        self.loop.r = np.delete(self.loop.r,0,axis=0)

        reflecty(self.loop)
        self.loop_ref = copy.deepcopy(self.loop)
        reflectxy(self.loop_ref)
        self.braid_ref = copy.deepcopy(self.braid)
        reflectxy(self.braid_ref)
        self.tail_ref = copy.deepcopy(self.tail)
        reflectxy(self.tail_ref)
        
        
        concatenate_elastica(self)
        self.end_to_end = 2*abs(self.r[0,1])
        
        # del_list = []
        # for i in range(self.r.shape[0]-1):
        #     dist = np.linalg.norm(self.r[i+1,:]-self.r[i,:])
        #     if dist < 0.001:
        #         print(i)
        #         print("jeb")
        #         del_list.append(i)

        self.npoints = self.r.shape[0]
        
        

    def get_twist(self):
        """
        Twist of entrire plectoneme is sum of individual twists
        Need to multiply by two for whole plectoneme as get_twist() methods
        only get the twist of half a segment
        
        """
        self.twist = self.loop.get_twist() + self.braid.get_twist() + self.tail.get_twist()
        self.twist *= 2
        return self.twist
    

    def get_writhe_fuller(self, r_ref, writhe_ref):
        writhe = writhe_fuller(self.r, r_ref, writhe_ref)
        self.writhe = writhe
        return self.writhe
    
    def get_writhe_analytic(self):
        writhe = 1
        writhe += self.n * np.sin(self.braid.theta)
        writhe *= -1
        return writhe 

    def get_link(self):
        self.writhe = get_writhe(self.r,self.npoints,self.length)
        self.link = self.get_twist() + self.writhe
        return self.link
    
    def get_link_fuller(self, r_ref, writhe_ref):
        self.writhe = writhe_fuller(self.r, r_ref, writhe_ref)
        self.link = self.get_twist() + self.writhe
        return self.link


    def plot_plectoneme(self, style= None):

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_ylabel("y")
        ax.set_xlabel("x")
        ax.set_zlabel("z")
        self.height = max(self.loop.r[0][2], self.tail.r[-1][1])
        ax.grid(False)
        ax.set_axis_off()

        ax.set_xlim(-self.height/2, self.height/2)
        ax.set_ylim(-self.height/2, self.height/2)
        ax.set_zlim(0, self.height)

        if style == "tubes":
            plot_tubes(self, ax, self.r0, divisions=50)

        elif style == "ribbon":
            plot_surface(self, ax, divisions=100)
            
        else:
            ax.plot(*self.r.T, color="red")

        plt.show()
        
    
    def get_length_discrete(self):
        self.length_disc = 0
        for i in range(0,self.npoints-1):
            self.length_disc += norm(self.r[i+1,:]-self.r[i,:])
        return self.length_disc
    
    def get_energy(self):
        self.energy = self.loop.get_energy_closed() + self.braid.get_energy() + self.tail.get_energy()
        self.energy *= 2
        # print(self.energy)
        self.energy += self.n * electrostatic_energy(self.r0, self.braid.L, self.braid.theta, self.salt_conc)
        self.energy_electrostatic = self.n * electrostatic_energy(self.r0, self.braid.L, self.braid.theta, self.salt_conc)
        # print(self.r0, self.n * electrostatic_energy(self.r0, self.braid.L, self.braid.theta, self.salt_conc))
        return self.energy
    
    def get_energy_discrete(self):
        self.energy = self.loop.get_energy_discrete() + self.braid.get_energy_discrete() + self.tail.get_energy_discrete()
        self.energy *= 2
        self.energy += self.n *  electrostatic_energy(self.r0, self.braid.L, self.braid.theta, self.salt_conc)
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
        print("Linking Number = " + str(self.get_link())+ "\n")


class Solenoid3D:
    def __init__(self, A, C, n, u0, radius, righthanded=True):
        self.A = A
        self.C = C
        self.n = n
        self.u0 = u0
        self.theta = np.arccos(self.u0) 
        self.radius = radius
        self.pitch = self.radius / np.tan(self.theta)
        self.L = np.pi*self.n*(self.pitch**2+self.radius**2)**0.5
        self.p_phi = self.A * np.pi * self.n /self.L
        self.p_psi = self.p_phi * self.u0
        self.righthanded = righthanded
        if righthanded:
            self.factor = 1
        else:
            self.factor = - 1

    def generate_euler_angles(self, npoints):
        self.npoints = npoints
        self.uvals = np.full(self.npoints, self.u0)
        self.s = np.linspace(0, self.L, self.npoints)
        self.ds = self.L/self.npoints
        self.phi = np.pi*self.n*self.s/self.L
        self.psi = (self.p_psi/self.C - np.pi*self.n*self.u0/self.L)*self.s

    def euler2cartesian(self, origin, phi0, psi0):
        self.phi += phi0
        self.psi += psi0
        self.normal, self.binormal, self.tangent = evolve_frame(
            self.triad0, self.psi, self.theta, self.phi)
        self.r = np.cumsum(self.tangent, axis=0)*self.ds
        if self.righthanded:
            reflecty(self)
        self.r += origin
        
    def generate_solenoid(self, npoints, phi0):
        """
        Generates positions, tangents and normals for a right/left handed solenoid of self.n turns
        Solenoids are aligned in the z-axis and 

        """
        self.npoints = npoints
        self.ds = self.L/self.npoints
        self.s = np.linspace(0, self.L, self.npoints)
        self.phi = self.factor*np.pi*self.n*self.s/self.L
        self.phi -= phi0
        self.psidot = (self.p_psi/self.C - self.factor*np.pi*self.n*self.u0/self.L)
        self.r = np.array([[ self.radius * np.cos(phi),  self.radius * np.sin(phi),  - self.factor * self.pitch * phi] for phi in self.phi])
        self.r = np.delete(self.r, -1, axis = 0)
        self.r = np.delete(self.r, 0, axis = 0)
        self.npoints -= 2
        self.tangent = np.array([[ - self.factor * self.radius * np.sin(phi),  self.factor * self.radius * np.cos(phi), - self.pitch] for phi in self.phi]/(self.radius**2+self.pitch**2)**0.5)
        
    def plot_solenoid(self, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(*self.r.T, color='red')
        ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
        set_axes_equal(ax)
        plot_surface(ax, divisions=divisions)
        plt.show()

    def get_twist(self):
        self.twist =  self.factor * self.n * self.u0 * self.A / (2 * self.C)
        return self.twist
        
    def get_energy(self):
        self.bend = self.A * self.L * (self.radius/(self.radius**2+self.pitch**2))**2 / 2
        self.twist_ene = 2*np.pi**2*self.C*(self.get_twist())**2/(self.L)
        self.energy = self.bend + self.twist_ene
        return self.energy
    
    def get_energy_discrete(self):
        bend = 0
        tangent = np.empty((self.npoints-1,3))
        for i in range(0,self.npoints-1):
            tangent[i,:] = (self.r[i+1,:]-self.r[i,:])/norm((self.r[i+1,:]-self.r[i,:]))
        dt = np.empty((self.npoints-2,3))
        for i in range(0,self.npoints-2):
            dt[i,:] = (tangent[i+1,:]-tangent[i,:])
        for i in range(0,self.npoints-2):
            bend += np.dot(dt[i,:],dt[i,:])
        
        self.length_disc = 0
        for i in range(0,self.npoints-2):
            self.length_disc += norm(self.r[i+1,:]-self.r[i,:])
        bend *= self.A / (2*self.ds)
        self.bend_disc = bend
        twist = self.p_psi**2 * (1/self.C) * self.L / 2
        self.twist_ene_disc = twist
        energy = twist + bend
        return energy

    def printvals(self):
        print("Length = " + str(self.L))
        if self.righthanded == True:
            print("Braid is a right handed solenoid")
        else:
            print("Braid is a left handed solenoid")
        print("Number of braids = " + str(self.n))
        print("Radius = " + str(self.radius))
        print("Pitch = " + str(self.pitch))
        print("Twist = " + str(self.get_twist())+ "\n")


class Elastica3D:

    def __init__(self, a, b, c, L, A, C, f, increase = True, overtwist = False ):
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
        self.L_complete = self.lam * (2/(self.a-self.c))**0.5 * K(self.m)

        self.p_psi = (p1 + p2)/2
        self.p_phi = (p1 - p2)/2
        
        if not overtwist:
            self.p_psi *= -1
            self.p_phi *= -1

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

        # self.normal, self.binormal, self.tangent = evolve_frame(
        #     self.triad0, self.psi, self.theta, self.phi)
        self.tangent = evolve_tangent(self.triad0, self.psi, self.theta, self.phi)
        self.r = np.cumsum(self.tangent, axis=0)*self.ds
        self.r += origin
        self.r = np.insert(self.r, 0, origin, axis=0)

    def plot_elastica(self, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(*self.r.T, color='red')

        ax.set_ylabel("y")
        ax.set_xlabel("x")
        ax.set_zlabel("z")
        set_axes_equal(ax)
        plot_surface(self, ax, divisions)
        plt.show()

    def get_twist(self):
        self.twist = self.L*self.p_psi/(2*np.pi*self.C)
        return self.twist
    
    def get_energy(self):
        self.energy = self.f * (self.a + self.b + self.c) * self.L
        self.energy += 2*(2*self.A*self.f*(self.a - self.c))**0.5*E(np.pi/2,self.m) - 2*self.a*self.f*self.L
        self.energy += self.p_psi**2 * (1/self.C - 1/self.A) * self.L / 2
        dell = self.f * (self.a + self.b + self.c) * self.L - self.p_psi**2 * (1/self.A) * self.L / 2
        self.twist_ene = self.p_psi**2 * (1/self.C) * self.L / 2
        self.bend = (self.energy - self.twist_ene - dell)/2
        self.pull = self.bend - dell
        return self.energy
    
    def get_energy_closed(self):
        self.pull = 0
        self.bend =  self.f * (self.a + self.b + self.c) * self.L - self.p_psi**2 * ( 1 / self.A) * self.L / 2
        self.twist = self.p_psi**2 * ( 1 / self.C) * self.L / 2
        self.energy = self.pull + self.bend + self.twist
        return self.energy
    
    def get_energy_discrete(self):
        pull = -self.f * abs(self.r[-1,1])
        bend = 0
        tangent = np.empty((self.npoints-1,3))
        for i in range(0,self.npoints-1):
            tangent[i,:] = (self.r[i+1,:]-self.r[i,:])/norm((self.r[i+1,:]-self.r[i,:]))
        dt = np.empty((self.npoints-2,3))
        for i in range(0,self.npoints-2):
            dt[i,:] = (tangent[i+1,:]-tangent[i,:])
        for i in range(0,self.npoints-2):
            bend += np.dot(dt[i,:],dt[i,:])
        bend *= self.A / (2*self.ds)
        twist = self.p_psi**2 * (1/self.C) * self.L / 2
        self.bend_disc = bend
        self.twist_disc = twist
        self.pull_disc = pull
        energy = twist + bend + pull
        return energy
        

        

    def printvals(self):
        print("L = " + str(self.L))
        print("a = " + str(self.a))
        print("b = " + str(self.b))
        print("c = " + str(self.c))
        print("m = " + str(self.m))
        print("Twist = " + str(self.get_twist()) + "\n")
        
        
class TailElastica3D:

    def __init__(self, theta0, thetaL, m, L, A, C, f, p_twist, overtwist = False ):
        self.theta0 = theta0
        self.thetaL = thetaL
        self.m = m
        self.A = A
        self.C = C
        self.f = f
        self.lam = (self.A/self.f)**0.5
        self.L = L
        self.p_twist = p_twist

    def generate_tail(self, npoints):

        self.npoints = npoints
        self.ds = self.L/self.npoints
        k = self.m**0.5
        self.s = np.linspace(0,self.L, self.npoints)
        
        ang0 = np.full(self.npoints,(self.theta0-np.pi)/2)
        
        if np.isclose(self.m,1):
            # self.theta = self.theta0 + (self.thetaL-self.theta0)*self.s/self.L
            self.theta = np.array([4*np.arctan( np.tan(self.thetaL/4) * np.exp((s-self.L)/self.lam)) for s in self.s])
        else:
            self.theta = np.pi + 2*am(self.s/(self.lam*k)+ F(ang0,self.m), self.m)
        # discrete version
        x = np.cumsum(np.cos(self.theta)) * self.ds 
        y = np.cumsum(np.sin(self.theta)) * self.ds 
        x[0] = 0
        y[0] = 0
        self.r = np.array([[0,x,y] for x,y in zip(x,y)])
        # self.r = np.array([[0, xfunc(s,self.theta0,self.thetaL,self.m,self.A,self.f),yfunc(s,self.theta0,self.thetaL,self.m,self.A,self.f)] for s in self.s])
        self.psi = self.p_twist * self.s / self. C
        self.phi = np.full(npoints,0)
        

    def euler2cartesian(self, origin, phi0=0, psi0=0):
        self.phi += phi0
        self.psi += psi0

        # self.normal, self.binormal, self.tangent = evolve_frame(
        #     self.triad0, self.psi, self.theta, self.phi)
        self.tangent = evolve_tangent(self.triad0, self.psi, self.theta, self.phi)
        # self.r = np.cumsum(self.tangent, axis=0)*self.ds

    def plot_elastica(self, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(*self.r.T, color='red')

        ax.set_ylabel("y")
        ax.set_xlabel("x")
        ax.set_zlabel("z")
        set_axes_equal(ax)
        plot_surface(self, ax, divisions)
        plt.show()

    def get_twist(self):
        self.twist = self.L*self.p_twist/(2*np.pi*self.C)
        return self.twist
    
    # def get_energy(self):
    #     self.energy = 
    #     return self.energy

    
    def get_energy_discrete(self):
        pull = -self.f * abs(self.r[-1,1])
        bend = 0
        tangent = np.empty((self.npoints-1,3))
        for i in range(0,self.npoints-1):
            tangent[i,:] = (self.r[i+1,:]-self.r[i,:])/norm((self.r[i+1,:]-self.r[i,:]))
        dt = np.empty((self.npoints-2,3))
        for i in range(0,self.npoints-2):
            dt[i,:] = (tangent[i+1,:]-tangent[i,:])
        for i in range(0,self.npoints-2):
            bend += np.dot(dt[i,:],dt[i,:])
        bend *= self.A / (2*self.ds)
        twist = self.p_twist**2 * (1/self.C) * self.L / 2
        self.bend_disc = bend
        self.twist_disc = twist
        self.pull_disc = pull
        energy = twist + bend + pull
        return energy
        

    def printvals(self):
        print("L = " + str(self.L))
        print("theta_c = " + str(self.thetac))
        print("m = " + str(self.m))
        print("Twist = " + str(self.get_twist()) + "\n")


class Twisted_Line3D:

    def __init__(self, length, A, C, f, link):
        self.length = length
        self.A = A
        self.C = C
        self.f = f
        self.link = link
        self.L  = length/2

    def generate_twisted_line(self, npoints):
        self.npoints = npoints
        self.s = np.linspace(0, self.length, npoints)
        self.psi = 2 * self.s * np.pi * self.link / self.length
        self.r = np.array([[ 0, s, 0 ] for s in self.s ])
        self.tangent = np.array([[0,1,0] for s in self.s])
        tangent0 = np.array([0,1,0])
        normal0 = np.array([1,0,0])
        binormal0 = np.array([0,0,1])
        self.normal = np.array([rot(psi, tangent0, normal0) for psi in self.psi])
        self.binormal = np.array([rot(psi, tangent0, binormal0) for psi in self.psi])
        
    def get_energy(self):
        ene =  - self.f * self.length
        ene += 2 * np.pi**2 * self.C * self.link **2 / self. length
        self.energy = ene
        return self.energy
    
    def plot_twisted_line(self, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlim(-self.length/2, self.length/2)
        ax.set_ylim(0, self.length)
        ax.set_zlim(-self.length/2, self.length/2)
        ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
        set_axes_equal(ax)
        ax.grid(False)
        ax.set_axis_off()
        plot_surface(self, ax, divisions)
        plt.show()
        
        
class Curl3D:

    def __init__(self, length, A, C, f, link):
        self.length = length
        self.A = A
        self.C = C
        self.f = f
        self.link = link
        self.L  = length/2
        if self.link <= 0:
            self.twist = - self.link + 1
        else:
            self.twist = self.link - 1
        if self.twist < 0:
            self.overtwist = False
        else:
            self.overtwist = True

    def generate_curl(self, npoints):
        c = -1
        b = 1
        a = 1
        
        m = (b-c)/(a-c)
        self.curl = Elastica3D(a, b, c, self.L, self.A, self.C, self.f, overtwist=self.overtwist)
        self.npoints = npoints
        n_points_loop = round(self.npoints/2)
        origin = np.array([0, 0, 0])
        normal0, binormal0, tangent0 = np.array([1, 0, 0]), np.array([0, 0, -1]) ,np.array([0, 1, 0])
        self.curl.triad0 = np.array([normal0, binormal0, tangent0])
        self.s = np.linspace(0, self.L, n_points_loop)
        self.curl.generate_euler_angles(n_points_loop)
        self.curl.psi = 2*np.pi*self.s*self.twist/self.length
        self.curl.euler2cartesian(origin, phi0 = 0)
        
        
    def get_energy(self):
        ene =  - self.f * self.length + 8 * (self.A * self.f)**0.5
        ene += 2 * np.pi**2 * self.C * self.twist **2 / self. length
        self.energy = ene
        return self.energy
    
    def plot_curl(self, divisions):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlim(-self.length/4, self.length/4)
        ax.set_ylim(-self.length/2, self.length/2)
        ax.set_zlim(-self.length/4, self.length/4)
        ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
        set_axes_equal(ax)
        ax.grid(False)
        ax.set_axis_off()
        self.curl_2 = copy.deepcopy(self.curl)
        reflectxy(self.curl_2)
        for elas in [self.curl, self.curl_2]:
            elas.r -= self.curl.r[0]
        plot_surface(self.curl,ax,divisions)
        plot_surface(self.curl_2,ax,divisions)
        plt.show()
        

def costheta(a, b, c, A, f, npoints, increase=True):
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
        if not increase:
            s = s[::-1]
        costheta = c + (b-c)*sn(s*((a-c)/2)**0.5/lam, m)**2

    costheta[0] = c

    return costheta

# @njit(fastmath=True)
def costheta_incomp(a, b, c, L, A, f, npoints, increase=True):

    m = (b-c)/(a-c)
    lam = (A/f)**0.5
    # if np.isclose(b, c):
    #     costheta = np.full(npoints, b)
    # else:
    s = np.linspace(0, L, npoints)
    if not increase:
        s = s[::-1]
    costheta = c + (b-c)*sn(s*((a-c)/2)**0.5/lam, m)**2

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

# @njit(fastmath=True)
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


# def xfunc(a, b, c, A, f, npoints):
#     lam = (A/f)**0.5
#     m = (b-c)/(a-c)
#     L = lam*(2/(a-c))**0.5*K(m)
#     ds = L/npoints
#     uvals = costheta(a, b, c, A, f, npoints)
#     sintheta = (1-uvals**2)**0.5
#     phi = phi_func(a, b, c, A, f, npoints)
#     val = np.cumsum(sintheta*np.sin(phi))*ds
#     return val


# def yfunc(a, b, c, A, f, npoints):
#     lam = (A/f)**0.5
#     m = (b-c)/(a-c)
#     k = m**0.5
#     L = lam*k*(2/(b-c))**0.5*K(m)
#     ds = L/npoints
#     uvals = costheta(a, b, c, A, f, npoints)
#     sintheta = (1-uvals**2)**0.5
#     phi = phi_func(a, b, c, A, f, npoints)
#     val = np.cumsum(sintheta*np.cos(phi))*ds
#     return val

def xfunc(s,theta0,thetaL,m,A,f):
    lam = (A/f)**0.5
    k = m**0.5
    theta = np.pi + 2*am(s/(lam*k) + F((theta0-np.pi)/2,m),m)
    xval = (2/m-1)*s - 2*lam/k*( E((theta-np.pi)/2,m)- E((theta0-np.pi)/2,m))
    return xval

def yfunc(s,theta0,thetaL,m,A,f):
    lam = (A/f)**0.5
    k = m**0.5
    u = s/(lam*k) + F((theta0-np.pi)/2,m)
    u0 = F((theta0-np.pi)/2,m)
    yval = 2 * lam/k * ( dn(u, m) - dn(u0, m))
    return yval



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

# @njit(fastmath=True)
def xend_incomp(a, b, c, L, A, f, npoints):
    ds = L/npoints
    uvals = costheta_incomp(a, b, c, L, A, f, npoints)
    sintheta = (1 - uvals **2) ** 0.5
    phi = phi_func_incomp(a, b, c, L, A, f, npoints)
    xend = np.sum( sintheta * np.sin(phi) ) * ds
    return xend

def zend_incomp(a, c, m, L, A, f, npoints):
    lam = (A/f)**0.5
    ds = L/npoints
    s = np.linspace(0, L, npoints)
    uvals = c + m*(a-c)*sn(s*((a-c)/2)**0.5/lam, m)**2
    zend = np.sum(uvals)*ds
    return zend


def tail_solver(b, c, A, f, npoints, r0):
    a = newton(lambda a: xend(a, b, c, A, f, npoints)-r0, 1.2)
    return a

# @njit(fastmath=True)
def funcy(a,c,L,lam,m):
    val = a*L - lam*(2*(a-c))**0.5*E(am(L*((a-c)/2)**0.5/lam, m), m)
    return val

def loop_func2(a, c, L, A, f, npoints):
    lam = (A/f)**0.5
    # m = fsolve(lambda m: a*L - lam*(2*(a-c))**0.5*E(am(L*((a-c)/2)**0.5/lam, m), m), 0.8)[0]
    m = fsolve(lambda m: funcy(a, c, L, lam, m),0.8)[0]
    b = c + ( a - c ) * m
    return xend_incomp(a, b, c, L, A, f, npoints) 



def loop_solver2(c, L, A, f, r0, niter):
    lam = (A/f)**0.5
    a = newton(lambda a: loop_func2(a, c, L, A, f, niter) + r0, 1.0001)
    # a = fsolve(lambda a: loop_func2(a, c, L, A, f, niter) + r0, 1.0001)[0]
    # m = fsolve(lambda m: a*L - lam*(2*(a-c))**0.5*E(am(L*((a-c)/2)**0.5/lam, m), m), 0.8)[0]
    m = fsolve(lambda m: funcy(a, c, L, lam, m),0.8)[0]
    b = c + (a-c)*m
    return a, b

# nemey_daniel = Plectoneme3D(length=500, A=50*kT, C=100*kT, f=1, r0 = 3, salt_conc = 100, n = 2, L_loop = 40, f_loop=1, p_twist = 0, righthanded = False)
# nemey_daniel.generate_plectoneme(1000)
# nemey_daniel.plot_plectoneme(style="ribbons")
# nemey_daniel.get_link()
# twist = nemey_daniel.twist
# writhe = nemey_daniel.writhe


nemey_daniel = Plectoneme3D(length=500, A=50*kT, C=100*kT, f=1, r0 = 4, salt_conc = 100, n = 3, L_loop = 40, f_loop=1, p_twist = 0, righthanded = False)
nemey_daniel.generate_plectoneme(1000)
nemey_daniel.plot_plectoneme(style="ribbons")

# xyzfile = "neme.xyz"
# with open(xyzfile,"w") as file:
#     n = nemey_daniel.r.shape[0]
#     file.write(str(n))
#     file.write("\n")
#     file.write("\n")
#     for i in range(nemey_daniel.r.shape[0]):
#         file.write("C " + str(nemey_daniel.r[i,0]) + " " + str(nemey_daniel.r[i,1]) + " " + str(nemey_daniel.r[i,2]) +"\n")
        
    
# xyzfile = "line.xyz"
# r_line = np.zeros((1000,3))
# for i in range(1000):
#     r_line[i,0] = 8*np.sin(i/75)
#     r_line[i,1] = 4*np.cos(i/100)
#     r_line[i,2] = i*300/1000

# with open("line.xyz","w") as file:
#     n = r_line.shape[0]
#     file.write(str(n))
#     file.write("\n")
#     file.write("\n")
#     for i in range(r_line.shape[0]):
#         file.write("C " + str(r_line[i,0]) + " " + str(r_line[i,1]) + " " + str(r_line[i,2]) +"\n")
