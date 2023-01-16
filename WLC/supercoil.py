#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 10:43:11 2023

Calculates the energies and configurations of supercoiled structures
and determines the elastic minima

@author: michaelselby
"""
import numpy as np
import time
import matplotlib.pyplot as plt
from three_dimensional_elastica import Plectoneme3D, Twisted_Line3D
import copy

kT = 4.114  # kT in pN nm at 298K


# nemey_daniel = Plectoneme3D(300, A=50*kT, C=100*kT, f=1.5, r0 = 2.5, n = 2, L_loop = 50, p_twist = 0, righthanded = False)
# nemey_daniel.generate_plectoneme(2000)
# nemey_daniel.plot_plectoneme(style="ribbon")
# nemey_daniel.printvals()

# L = 300 # length in nm
nbp = 600
L = 0.34 * nbp
Lk0 = nbp / 10.55

sigma_vals = np.linspace( -0.03 , -0.08, 20)
f_vals = np.linspace( 0.5, 2, 20)

init_dict = {}
init_dict["lowest energy"] = 10**10

C = 100 * kT

dict_of_dicts = {}
for sigma in sigma_vals:
    deltaLk = sigma * Lk0
    # deltaLk += 15*np.pi/180
    for force in f_vals:
        dic = copy.deepcopy(init_dict)
        dic["Link"] = deltaLk
        dic["Nemey"] = False
        ene_line = - force * L + 2 * np.pi**2 * C * deltaLk ** 2 / L
        dic["Twisted Line Energy"] = ene_line
        dict_of_dicts[(sigma,force)] = dic
        

L_vals = np.linspace(10,100,50)
r_vals = np.linspace(2,10,50)
n_vals = [1,2,3,4]

# L_vals = [70]
# f_vals = [1]
# n_vals = [3]
npoints = 500

n_nemes = len(r_vals) * len(L_vals) * len(n_vals) * len(f_vals)


salt_conc = 100

points = []
list2 = []
st = time.time()

for f in f_vals:
    # f_loop_vals = np.linspace(0.5*f,f,2)
    f_loop_vals = [f]
    for n in n_vals:
        fuller = False
        for r in r_vals:
            min_dot = 1
            neme_exists = False
            for L_loop in L_vals:
                for f_loop in f_loop_vals:
                    neme0 = Plectoneme3D(length=L, A=50*kT, C=100*kT, f=f, r0 = r, salt_conc = salt_conc , n = n, L_loop = L_loop, f_loop = f_loop, p_twist = 0, righthanded = False)
                    neme0.generate_plectoneme(npoints)
                    dot = np.dot(neme0.loop.tangent[-1],neme0.braid.tangent[0])
                    if dot < min_dot and neme0.neme_bool:
                        min_dot = dot
                        neme = neme0
                        neme_exists = True
                    
            if neme_exists:
                
                if not fuller:
                    deltaLk_zero = neme.get_link()
                    fuller = True
                    # r_ref = neme.r
                    # writhe_ref = neme.writhe
                else:
                    fuller = True
                    deltaLk_zero = neme.get_link_fuller(r_ref, writhe_ref)
                    # r_ref = neme.r
                    # writhe_ref = neme.writhe
                    
                r_ref = neme.r
                writhe_ref = neme.writhe
                energy = neme.get_energy_discrete()
                
                delta = neme.braid.theta  
                
                for sigma in sigma_vals:
                    dic = dict_of_dicts[(sigma,f)]
                    deltaLk = sigma * Lk0
                    twist_diff = deltaLk - deltaLk_zero
                    p_twist = neme.C * twist_diff * np.pi / neme.tail.L
                    neme_ene = energy + neme.tail.L * p_twist**2/(neme.C)
                    nemey = False

                    if neme_ene < dic["lowest energy"]:
                        dic["lowest energy"] = neme_ene
                        if energy < dic["Twisted Line Energy"]:
                            nemey = True
                        dic["n"] = n
                        dic["Nemey"] = nemey
                        dic["Radius"] = round(r, 4)
                        dic["z/L"] = round(neme.end_to_end / L, 4)
                    
et = time.time()
elapsed_time = et - st
nemes_per_second =  n_nemes / elapsed_time
print("Generated ", n_nemes, " nemes at a rate of ", nemes_per_second, " nemes per second" )
   

file_name = "neme_data.dat"
template = "{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n"
fields = ["Sigma", "Force", "Nemey", "nbraids", "Radius", "z/L" ]
with open(file_name, "w") as file:
    file.write(template.format(*fields))
    for key, dic in dict_of_dicts.items():
        file.write(template.format(round(key[0], 4), round(key[1],4), dic["Nemey"], dic["n"], dic["Radius"], dic["z/L"]))




