#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 11:36:12 2022

analysis of T_motif simulations

@author: michaelselby
"""
import numpy as np
import pdb_miek
import tools
import tools_fast_math
import matplotlib.pyplot as plt
from molecular_contour import get_molecular_contour

direc = "/Users/michaelselby/OneDrive - Nexus365/DPhil/Behnam/T_motifs_michael/"
polarity = 3
bulge = 3
direc += str(polarity) + "E/" + str(bulge) + "/"
crd = direc + str(polarity)+"E_"+str(bulge)+"_longrun_nw.crd"
top = direc + str(polarity)+"E_"+str(bulge)+".top"

def t_motif_angles(conf):
    from skspatial.objects import Line

    strand1, strand2, strand3, strand4 = conf.strand_list
    strand1.reverse_strand()
    strand3.reverse_strand()
    
    strandApos = np.array(strand4.get_atom_list("C1'"))
    strandBpos = np.array(strand3.get_atom_list("C1'")[:strand4.nres])
    buff = 3
    bp = strand4.nres
    rb = get_molecular_contour(bp, strandApos, strandBpos, linear=True)
    rb = rb[buff:-buff,:]
    lineb = Line.best_fit(rb)
    gradb = lineb.direction
    gradb /= np.linalg.norm(gradb)
    
    length = 15
    strandApos = np.array(strand1.get_atom_list("C1'")[-length:])
    strandBpos = np.array(strand2.get_atom_list("C1'")[strand2.nres-length:])

    ra = get_molecular_contour(length, strandApos, strandBpos, linear=True)
    ra = ra[buff:-buff,:]
    linea = Line.best_fit(ra)
    grada = linea.direction
    grada /= np.linalg.norm(grada)
    
    strandApos = np.array(strand1.get_atom_list("C1'")[:length-strand1.nres])
    strandBpos = np.array(strand2.get_atom_list("C1'")[:length])

    rc = get_molecular_contour(length, strandApos, strandBpos, linear=True)
    rc = rc[buff:-buff,:]
    linec = Line.best_fit(rc)
    gradc = linec.direction
    gradc /= np.linalg.norm(gradc)

    # make sure all vectors are pointing in same direction by checking againt centre to end-point vector
    # mid = np.array(strand1.get_atom_list("C1'")[int(strand1.nres/2)])
    mid = np.array(strand4.get_atom_list("C1'")[-1])
    a_end = strand2.get_atom_list("C1'")[0]
    c_end = strand2.get_atom_list("C1'")[-1]
    b_end = strand4.get_atom_list("C1'")[-1]
    
    for end, grad in zip([a_end,b_end,c_end], [grada,gradb,gradc]):
        if np.dot(end - mid, grad) < 0:
            grad *= -1


    
    cosalpha = np.dot(grada,gradc)
    alpha = np.arccos(cosalpha)*180/np.pi
    
    cosbeta = np.dot(gradb,gradc)
    beta = np.arccos(cosbeta)*180/np.pi
    
    cosphi = np.dot(-grada,gradb)
    phi = np.arccos(cosphi)*180/np.pi 
    
    # use sign convention outlined in Behnam paper
    sign = np.dot( np.cross(grada,gradc), np.cross(gradb,gradc))
    if sign < 0:
        alpha = abs(alpha)
        beta = abs(beta)
    else:
        print("boon")
        alpha = - abs(alpha)
        beta =  - abs(beta)
        
    sign = np.dot( np.cross(np.cross(grada,gradc), np.cross(gradb,gradc)), gradc)
    if sign < 0:
        phi = abs(phi)
    else:
        phi = - abs(phi)
        
    return alpha, beta, phi

traj = pdb_miek.trajectory(crd,circular=False)
res = traj.load_topology(top)

traj.process_configurations([("angles",t_motif_angles)], max_steps=2000)
# traj.process_configurations([("angles",t_motif_angles)])
dic = traj.quant_dict
conf_test = traj.last_conf

alpha, beta, phi = t_motif_angles(conf_test)

alpha_list, beta_list, phi_list = [], [], []

for alpha, beta, phi in traj.quant_dict["angles"]:
    alpha_list.append(alpha)
    beta_list.append(beta)
    phi_list.append(phi)

file_name = direc + "angles.txt"
with open(file_name, "w") as file:
    for alpha, beta, phi in zip(alpha_list,beta_list,phi_list):
        file.write("{:3.8f} {:3.8f} {:3.8f} \n".format(alpha, beta, phi))

av_alpha, std_alpha = np.mean(alpha_list), np.std(alpha_list)
av_beta, std_beta = np.mean(beta_list), np.std(beta_list)
av_phi, std_phi = np.mean(phi_list), np.std(phi_list)

# strand1, strand2, strand3, strand4 = conf_test.strand_list

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlim(0,100)
# ax.set_ylim(0,100)
# ax.set_zlim(0,100)
# colour_i = 0
# colour_list = ["blue", "red", "green", "purple"]
# for strand in conf_test.strand_list:
#     for base in strand.base_list:
#         ax.scatter(base.bb[0],base.bb[1],base.bb[2], color=colour_list[colour_i], marker="x")
#     ax.scatter(strand.base_list[0].bb[0],strand.base_list[0].bb[1],strand.base_list[0].bb[2], color="black")
#     ax.scatter(strand.base_list[-1].bb[0],strand.base_list[-1].bb[1],strand.base_list[-1].bb[2], color="yellow")
#     colour_i += 1
    
# mid = np.array(conf_test.strand_list[0].get_atom_list("C1'")[int(strand1.nres/2)])
# ax.scatter(mid[0], mid[1], mid[2], color="orange")
    
# for end in [a_end, b_end, c_end]:
#     ax.scatter(end[0],end[1],end[2], color="black")

# for r in [ra,rb,rc]:
#     for vec in r:
#         ax.scatter(vec[0],vec[1],vec[2], color="orange", marker='x')
        
# for r in [r1,r2]:
#     for vec in r:
#         ax.scatter(vec[0],vec[1],vec[2], color="black")
        
# for line in [linea,lineb,linec]:
#     points = [np.array(line.point + line.direction * s) for s in np.linspace(0,50,100)]
#     for point in points:
#         ax.scatter(point[0],point[1],point[2], color="black")
              
# plt.show()