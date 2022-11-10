#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 11:36:12 2022

analysis of T_motif simulations

@author: michaelselby
"""
import sys
import numpy as np
import pdb_miek
import tools
import tools_fast_math
import matplotlib.pyplot as plt
from molecular_contour import get_molecular_contour

direc = sys.argv[1]
polarity = sys.argv[2]
bulge = sys.argv[3]


direc += str(polarity) + "E/" + str(bulge) + "/"
crd = direc + str(polarity)+"E_"+str(bulge)+"_longrun_nw.crd"
top = direc + str(polarity)+"E_"+str(bulge)+".top"
# crd = direc + "test.crd"

def t_motif_angles(conf):
    from skspatial.objects import Line
    

    strand1, strand2, strand3, strand4 = conf.strand_list
    mid = np.array(strand1.get_atom_list("C1'")[int(strand1.nres/2)])
    
    end_base_1 = strand1.base_list[0]
    end_base_2 = strand1.base_list[-1]


    # put strands 1 and 3 in 3' -> 5' direction
    # put strands 2 and 4 in 5' -> 3' direction
    if end_base_1.resname[2] == "3":
        strand2.reverse_strand()
        strand4.reverse_strand()
    else:
        strand1.reverse_strand()
        strand3.reverse_strand()
        
    
    a_end = end_base_1.get_atom("C1'")
    c_end = end_base_2.get_atom("C1'")
    
    # a_end = end_base_2.get_atom("C1'")
    # c_end = end_base_1.get_atom("C1'")
    
    strandApos = np.array(strand4.get_atom_list("C1'"))
      
    b_end1 = strand4.get_atom_list("C1'")[0]
    b_end2 = strand4.get_atom_list("C1'")[-1]
    dist1 = np.linalg.norm(mid-b_end1)
    dist2 = np.linalg.norm(mid-b_end2)
    
    # if dist1 > dist 2, that means strand 4s 5' end is further away from T-motif centre
    # hence we want strandBpos first strand4.nres base-pairs 
    
    if dist1 > dist2:
        b_end = b_end1
        strandBpos = np.array(strand3.get_atom_list("C1'")[-strand4.nres:])
    else:
        b_end = b_end2
        strandBpos = np.array(strand3.get_atom_list("C1'")[:strand4.nres])

            
    buff = 3
    bp = strand4.nres
    rb = get_molecular_contour(bp, strandApos, strandBpos, linear=True)
    rb = rb[buff:-buff,:]
    lineb = Line.best_fit(rb)
    gradb = lineb.direction
    gradb /= np.linalg.norm(gradb)
    
    length = 15
    strandApos = np.array(strand1.get_atom_list("C1'")[:length])
    strandBpos = np.array(strand2.get_atom_list("C1'")[:length])
    
    
    r1 = get_molecular_contour(length, strandApos, strandBpos, linear=True)
    r1 = r1[buff:-buff,:]
    line1 = Line.best_fit(r1)
    grad1 = line1.direction
    grad1 /= np.linalg.norm(grad1)
    
    strandApos = np.array(strand1.get_atom_list("C1'")[strand1.nres-length:])
    strandBpos = np.array(strand2.get_atom_list("C1'")[strand2.nres-length:])



    r2 = get_molecular_contour(length, strandApos, strandBpos, linear=True)
    r2 = r2[buff:-buff,:]
    line2 = Line.best_fit(r2)
    grad2 = line2.direction
    grad2 /= np.linalg.norm(grad2)
    

    # ra, grada = r1, grad1
    # rc, gradc = r2, grad2
    rc, gradc = r1, grad1
    ra, grada = r2, grad2


    # make sure all vectors are pointing in same direction by checking againt centre to end-point vector 
    
    for end, grad in zip([a_end,b_end,c_end], [grada,gradb,gradc]):
        # print(np.dot((end - mid)/(np.linalg.norm(end - mid)), grad))
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
        alpha = - abs(alpha)
        beta =  - abs(beta)
        
    sign = np.dot( np.cross(np.cross(grada,gradc), np.cross(gradb,gradc)), gradc)
    if sign < 0:
        phi = abs(phi)
    else:
        phi = - abs(phi)
        
    # print(alpha,beta,phi)
        
    return alpha, beta, phi

traj = pdb_miek.trajectory(crd,circular=False)
traj.polarity = polarity
res = traj.load_topology(top)

# traj.process_configurations([("angles",lambda x: print("jubz"))], max_steps=3)
traj.process_configurations([("angles",t_motif_angles)], max_steps=1000)
# traj.process_configurations([("angles",t_motif_angles)])
dic = traj.quant_dict

alpha_list, beta_list, phi_list = [], [], []

for alpha, beta, phi in traj.quant_dict["angles"]:
    alpha_list.append(alpha)
    beta_list.append(beta)
    phi_list.append(phi)

file_name = direc + "angles_test.txt"
with open(file_name, "w") as file:
    for alpha, beta, phi in zip(alpha_list,beta_list,phi_list):
        file.write("{:3.8f} {:3.8f} {:3.8f} \n".format(alpha, beta, phi))

av_alpha, std_alpha = np.mean(alpha_list), np.std(alpha_list)
av_beta, std_beta = np.mean(beta_list), np.std(beta_list)
av_phi, std_phi = np.mean(phi_list), np.std(phi_list)

# chun_alpha = np.mean([abs(a) for a in alpha_list])





def T_motif_angles_debug(conf):
    
    
    from skspatial.objects import Line
    
    strand1,strand2,strand3,strand4 = conf.strand_list
    
    mid = np.array(strand1.get_atom_list("C1'")[int(strand1.nres/2)])
    

    end_base_1 = strand1.base_list[0]
    end_base_2 = strand1.base_list[-1]
    # put strands 1 and 3 in 3' -> 5' direction
    # put strands 2 and 4 in 5' -> 3' direction
    if end_base_1.resname[2] == "3":
        strand2.reverse_strand()
        strand4.reverse_strand()
    else:
        strand1.reverse_strand()
        strand3.reverse_strand()
        
        
    
    # a_end = end_base_2.get_atom("C1'")
    # c_end = end_base_1.get_atom("C1'")
    
    a_end = end_base_1.get_atom("C1'")
    c_end = end_base_2.get_atom("C1'")
    print(end_base_2.resname)
    
    strandApos = np.array(strand4.get_atom_list("C1'"))
      
    b_end1 = strand4.get_atom_list("C1'")[0]
    b_end2 = strand4.get_atom_list("C1'")[-1]
    dist1 = np.linalg.norm(mid-b_end1)
    dist2 = np.linalg.norm(mid-b_end2)
    
    # if dist1 > dist 2, that means strand 4s 5' end is further away from T-motif centre
    # hence we want strandBpos first strand4.nres base-pairs 
    if dist1 > dist2:
        b_end = b_end1
        strandBpos = np.array(strand3.get_atom_list("C1'")[-strand4.nres:])
    else:
        b_end = b_end2
        strandBpos = np.array(strand3.get_atom_list("C1'")[:strand4.nres])

    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0,100)
    ax.set_ylim(0,100)
    ax.set_zlim(0,100)
     
    for r in [strandApos,strandBpos]:
        for vec in r:
            ax.scatter(vec[0],vec[1],vec[2], color="blue", marker='x')

            
    buff = 3
    bp = strand4.nres
    rb = get_molecular_contour(bp, strandApos, strandBpos, linear=True)
    rb = rb[buff:-buff,:]
    lineb = Line.best_fit(rb)
    gradb = lineb.direction
    gradb /= np.linalg.norm(gradb)
    
    length = 15
    strandApos = np.array(strand1.get_atom_list("C1'")[:length])
    strandBpos = np.array(strand2.get_atom_list("C1'")[:length])
    for r in [strandApos,strandBpos]:
        for vec in r:
            ax.scatter(vec[0],vec[1],vec[2], color="green", marker='x')
    
    r1 = get_molecular_contour(length, strandApos, strandBpos, linear=True)
    r1 = r1[buff:-buff,:]
    line1 = Line.best_fit(r1)
    grad1 = line1.direction
    grad1 /= np.linalg.norm(grad1)
    
    strandApos = np.array(strand1.get_atom_list("C1'")[strand1.nres-length:])
    strandBpos = np.array(strand2.get_atom_list("C1'")[strand2.nres-length:])
    for r in [strandApos,strandBpos]:
        for vec in r:
            ax.scatter(vec[0],vec[1],vec[2], color="purple", marker='x')
            


    r2 = get_molecular_contour(length, strandApos, strandBpos, linear=True)
    r2 = r2[buff:-buff,:]
    line2 = Line.best_fit(r2)
    grad2 = line2.direction
    grad2 /= np.linalg.norm(grad2)
    
    
    # ra, grada, linea = r1, grad1, line1
    # rc, gradc, linec = r2, grad2, line2
    
    rc, gradc, linec = r1, grad1, line1
    ra, grada, linea = r2, grad2, line2
    # make sure all vectors are pointing in same direction by checking againt centre to end-point vector   
    
    for end, grad in zip([a_end,b_end,c_end], [grada,gradb,gradc]):
        print(np.dot((end - mid)/np.linalg.norm(end - mid), grad))
        if np.dot(end - mid, grad) < 0:
            grad *= -1
            
    print(grada,gradb,gradc)

    
    cosalpha = np.dot(grada,gradc)
    alpha = np.arccos(cosalpha)*180/np.pi
    
    cosbeta = np.dot(gradb,gradc)
    beta = np.arccos(cosbeta)*180/np.pi
    
    cosphi = np.dot(-grada,gradb)
    phi = np.arccos(cosphi)*180/np.pi 
    
    # use sign convention outlined in Behnam paper
    sign = np.dot( np.cross(grada,gradc), np.cross(gradb,gradc))
    print(np.cross(grada,gradc))
    print(np.cross(gradb,gradc))
    print(sign)
    if sign < 0:
        alpha = abs(alpha)
        beta = abs(beta)
    else:
        # print("blep")
        alpha = - abs(alpha)
        beta =  - abs(beta)
        
    sign = np.dot( np.cross(np.cross(grada,gradc), np.cross(gradb,gradc)), gradc)
    if sign < 0:
        phi = abs(phi)
    else:
        phi = - abs(phi)
    
    
    colour_list = ["blue", "red", "green", "purple"]
    ax.scatter(mid[0],mid[1],mid[2],c="yellow",s=50)
    
    for name, end in zip(["a","b","c"],[a_end, b_end, c_end]):
        ax.scatter(end[0],end[1],end[2], color="black")
        ax.text(end[0],end[1],end[2], name)
    
    for r in [ra,rb,rc]:
        for vec in r:
            ax.scatter(vec[0],vec[1],vec[2], color="red", marker='x')
            
    
            
    for line in [linea,lineb,linec]:
        points = [np.array(line.point + line.direction * s) for s in np.linspace(-25,25,100)]
        for point in points:
            ax.scatter(point[0],point[1],point[2], color="black")
                  
    plt.show()
    
    return alpha, beta, phi

conf_test = traj.last_conf

alpha, beta, phi = T_motif_angles_debug(conf_test)

