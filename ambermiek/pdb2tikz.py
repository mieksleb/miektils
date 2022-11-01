#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 10:14:35 2021

This script is for producing Tikz pictures from pdb files

@author: michaelselby
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import copy
import tikzplotlib
import pdb_miek

sys.path.append('/Users/michaelselby/Documents/DPhil/miektils/ambermiek/linear/')

# pdb_file = open("TT_40/MTTD/dv40tt_md9_nw.pdb","r")
# pdb_file = open("TT_40/MDNA/dv40t_md9_nw.pdb","r")
pdb_file = "/Users/michaelselby/Documents/DPhil/miektils/ambermiek/linear/dv22t.pdb"
# pdb_file = open("bending/MDNA/short.pdb","r")
# pdb_file = open("bending/MDNA/dv40t_nw.pdb","r")


pdb = open(pdb_file,"r")
traj = pdb_miek.trajectory(pdb_file,circular=False)
traj.add_configurations(scale=0.1)
config = traj.config_list[0]
strand1 = config.strand_list[0]
strand2 = config.strand_list[1]
base_list1 = strand1.base_list
base_list2 = strand2.base_list


bp = 22
natoms = 2*bp
base = 20
limit = 100
n_color = "blue"
c_color = "cyan"
h_color = "white"
p_color = "orange"
o_color = "red"



def rot(theta, axis, vec):  #rotates vector B about axis A using Euler-Rodrigues formula
    return vec*np.cos(theta) + np.sin(theta)*np.cross(axis, vec) + axis*(np.dot(axis,vec))*(1-np.cos(theta))


def np2coord(arr):
    return "("+str(arr[0])+","+str(arr[1])+","+str(arr[2])+")" 

def coords2string(arrays):
    string = np2coord(arrays[0])
    for i in range(1,len(arrays)):
        string += " --  " + np2coord(arrays[i])
    return string

def texline(base_dict,atom1,atom2):
    """
    Given 2 atom names, will draw the bond between atom1 and atom2, colouring 
    the bond half and half depending on the atom types
    """
    type1 = list(atom1)[0]
    type2 = list(atom2)[0]

    if type1=="C":
        color1 = c_color
    elif type1=="H":
        color1 = h_color
    elif type1=="N":
        color1 = n_color
    elif type1=="P":
        color1 = p_color
    elif type1=="O":
        color1 = o_color
    if type2=="C":
        color2 = c_color
    elif type2=="H":
        color2 = h_color
    elif type2=="N":
        color2 = n_color
    elif type2=="P":
        color2 = p_color
    elif type2=="O":
        color2 = o_color

    centre = (base_dict[atom1]+base_dict[atom2])/2
    # string = "\draw[" + color1 + "] "+np2coord(base_dict[atom1]) + " -- " + np2coord(centre)+";"+"\n"
    # string += "\draw[" + color2 + "] "+np2coord(base_dict[atom2]) + " -- " + np2coord(centre)+";"+"\n"
    
    string = "\draw[line width=0.1mm," + color1 + "] " +np2coord(base_dict[atom1]) + " -- " + np2coord(centre)+";"+"\n"
    string += "\draw[line width=0.1mm," + color2 + "] " +np2coord(base_dict[atom2]) + " -- " + np2coord(centre)+";"+"\n"
    string += "\\shade[ball color="+color1+"] "+np2coord(base_dict[atom1]) +" circle [radius=1pt];\n"
    string += "\\shade[ball color="+color2+"] "+np2coord(base_dict[atom2]) +" circle [radius=1pt];\n"
    # string += "\\node[circle,inner sep=0pt,minimum size=2pt,fill="+color1+"] at "+np2coord(base_dict[atom1]) +"{};"+"\n"
    # string += "\\node[circle,inner sep=0pt,minimum size=2pt,fill="+color2+"] at "+np2coord(base_dict[atom2])+"{};"+"\n"
    return string

def texline_different(base_dict1,base_dict2,atom1,atom2):
    """
    Given 2 atom names, will draw the bond between atom1 and atom2, colouring 
    the bond half and half depending on the atom types
    """
    type1 = list(atom1)[0]
    type2 = list(atom2)[0]

    if type1=="C":
        color1 = c_color
    elif type1=="H":
        color1 = h_color
    elif type1=="N":
        color1 = n_color
    elif type1=="P":
        color1 = p_color
    elif type1=="O":
        color1 = o_color
    if type2=="C":
        color2 = c_color
    elif type2=="H":
        color2 = h_color
    elif type2=="N":
        color2 = n_color
    elif type2=="P":
        color2 = p_color
    elif type2=="O":
        color2 = o_color

    centre = (base_dict1[atom1]+base_dict2[atom2])/2
    string = "\draw[line width=0.1mm," + color1 + "] "+np2coord(base_dict1[atom1]) + " -- " + np2coord(centre)+";"+"\n"
    string += "\draw[line width=0.1mm," + color2 + "] "+np2coord(base_dict2[atom2]) + " -- " + np2coord(centre)+";"+"\n"
    return string

def texlabel(base_dict, base_num = True):
    """
    Given a base directory, will label each atom, and optionally the base number if base_num = True
    """
    string = "" 
    for key,value in base_dict.items():
        if isinstance(value,np.ndarray):
            type1 = list(key)[0]
            if type1=="C":
                color1 = c_color
            elif type1=="H":
                color1 = h_color
            elif type1=="N":
                color1 = n_color
            elif type1=="P":
                color1 = p_color
            elif type1=="O":
                color1 = o_color
                
            # only draw non-hydrogen atoms
            if type1!="H":         
                string += "\\node[circle,inner sep=0pt,minimum size=20pt,fill="+color1+",label=below:{"+key+"}] at "+np2coord(value) +" {};"+"\n"
    
    if base_num==True:
        base_number = base_dict["res_number"]
        points = [] 
        if base_dict["resname"]=="DT" or base_dict["resname"]=="TT" or base_dict["resname"]=="DC":  
            Ydict = base_dict
            # add points which are part of pyrimidine ring
            points.append(Ydict["N1"])
            points.append(Ydict["C6"])
            points.append(Ydict["C5"])
            points.append(Ydict["C4"])
            points.append(Ydict["N3"])
            points.append(Ydict["C2"])

        else:
            Rdict = base_dict
            # add points which are part of purine ring
            points.append(Rdict["N1"])
            points.append(Rdict["C2"])
            points.append(Rdict["N3"])
            points.append(Rdict["C6"])
            points.append(Rdict["C5"])
            points.append(Rdict["C4"])

            
        com = 0
        for point in points:
            com += point
        com /= len(points)
        string += "\\node[label=below:{"+str(base_number)+"}] at "+np2coord(com) +" {};"+"\n"
        
    return string


def base2tex(base_dict):
    threetofive = False
    resname  = base_dict["resname"]
    split = list(resname)
    res_type = split[1]
    if len(split) > 2:
        terminal = True
        term_num = split[2]
    
    # add sugar-phosphate backbone
    string = texline(base_dict,"P","OP2")
    string += texline(base_dict,"P","OP1")
    string += texline(base_dict,"C5'","C4'") 
    string += texline(base_dict,"C4'","C3'")
    string += texline(base_dict,"C3'","C2'")
    string += texline(base_dict,"C3'","O3'")
    string += texline(base_dict,"C2'","C1'")
    string += texline(base_dict,"C1'","O4'")
    string += texline(base_dict,"C4'","O4'")
    
    if res_type=="T":
        string += texline(base_dict,"N1","C1'")
        string += texline(base_dict,"N1","C6")
        string += texline(base_dict,"C6","C5")
        string += texline(base_dict,"C5","C7")
        string += texline(base_dict,"C5","C4")
        string += texline(base_dict,"C4","N3")
        string += texline(base_dict,"C4","O4")
        string += texline(base_dict,"N3","C2")
        string += texline(base_dict,"C2","O2")
        string += texline(base_dict,"C2","N1")
        string += texline(base_dict,"P","C5'")  
        
    elif res_type=="A":
        
        string += texline(base_dict,"N9","C8")
        string += texline(base_dict,"N9","C4")
        string += texline(base_dict,"N7","C8")
        string += texline(base_dict,"N7","C5")
        string += texline(base_dict,"C4","C5")
        string += texline(base_dict,"C6","C5")
        string += texline(base_dict,"C6","N6")
        string += texline(base_dict,"C6","N1")
        string += texline(base_dict,"C2","N1")
        string += texline(base_dict,"C2","N3")
        string += texline(base_dict,"C4","N3")
        string += texline(base_dict,"P","C5'")
        string += texline(base_dict,"C1'","N9")
        
    elif res_type=="G":
        string += texline(base_dict,"P","C5'")   
        string += texline(base_dict,"C1'","N9")
        string += texline(base_dict,"C4","N9")
        string += texline(base_dict,"C8","N9")
        string += texline(base_dict,"C8","N7")
        string += texline(base_dict,"C5","N7")
        string += texline(base_dict,"C6","N1")
        string += texline(base_dict,"C6","O6")
        string += texline(base_dict,"C5","C6")
        string += texline(base_dict,"C4","C5")
        string += texline(base_dict,"N3","C4")
        string += texline(base_dict,"C2","N1")
        string += texline(base_dict,"C2","N3")
        string += texline(base_dict,"C2","N2")
    
    elif res_type=="C":
        string += texline(base_dict,"P","C5'") 
        string += texline(base_dict,"C1'","N1")
        string += texline(base_dict,"N1","C6")
        string += texline(base_dict,"C6","C5")
        string += texline(base_dict,"C5","C4")
        string += texline(base_dict,"C4","N4")
        string += texline(base_dict,"C4","N3")
        string += texline(base_dict,"N3","C2")
        string += texline(base_dict,"C2","O2")
        string += texline(base_dict,"C2","N1")

            
    return string

        


fig2 = plt.figure()

ax = fig2.add_subplot(111, projection='3d')

step = 0

            
        
            
                
"""
Now we have all the information we can go to work

"""

# smallx = 0
# smally = 0
# smallz = 0

# for dic in base_dict_list:
#     for key,value in dic.items():
#         if isinstance(value,np.ndarray):
#             x,y,z = value[0],value[1],value[2]
#         if x < smallx:
#             smallx = x
#         if y < smally:
#             smally = y
#         if y < smallz:
#             smallz = z

# num = 20
# origin = copy.deepcopy(base_dict_list[num]["O3'"])
# for dic in base_dict_list:
#     for key,value in dic.items():
#         if isinstance(value,np.ndarray):
#             dic[key] -= origin
            
  

# in tikz, the z-axis points out of the screen, the y-axis points upwards and x-axis points left to right          
# we wish to rotate the coordinates such that your chosen axis becomes y-axis
            
# vec = base_dict_list[60]["P"]-base_dict_list[61]["P"]
vec = np.array([0,0,1])
# vec /= np.linalg.norm(vec)
axis = np.array([0,1,0])
angle = np.arccos(np.dot(axis,vec))
cross = np.cross(axis,vec)
# cross /= np.linalg.norm(cross)

# for dic in base_dict_list:
#     for key,value in dic.items():
#         if isinstance(value,np.ndarray):
#             dic[key] = rot(-angle,cross,dic[key])


            
            
   
tex_file = open("dna.tex","w")
tex_file.write("\\documentclass[10pt]{article}"+ "\n"+\
"\\usepackage[usenames]{color} %used for font color"+ "\n"+\
"\\usepackage{amssymb} %maths"+ "\n"+\
"\\usepackage{amsmath} %maths"+ "\n"+\
"\\usepackage{tikz}"+ "\n"+\
"\\usepackage{pgfplots}"+ "\n"+\
"\\usepackage{tikz-3dplot}"+ "\n"+\
"\\tdplotsetmaincoords{0}{0}"+"\n"+\
"\\tdplotsetrotatedcoords{90}{90}{90}" +"\n"+\
"\\begin{document}"+"\n"+\
"\\tikzset{font=\\fontsize{2}{2.5}\selectfont}"+"\n"+\
"\\begin{tikzpicture}[tdplot_rotated_coords,scale=1]"+"\n")
         
            
# select the base you wish to turn into tex
# for i in range(2,bp-1):
#     my_base_dict = base_dict_list[i]
#     stringus = base2tex(my_base_dict)
#     tex_file.write(stringus)
    


base_min = 0
join_prev = False
base_dict_old = {}
for base in base_list1:
    stringus = ""
    base_dict = base.base_dict
    chum = base_dict["res_number"]
    if chum==1 or chum==22:
        0
    else:
        if join_prev:
            stringus += texline_different(base_dict,base_dict_old,"P","O3'")
        stringus += base2tex(base_dict)
        pingus = texlabel(base_dict)
        tex_file.write(stringus)
        # tex_file.write(pingus)
        join_prev = True
        base_dict_old = base_dict
        
join_prev = False
base_dict_old = {}     
for base in base_list2:
    stringus=""
    base_dict = base.base_dict
    chum = base_dict["res_number"]
    if chum==44 or chum==23:
        0
    else:
        if join_prev:
            stringus += texline_different(base_dict,base_dict_old,"P","O3'")
        stringus += base2tex(base_dict)
        pingus = texlabel(base_dict)
        tex_file.write(stringus)
        # tex_file.write(pingus)
        join_prev = True
        base_dict_old = base_dict




pdb.close()
tex_file.write("\end{tikzpicture}")
tex_file.write("\end{document}")
tex_file.close()


