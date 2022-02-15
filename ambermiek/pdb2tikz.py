#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 10:14:35 2021

This script is for producing Tikz pictures from pdb files

@author: michaelselby
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
import tikzplotlib



# pdb_file = open("TT_40/MTTD/dv40tt_md9_nw.pdb","r")
pdb_file = open("TT_40/MDNA/dv40t_md9_nw.pdb","r")
# pdb_file = open("bending/MDNA/short.pdb","r")
# pdb_file = open("bending/MDNA/dv40t_nw.pdb","r")




bp = 40
natoms = 2*bp
base = 20
limit = 100
# f = open("roll.dat","w")
n_color = "green"
c_color = "red"
h_color = "white"
p_color = "orange"
o_color = "blue"



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
    string = "\draw[thick," + color1 + "] "+np2coord(base_dict[atom1]) + " -- " + np2coord(centre)+";"+"\n"
    string += "\draw[thick," + color2 + "] "+np2coord(base_dict[atom2]) + " -- " + np2coord(centre)+";"+"\n"
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
    string = "\draw[thick," + color1 + "] "+np2coord(base_dict1[atom1]) + " -- " + np2coord(centre)+";"+"\n"
    string += "\draw[thick," + color2 + "] "+np2coord(base_dict2[atom2]) + " -- " + np2coord(centre)+";"+"\n"
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
                string += "\\node[circle,inner sep=0pt,minimum size=5pt,fill="+color1+",label=below:{"+key+"}] at "+np2coord(value) +" {};"+"\n"
    
    if base_num==True:
        base_number = base_dict["base_number"]
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
    resname  = base_dict["resname"]
    if resname=="DT":
        string = texline(base_dict,"P","O5'")
        string += texline(base_dict,"P","OP2")
        string += texline(base_dict,"P","OP1")
        string += texline(base_dict,"P","C5'")
        string += texline(base_dict,"C5'","C4'") 
        string += texline(base_dict,"C4'","C3'")
        string += texline(base_dict,"C3'","C2'")
        string += texline(base_dict,"C3'","O3'")
        string += texline(base_dict,"C2'","C1'")
        string += texline(base_dict,"C1'","N1")
        string += texline(base_dict,"C1'","O4'")
        string += texline(base_dict,"C4'","O4'")
        string += texline(base_dict,"N1","C6")
        string += texline(base_dict,"C6","C5")
        string += texline(base_dict,"C5","C7")
        string += texline(base_dict,"C5","C4")
        string += texline(base_dict,"C4","N3")
        string += texline(base_dict,"C4","O4")
        string += texline(base_dict,"N3","C2")
        string += texline(base_dict,"C2","O2")
        string += texline(base_dict,"C2","N1")
    elif resname=="TT":
        string = texline(base_dict,"P","O5'")
        string += texline(base_dict,"P","OP2")
        string += texline(base_dict,"P","OP1")
        string += texline(base_dict,"P","C5'")
        string += texline(base_dict,"C5'","C4'")  
        string += texline(base_dict,"C4'","C3'")
        string += texline(base_dict,"C3'","C2'")
        string += texline(base_dict,"C3'","O3'")
        string += texline(base_dict,"C2'","C1'")
        string += texline(base_dict,"C1'","N1")
        string += texline(base_dict,"C1'","O4'")
        string += texline(base_dict,"C4'","O4'")
        string += texline(base_dict,"N1","C6")
        string += texline(base_dict,"C6","C5")
        string += texline(base_dict,"C5","C7")
        string += texline(base_dict,"C5","C4")
        string += texline(base_dict,"C4","N3")
        string += texline(base_dict,"C4","O4")
        string += texline(base_dict,"N3","C2")
        string += texline(base_dict,"C2","O2")
        string += texline(base_dict,"C2","N1")
    elif resname=="DA":
        string = texline(base_dict,"P","O5'")
        string += texline(base_dict,"P","OP2")
        string += texline(base_dict,"P","OP1")
        string += texline(base_dict,"P","C5'")
        string += texline(base_dict,"C5'","C4'")  
        string += texline(base_dict,"C4'","C3'")
        string += texline(base_dict,"C3'","C2'")
        string += texline(base_dict,"C3'","O3'")
        string += texline(base_dict,"C2'","C1'")
        string += texline(base_dict,"O4'","C1'")
        string += texline(base_dict,"O4'","C4'")
        string += texline(base_dict,"C3'","O3'")
        string += texline(base_dict,"C1'","N9")
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
            
        
    return string

        


fig2 = plt.figure()

ax = fig2.add_subplot(111, projection='3d')

step = 0

base_dict_list = []
base_dict = {}
base_dict = dict.fromkeys(["base_number","resname","P","OP1","OP2","O5'","C5'","H5'","H5''","N9","HO5'","N7","N6",\
                            "C4'","H4'","O4'","C1'","H1'","N1","H61","H62","O5'","H5'","C8","H8", \
                            "C6","H6","C5","C7","H71","H72","H73","C4","O4","N3","H3","C2","O2",\
                                "C3'","H3'","C2'","H2'","H2''","O3'"])
bpimatch = 1
for line in pdb_file:
    split_line = line.split()
    word = split_line[0]
    if (word=="MODEL"):
        count = 0
        frame = int(split_line[1])
    elif (word=="ATOM"):
        atom = split_line[2]
        resname = split_line[3]
        bpi = int(split_line[4])
        x,y,z = float(split_line[5]),float(split_line[6]),float(split_line[7])


        
        if bpi==bpimatch:
            base_dict["resname"] = resname
            base_dict["base_number"] = bpi
            keys = base_dict.keys()
            if atom in keys:
                base_dict[atom] = np.array([x,y,z])
                
            
        else:

            base_dict_list.append(base_dict)
            base_dict = {}
            base_dict = dict.fromkeys(["base_number","resname","P","OP1","OP2","O5'","C5'","H5'","HO5'","H5''","N9","N7","N6",\
                                    "C4'","H4'","O4'","C1'","H1'","N1","H61","H62","O5'","H5'","C8","H8", \
                                    "C6","H6","C5","C7","H71","H72","H73","C4","O4","N3","H3","C2","O2",\
                                        "C3'","H3'","C2'","H2'","H2''","O3'"])
            base_dict["resname"] = resname
            keys = base_dict.keys()
            if atom in keys:
                base_dict[atom] = np.array([x,y,z])
            bpimatch = bpi
        

        
        
        
        
    else:
        if (word=="ENDMDL"):   
            base_dict_list.append(base_dict)  

            step += 1
            # if step>limit:
            break
            
        
            
                
"""
Now we have all the information we can go to work

"""

smallx = 0
smally = 0
smallz = 0

for dic in base_dict_list:
    for key,value in dic.items():
        if isinstance(value,np.ndarray):
            x,y,z = value[0],value[1],value[2]
        if x < smallx:
            smallx = x
        if y < smally:
            smally = y
        if y < smallz:
            smallz = z

num = 20
origin = copy.deepcopy(base_dict_list[num]["O3'"])
for dic in base_dict_list:
    for key,value in dic.items():
        if isinstance(value,np.ndarray):
            dic[key] -= origin
            
  

# in tikz, the z-axis points out of the screen, the y-axis points upwards and x-axis points left to right          
# we wish to rotate the coordinates such that your chosen axis becomes y-axis
            
vec = base_dict_list[60]["P"]-base_dict_list[61]["P"]
vec /= np.linalg.norm(vec)
axis = np.array([0,1,0])
angle = np.arccos(np.dot(axis,vec))
cross = np.cross(axis,vec)
cross /= np.linalg.norm(cross)

for dic in base_dict_list:
    for key,value in dic.items():
        if isinstance(value,np.ndarray):
            dic[key] = rot(-angle,cross,dic[key])


            
            
   
tex_file = open("dna.tex","w")
tex_file.write("\\tdplotsetmaincoords{0}{0}"+"\n"\
               +"\\tdplotsetrotatedcoords{0}{0}{0}" +"\n"+\
                   "\\begin{tikzpicture}[tdplot_rotated_coords]"+"\n")
         
            
# select the base you wish to turn into tex
# for i in range(2,bp-1):
#     my_base_dict = base_dict_list[i]
#     stringus = base2tex(my_base_dict)
#     tex_file.write(stringus)
    


old_base_dict = {}
base_min = 34
base_max = 35
base_list = [34,46]
for i in range(0,len(base_dict_list)):

    my_base_dict = base_dict_list[i]
    
    if my_base_dict["base_number"] in base_list:
        stringus = base2tex(my_base_dict)
        pingus = texlabel(my_base_dict)
        tex_file.write(stringus)
        tex_file.write(pingus)
    
    # if i>base_min:
    #     bond_string = texline_different(old_base_dict,my_base_dict,"O3'","P")
    #     tex_file.write(bond_string)
    #     if old_base_dict["resname"]=="TT" and my_base_dict["resname"]=="TT":
    #         print('chubs')
    #         bond_string1 = texline_different(old_base_dict,my_base_dict,"C5","C5")
    #         tex_file.write(bond_string1)
    #         bond_string2 = texline_different(old_base_dict,my_base_dict,"C6","C6")
    #         tex_file.write(bond_string2)
        
    old_base_dict = my_base_dict



                    
            
            
            
            
pdb_file.close()
tex_file.write("\end{tikzpicture}")
tex_file.close()


