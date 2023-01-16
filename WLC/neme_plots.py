#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 20:59:52 2023

@author: michaelselby
"""
import matplotlib.pyplot as plt


lines = open("neme_data.dat", "r").readlines()

keys = lines[0].split()


dict_list = []
for line in lines[1:]:
    quants = line.split()
    dic = {}
    for i, key in enumerate(keys):
        dic[key] = quants[i]
    dict_list.append(dic)


for dic in dict_list:
    
    x = float(dic["Sigma"])
    y = float(dic["Force"])
    
    if int(dic["Nemey"]) == 1:
        col = "red"
        label = "Plectoneme"
        plt.scatter( x , y, color="red")
        last_neme = x , y 
    else:
        col = "blue"
        label = "BDNA"
        plt.scatter( x , y, color="blue")
        last_bdna = x, y
        
plt.scatter( last_neme[0] , last_neme[1], label="Plectoneme", color="red")
plt.scatter( last_bdna[0] , last_bdna[1], label="BDNA", color="blue")
plt.legend()
    

plt.show()