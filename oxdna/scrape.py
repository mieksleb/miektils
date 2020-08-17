#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 23:22:19 2019


python script which defines classes which grab info from input and force files and stores 

@author: michaelselby
"""
import re 
import tools

#class of input parameters
class INPUT:
    def __init__(self, input_file):
        self.print_group = []
        for i in range(len((input_file))):
            if input_file[i].startswith("steps"):
                self.steps = re.findall(r"\-?\d*\.?\d*[e]\-?\d*(?!.*=)", input_file[i])[0]
            elif "direction" in input_file[i]:
                self.direction = tools.get_array(re.findall(r"\d.?\d*", input_file[i]))
            elif "origin" in input_file[i]:
                self.origin = tools.get_array(re.findall(r"\d.?\d*", input_file[i]))
            elif "print_group" in input_file[i]:
                self.print_group.append(input_file[i].split()[-1])
            elif "print_conf_interval" in input_file[i]:
                self.print_conf_interval = input_file[i].split()[-1]
                
                
                
                
#class of external torque force parameters
class twist:
    def __init__(self, twist_snippet):
        print
        for i in range(len(twist_snippet)):
            if "particle" in twist_snippet[i]:
                self.particle = int(re.findall(r"\d+", twist_snippet[i])[0])
            elif "pos0" in twist_snippet[i]:
                self.pos0 = tools.get_array(re.findall(r"\-?\d*\.\d*(?!.*=)", twist_snippet[i]))
            elif "stiff" in twist_snippet[i]:
                self.stiff = float(re.findall(r"\-?\d+\.?\d*", twist_snippet[i])[0])
            elif "rate" in twist_snippet[i]:
                self.rate = re.findall(r"\-?\d*\.?\d*[e]\-?\d*(?!.*=)", twist_snippet[i])[0]
            elif "base" in twist_snippet[i]:
                self.base = float(re.findall(r"\-?\d+", twist_snippet[i])[0])
            elif "center" in twist_snippet[i]:
                self.center = tools.get_array(re.findall(r"\-?\d\.?\d*", twist_snippet[i]))
            elif "axis" in twist_snippet[i]:
                self.axis = tools.get_array(re.findall(r"\-?\d\.?\d*", twist_snippet[i]))
            elif "mask" in twist_snippet[i]:
                self.mask = tools.get_array(re.findall(r"\-?\d\.?\d*", twist_snippet[i]))
            elif "group_name" in twist_snippet[i]:
                self.group = twist_snippet[i].split()[-1]
        self.omega = tools.pi2*float(self.rate)
                