#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 23:22:19 2019


python script which defines classes which grab info from input and force files and stores 

@author: michaelselby
"""
import re 
import numpy as np

def get_array(line):                #turns a line in a trajectpory file into a numpy array
    return np.array([line[0][:-1], line[1][:-1], line[2]])

float_re = r"\-?\d+\.?\d*"
float_file = r"\-?\w+\.?\w*"

#class of input parameters
class oxdnaINPUT:    
    def scrape(self, input_file): 
        self.print_group = []
        for i in range(len((input_file))):
            if input_file[i].startswith("steps"):
                self.steps = re.findall(r"\-?\d*\.?\d*[e]\-?\d*(?!.*=)", input_file[i])[0]              
            elif "T" in input_file[i]:
                self.temperature = re.findall(float_re,input_file[i])
            elif "dt" in input_file[i]:
                self.dt = re.findall(float_re,input_file[i])
            elif "verlet skin" in input_file[i]:
                self.verlet_skin = re.findall(float_re,input_file[i])
            elif "salt_conc" in input_file[i]:
                self.salt_conc = re.findall(float_re,input_file[i])
            elif "direction" in input_file[i]:
                self.direction = get_array(re.findall(r"\d.?\d*", input_file[i]))
            elif "origin" in input_file[i]:
                self.origin = get_array(re.findall(r"\d.?\d*", input_file[i]))
            elif "print_group" in input_file[i]:
                self.print_group.append(input_file[i].split()[-1])
            elif "print_conf_interval" in input_file[i]:
                self.print_conf_interval = input_file[i].split()[-1]
            elif "topology" in input_file[i]:
                self.top_file = re.findall(float_re,input_file[i])
            elif "conf_file" in input_file[i]:
                self.init_conf_file = re.findall(float_re,input_file[i])
            elif "last_conf_file" in input_file[i]:
                self.final_conf_file = re.findall(float_re,input_file[i])


    def write(self, input_file_name):
        input_file = open(input_file_name,"w")
        self.print_group = []
        for i in range(len((input_file))):
            if input_file[i].startswith("steps"):
                self.steps = re.findall(r"\-?\d*\.?\d*[e]\-?\d*(?!.*=)", input_file[i])[0]              
            elif "T" in input_file[i]:
                self.temperature = re.findall(float_re,input_file[i])
            elif "dt" in input_file[i]:
                self.dt = re.findall(float_re,input_file[i])
            elif "verlet skin" in input_file[i]:
                self.verlet_skin = re.findall(float_re,input_file[i])
            elif "salt_conc" in input_file[i]:
                self.salt_conc = re.findall(float_re,input_file[i])
            elif "direction" in input_file[i]:
                self.direction = get_array(re.findall(r"\d.?\d*", input_file[i]))
            elif "origin" in input_file[i]:
                self.origin = get_array(re.findall(r"\d.?\d*", input_file[i]))
            elif "print_group" in input_file[i]:
                self.print_group.append(input_file[i].split()[-1])
            elif "print_conf_interval" in input_file[i]:
                self.print_conf_interval = input_file[i].split()[-1]
                
                
     
        
class amberINPUT:
    def __init__(self,dt=None,temperature=None,salt_conc=None,steps=None,init_temp=None,rest_wt=None,restraint_mask=None,ntx=None,irest=None):
        self.dt = dt
        self.temperature = temperature
        self.salt_conc = salt_conc
        self.steps = steps
        self.init_temp = init_temp
        self.restraint_mask = restraint_mask
        self.rest_wt = rest_wt
        self.ntx = ntx
        self.irest=irest
    def write_md(self,input_file_name): 
            input_file = open(input_file_name,"w")
            input_file.write(
                "MD on the water and ions about the DNA"+"\n"
                " &cntrl"+"\n"
                "\t"+"ntf=2, ntc=2, ntb=0, ,cut=1000.,"+"\n"
                "\t"+"nstlim="+str(self.steps)+", dt="+str(self.dt)+", ntpr=5000,"+"\n"
                "\t"+"ntwx=5000,ntwe=5000,ntwr=5000,"+"\n"
                "\t"+"temp0="+str(self.temperature)+", ntt=3,gamma_ln=0.01,"+"\n"
                "\t"+"imin=0, irest="+str(self.irest)+", ntx="+str(self.ntx)+","+"\n"
                "\t"+"igb=8, gbsa=0, saltcon="+str(self.salt_conc)+","+"\n"
                "\t"+"ntr=0,"+"\n"
                " &end"
                )
            input_file.close()
    def write_equib(self,input_file_name): 
            input_file = open(input_file_name,"w")
            input_file.write(
                "MD (equilibration) on the water and ions about the DNA"+"\n"
                " &cntrl"+"\n"
                "\t"+"ntf=2, ntc=2, ntb=0, ,cut=1000.,"+"\n"
                "\t"+"nstlim="+str(self.steps)+", dt="+str(self.dt)+", ntpr=5000,"+"\n"
                "\t"+"tempi="+str(self.init_temp)+", temp0="+str(self.temperature)+","+"\n"
                "\t"+"ntt=3,gamma_ln=0.01,"+"\n"
                "\t"+"imin=0, irest="+str(self.irest)+", ntx="+str(self.ntx)+","+"\n"
                "\t"+"igb=8, gbsa=0, saltcon="+str(self.salt_conc)+","+"\n"
                "\t"+"ntr=1,"+"\n"
                "\t"+"restraint_wt="+str(self.rest_wt)+","+"\n"
                "\t"+"restraintmask="+str(self.restraint_mask)+"\n"
                " &end"
                )
            input_file.close()
    def write_min(self,input_file_name): 
            input_file = open(input_file_name,"w")
            input_file.write(
                "Initial minimization of the whole system"+"\n"
                " &cntrl"+"\n"
                "  imin=1,"+"\n"
                "  maxcyc=10000,"+"\n"
                "  ncyc=2500,"+"\n"
                "  ntb=0,"+"\n"
                "  ntp=0,"+"\n"
                "  ntr=0,"+"\n"
                "  igb=1,"+"\n"
                " /"
                )
            input_file.close()



        
class_test = amberINPUT(dt=0.002,temperature=300,salt_conc=0.2, steps = 50000).write(input_file_name="chumble.dat")       

        
        
             
#class of external torque force parameters
class twist:
    def __init__(self, twist_snippet):
        for i in range(len(twist_snippet)):
            if "particle" in twist_snippet[i]:
                self.particle = int(re.findall(r"\d+", twist_snippet[i])[0])
            elif "pos0" in twist_snippet[i]:
                self.pos0 = get_array(re.findall(r"\-?\d*\.\d*(?!.*=)", twist_snippet[i]))
            elif "stiff" in twist_snippet[i]:
                self.stiff = float(re.findall(r"\-?\d+\.?\d*", twist_snippet[i])[0])
            elif "rate" in twist_snippet[i]:
                self.rate = re.findall(r"\-?\d*\.?\d*[e]\-?\d*(?!.*=)", twist_snippet[i])[0]
            elif "base" in twist_snippet[i]:
                self.base = float(re.findall(r"\-?\d+", twist_snippet[i])[0])
            elif "center" in twist_snippet[i]:
                self.center = get_array(re.findall(r"\-?\d\.?\d*", twist_snippet[i]))
            elif "axis" in twist_snippet[i]:
                self.axis = get_array(re.findall(r"\-?\d\.?\d*", twist_snippet[i]))
            elif "mask" in twist_snippet[i]:
                self.mask = get_array(re.findall(r"\-?\d\.?\d*", twist_snippet[i]))
            elif "group_name" in twist_snippet[i]:
                self.group = twist_snippet[i].split()[-1]

                