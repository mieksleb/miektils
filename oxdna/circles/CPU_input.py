#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:10:11 2019

@author: michaelselby
"""
import sys

"""
A simple script to generate input files, this is the most general class and other more sepcific versions 
can be loaded in for repeat generations

"""
input_file = open("input.dat", "w")

topology_file = open(sys.argv[1], "r").readlines()
conf_file = open(sys.argv[2], "r").readlines()
trajectory_file = open(sys.argv[3],"w")

class Input:
    def __init__(self, backend="CPU", steps="1e3", diff="2.50", thermostat="john", temp="300K", dt="0.005", salt_conc="0.5", top_file="generate.top", conf_file="generate.dat", last_conf_file="last_conf.dat", traj_file="trajectory.dat", energy_file="energy.dat", force_file="forces.dat"):
        self.backend = str(backend)
        self.steps = str(steps)
        self.diff = str(diff)
        self.thermostat = str(thermostat)
        self.temp = str(temp)
        self.dt = str(dt)
        self.salt_conc = str(salt_conc)
        self.top_file = top_file
        self.conf_file = str(conf_file)
        self.last_conf_file = str(last_conf_file)
        self.traj_file = str(traj_file)
        self.energy_file = str(energy_file)
        self.force_file = str(force_file)
    def __repr__(self):
        return str(   
"####    PROGRAM PARAMETERS  ####"+"\n"+
"backend = "+self.backend+"\n"+
"backend_precision = double"+"\n"+
"#CUDA_list = verlet"+"\n"+
"#use_edge = 1"+"\n"+
"interaction_type = DNA2"+"\n"+
"\n"+
"\n"+
"####    SIM PARAMETERS    ####"+"\n"+
"steps = "+self.steps+"\n"+
"newtonian_steps = 103"+"\n"+
"diff_coeff ="+self.diff+"\n"+
"#pt = 0.1"+"\n"+
"thermostat ="+self.thermostat+"\n"+
"\n"+
"T = "+self.temp+"\n"+
"dt = "+self.dt+"\n"+
"verlet_skin = 0.05"+"\n"+
"\n"+
"salt_concentration = "+self.salt_conc+"\n"+
"\n"+
"####    INPUT / OUTPUT    ####"+"\n"+
"topology = "+self.top_file+"\n"+
"conf_file = "+self.conf_file+"\n"+
"last_conf_file = "+self.last_conf_file+"\n"+
"trajectory_file = "+self.traj_file+"\n"+
"refresh_vel = 1"+"\n"+
"log_file = log.dat"+"\n"+
"no_stdout_energy = 0"+"\n"+
"restart_step_counter = 1"+"\n"+
"energy_file = "+self.energy_file+"\n"
"print_conf_interval = 1e3"+"\n"+
"print_energy_every = 1e3"+"\n"+
"time_scale = linear"+"\n"+
"external_forces = 1"+"\n"+
"external_forces_file = "+self.force_file+"\n"
)

chinput = Input(backend="CPU", steps="1e6", diff=2.50, thermostat="john", top_file = sys.argv[1], conf_file = sys.argv[2], traj_file = sys.argv[3])   

input_file.write(chinput.__repr__())
#print(chinput.__repr__())
input_file.close()
    
