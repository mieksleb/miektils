#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:40:00 2019

A python script which will define all the classes to be used with miektils

@author: michaelselby
"""
import tools


    
class Input:
    def __init__(self, backend="CPU", steps="1e3", diff="2.50", thermostat="john", temp="300K", dt="0.05", salt_conc="1", top_file="generate.top", conf_file="generate.dat", last_conf_file="last_conf.dat", traj_file="trajectory.dat", energy_file="energy.dat", force_file="forces.dat"):
        self.backend = str(backend)
        self.steps = str(steps)
        self.diff = str(diff)
        self.thermostat = str(thermostat)
        self.temp = str(temp)
        self.dt = str(dt)
        self.salt_conc = str(salt_conc)
        self.top_file = str(top_file)
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
        
        
        
        
        
class ForceTwist:
	def __init__(self, particle, pos0, stiff, rate, base, axis, center, mask, group_name):
		self.particle = str(particle)
		self.pos0 = tools.vector_string(pos0)
		self.stiff = str(stiff)
		self.rate = str(rate)
		self.base = str(base)
		self.axis = tools.vector_string(axis)
		self.center = tools.vector_string(center)
		self.mask = tools.vector_string(mask)
		self.group_name = str(group_name)
	def __str__(self):
		return str(
"{"+"\n"+
"  type = twist\n"+
"  particle = "+self.particle+"\n"+
"  pos0 = "+self.pos0+"\n"+
"  stiff = "+self.stiff+"\n"+ 
"  rate = "+self.rate+"\n"+
"  base = "+self.base+"\n"+
"  axis = "+self.axis+"\n"+
"  center = "+self.center+"\n"+
"  mask = "+self.mask+"\n"+
"  group_name = "+self.group_name+"\n"+
"}"+"\n"+"\n")  
                       
    

    
class ForceString:
	def __init__(self,particle,F0="1.",rate="0.", direction = [1,0,0]):
		self.particle = str(particle)
		self.F0 = str(F0) #1 [SU] = 48.63 pN in oxDNA = 49.3 pN in oxRNA
		self.rate = str(rate)
		self.dir = tools.vector_string(direction)
	def __repr__(self):
		return str(
"{\n"+
"  type = string\n"+
"  particle = "+self.particle+"\n"+
"  F0 = "+self.F0+"\n"+
"  rate = "+self.rate+"\n"+
"  dir = "+self.dir+"\n"+
"}\n\n")  
    
    
    
    
    
class ForceTrap:
	def __init__(self,particle,pos0,stiff="1.0",rate="0.",direction="1., 0., 0.",group_name = []):
		self.particle = str(particle)
		self.pos0 = str(pos0)
		self.stiff = str(stiff)
		self.rate = str(rate)
		self.dir = str(direction)
		self.group_name = str(group_name)
	def __repr__(self):
		return str(
"{\n"+
"  type = trap\n"+
"  particle = "+self.particle+"\n"+
"  pos0 = "+self.pos0+"\n"+
"  stiff = "+self.stiff+"\n"+ #1 = 48.63 pN/[SU Length] in OxDNA = 49.3 pN/[SU Length] in OxRNA
"  rate = "+self.rate+"\n"+
"  dir = "+self.dir+"\n"+
"  group_name = "+self.group_name+"\n" * (self.group_name != str([]))+
"}")  

class ForcePlaneMoving:
	def __init__(self,particle,stiff="1.0",direction="1., 0., 0.",ref_particle="0"):
		self.particle = str(particle)
		self.stiff = str(stiff)
		self.dir = str(direction)
		self.ref_particle = str(ref_particle)
	def __repr__(self):
		return str(
"{\n"+
"  type = repulsion_plane_moving\n"+
"  particle = "+self.particle+"\n"+
"  stiff = "+self.stiff+"\n"+ #1 = 48.63 pN/[SU Length] in OxDNA model = 49.3 pN/[SU Length] in OxRNA
"  dir = "+self.dir+"\n"+
"  ref_particle = "+self.ref_particle+"\n"+
"}")  
    
    