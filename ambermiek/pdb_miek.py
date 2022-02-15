#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 11:40:16 2022

This python script defines the classes trajectory and configuration

@author: michaelselby
"""

import numpy as np

list1 = ["resname","P","OP1","OP2","O5'","C5'","H5'","HO5'","H5''","N9","N7","N6",\
                                            "C4'","H4'","O4'","C1'","H1'","N1","H61","H62","O5'","H5'","C8","H8", \
                                            "C6","H6","C5","C7","H71","H72","H73","C4","O4","N3","H3","C2","O2",\
                                                "C3'","H3'","C2'","H2'","H2''","O3'"]


class trajectory:
    def __init__(self, file_name, circular):
        self.file_name = file_name
        # self.get_steps()
        self.circular = circular
        
    def get_steps(self):
        lines = open(self.file_name,"r").readlines()
        rev_lines = reversed(lines)
        for line in rev_lines:
            split_line = line.split()
            if split_line[0]=="MODEL":
                steps = int(split_line[1])
                break
        self.steps = steps
        
       
    def add_configurations(self,max_steps=100000):
        """
        the method appends configuration objects to the config_list

        Returns
        -------
        None.

        """
        self.config_list = []
        with open(self.file_name,"r") as file:
            for line in file:
                # base_dict = {}
                # base_dict = dict.fromkeys(list1)
                split_line = line.split()
                word = split_line[0]
                if word=="MODEL":
                    bimatch = 1
                    step = int(split_line[1])
                    if step > max_steps:
                        break
                    base_dict = {}
                    base_dict = dict.fromkeys(list1)
                    config = configuration(step) # initalize a configuration class for new time step
                    strand1 = strand(self.circular)
                elif (word=="ATOM"):
                    atom = split_line[2]
                    resname = split_line[3]
                    bi = int(split_line[4])
                    x,y,z = float(split_line[5]),float(split_line[6]),float(split_line[7])
                    
                    if bi==bimatch:   # if base index matches, add atom to the base dictionary
                        base_dict["resname"] = resname
                        keys = base_dict.keys()
                        if atom in keys:
                            base_dict[atom] = np.array([x,y,z])
      
                    else:                 
                        base1 = base(base_dict)
                        strand1.add_base(base1) # add base to the configuration class
                        base_dict = {}          # reset the base dictionary
                        base_dict = dict.fromkeys(list1)
                        base_dict["resname"] = resname
                        keys = base_dict.keys()
                        if atom in keys:
                            base_dict[atom] = np.array([x,y,z])
                        bimatch = bi # increments the base index by one
                elif word=="TER":
                    base1 = base(base_dict)
                    strand1.add_base(base1) # add base to the configuration class
                    config.add_strand(strand1)             
                    strand1 = strand(self.circular)
                    bimatch += 1 # this just means next strand doesn't pick up last residue of last strand
                elif word=="ENDMDL":
                    self.config_list.append(config)
                else:
                    pass

                


class strand:
    def __init__(self,circular):
        self.base_list = [] # list of base classes
        self.nres = 0
        self.circular = circular
        
    def add_base(self,base):
        self.base_list.append(base)
        self.nres += 1
        
    def get_bb_list(self):
        self.bb_list = []
        for base in self.base_list:
            self.bb_list.append(base.bb)
            
    def get_base_centres(self):
        self.base_centres = []
        for base in self.base_list:
            base.get_base_centre()
            self.base_centres.append(base.base_centre)
        
        
        

        

class base:
    def __init__(self,base_dict):
        self.base_dict = base_dict
        self.bb = self.base_dict["P"]
        self.resname = self.base_dict["resname"]
        
    def get_base_centre(self):
        n1,c6,c5,c4,n3,c2 = self.base_dict["N1"],self.base_dict["C6"],self.base_dict["C5"],self.base_dict["C4"],self.base_dict["N3"],self.base_dict["C2"]
        self.base_centre = (n1+c6+c5+c4+n3+c2)/6

            

        
class configuration:
    '''
    decribes a single time-step configuration
    '''
    def __init__(self,step):
        self.step = step
        self.strand_list = [] # list of base classes
        self.nstrands = 0  # number of strands
        
    def add_strand(self,strand):
        self.strand_list.append(strand)
        self.nstrands +=1
        
       

    def get_mid_list(self):
        '''
        Returns
        -------
        mid_list, a list of the midpoints of a duplex. Can only be called after get_bb_list()

        '''
        self.mid_list = []
        strand1 = self.strand_list[0]
        strand2 = self.strand_list[1]
        strand1.get_bb_list()
        strand2.get_bb_list()
        for i in range(strand1.nres):
            mid = (strand1.bb_list[i] + strand2.bb_list[-i]) /2
            self.mid_list.append(mid)
        
        
        # for i in range(strand1.nres):
        #     mid = 0
        #     for strand in self.strand_list:
        #         strand.get_bb_list()
        #         mid += strand.bb_list[i]
        #     mid /= self.nstrands
        #     self.mid_list.append(mid)

    def get_duplex_centre(self, helical=True,alpha=0):
        '''
        Returns
        -------
        mid_list, a list of the midpoints of a duplex. Can only be called after get_bb_list()

        '''
        self.mid_list = []
        strand1 = self.strand_list[0]
        strand2 = self.strand_list[1]
        strand1.get_base_centres()
        strand2.get_base_centres()
        for i in range(strand1.nres):
            mid = (strand1.base_centres[i] + strand2.base_centres[-i]) /2
            self.mid_list.append(mid)

        
