#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 11:40:16 2022

This python script defines the classes trajectory and configuration

@author: michaelselby
"""

import numpy as np
import tools
import tools_fast_math 
import miek_tools
import copy
import re



list1 = ["resname","res_number","P","OP1","OP2","O5'","C5'","H5'","HO5'","H5''","N9","N7","N6",\
                                            "C4'","H4'","O4'","C1'","H1'","N1","H61","H62","O5'","H5'","C8","H8", \
                                            "C6","H6","C5","C7","H71","H72","H73","C4","O4","N3","H3","C2","O2",\
                                                "C3'","H3'","C2'","H2'","H2''","O3'","O6","N2","N4"]

    
def get_atoms(resname):
    
    if resname == "DT" or resname == "TT":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                     "C1'","H1'","N1","C6","H6","C5","C7","H71","H72","H73","C4",\
                     "O4","N3","H3","C2","O2","C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 32
                             
    elif resname == "DA":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                      "C1'","H1'","N9","C8","H8","N7","C5","C6","N6","H61","H62",\
                      "N1","C2","H2","N3","C4","C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 32
    elif resname == "DG":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                     "C1'","H1'","N9","C8","H8","N7","C5","C6","O6","N1","H1","C2",\
                     "N2","H21","H22","N3","C4","C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 33
            
    elif resname == "DC":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'",\
                     "H1'","N1","C6","H6","C5","H5","C4","N4","H41","H42","N3","C2","O2",\
                     "C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 30
        
    elif resname == "DT5":
        atom_list = ["HO5'","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'",\
                     "H1'","N1","C6","H6","C5","C7","H71","H72","H73","C4","O4",\
                     "N3","H3","C2","O2","C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 30
            
    elif resname == "DT3":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                     "C1'","H1'","N1","C6","H6","C5","C7","H71","H72","H73","C4",\
                     "O4","N3","H3","C2","O2","C3'","H3'","C2'","H2'","H2''","O3'","HO3'"]
        natoms = 33
            
    elif resname == "DG5":
        atom_list = ["HO5'","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'",\
                     "N9","C8","H8","N7","C5","C6","O6","N1","H1","C2","N2","H21",\
                     "H22","N3","C4","C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 31
        
    elif resname == "DG3":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                     "C1'","H1'","N9","C8","H8","N7","C5","C6","O6","N1","H1","C2",\
                     "N2","H21","H22","N3","C4","C3'","H3'","C2'","H2'","H2''","O3'","HO3'"]
        natoms = 34
        
    elif resname == "DC5":
        atom_list = ["HO5'","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'",\
                     "N1","C6","H6","C5","H5","C4","N4","H41","H42","N3","C2","O2","C3'",\
                     "H3'","C2'","H2'","H2''","O3'"]
        natoms = 28
        
    elif resname == "DC3":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                     "C1'","H1'","N1","C6","H6","C5","H5","C4","N4","H41","H42",\
                     "N3","C2","O2","C3'","H3'","C2'","H2'","H2''","O3'","HO3'"]
        natoms = 31
        
    elif resname == "DA5":
        atom_list = ["HO5'","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'",\
                     "H1'","N9","C8","H8","N7","C5","C6","N6","H61","H62","N1",\
                     "C2","H2","N3","C4","C3'","H3'","C2'","H2'","H2''","O3'"]
        natoms = 30
        
    elif resname == "DA3":
        atom_list = ["P","OP1","OP2","O5'","C5'","H5'","H5''","C4'","H4'","O4'",\
                     "C1'","H1'","N9","C8","H8","N7","C5","C6","N6","H61","H62",\
                     "N1","C2","H2","N3","C4","C3'","H3'","C2'","H2'","H2''","O3'","HO3'"]
        natoms = 33
    else:
        print("Error: residue not recognised")
        
    return atom_list, natoms
        

    
    

class trajectory:
    def __init__(self, file_name, circular):
        self.file_name = file_name
        print("Processing " + self.file_name)
        self.circular = circular
        self.file_type = self.get_file_type()
        
    def get_steps_pdb(self):
        lines = open(self.file_name,"r").readlines()
        rev_lines = reversed(lines)
        for line in rev_lines:
            split_line = line.split()
            if split_line[0]=="MODEL":
                steps = int(split_line[1])
                break
        self.steps = steps
        return steps
    
    def get_file_type(self):
        split_name = self.file_name.split(".")
        suffix = split_name[-1]
        if suffix=="pdb":
            file_type = "pdb"
        else:
            file_type = "amber"
        return file_type
    
    def load_sequence(self,seq_file):
        line = open(seq_file,"r").readline()
        seq = list(line)
        seq = [x for x in seq if x.strip()]
        self.seq = seq
        return seq
    
    def load_topology(self,top_file):
        self.strand_list = []
        '''
        Loads in amber topology, returning the list of residues

        '''
        lines = open(top_file,"r").readlines()
        res_bool = False
        for line in lines:
            if line.strip() == "%FLAG RESIDUE_LABEL":
                residues = []
                res_bool = True
            elif line.strip() == "%FLAG RESIDUE_POINTER":
                res_bool = False
            if res_bool:
                new_residues = [res for res in line.split()]
                residues += new_residues
                
        residues.pop(0)
        residues.pop(0)
        residues.pop(0)
                
        self.residues = residues
        self.nres = len(self.residues)
        
        sub_res = []
        self.direction = residues[0][2]
        if self.direction=="5":
            self.direction += "3"
        else:
            self.direction += "5"
            
        for res in residues:
            if len(res) == 3 and res[2]==self.direction[1]:
                sub_res.append(res)
                self.strand_list.append(strand(self.circular, sub_res))
                sub_res = []
            else:    
                sub_res.append(res)
                
        self.nstrand = len(self.strand_list)
        
        print("Amber topolgy loaded.\n")
        print( str(self.nres) + " residues and " +str(self.nstrand) + " strands detected")
        for i, strandy in enumerate(self.strand_list):
            print("Strand " + str(i+1) + " has " + str(len(strandy.res)) + " residues.")
            
                
        self.residues = residues
            
        return residues


        
       
 
    def process_configurations(self, list_of_tups, max_steps=100000):
        """
        

        Parameters
        ----------
        list_of_tups : list
            list contains tuples (name, func) where func calculates a quantity given a quantity name
             

        Returns
        -------
        None.

        """
        
        self.quant_dict = {}
        for tup in list_of_tups:
            self.quant_dict[tup[0]] = []
            
            

            
        if self.file_type == "pdb":
            strandname = True # sometimes strand labelling is extra column which we shall ignore
            with open(self.file_name,"r") as file:
                for line in file:
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
                        if strandname:
                            split_line.pop(4)
                        atom = split_line[2]
                        resname = split_line[3]
                        bi = int(split_line[4])
                        try:
                            x,y,z = float(split_line[5]),float(split_line[6]),float(split_line[7])
                        except:
                            subline = split_line[5]
                            sep = subline.index("-") 
                            sep = [i for i, x in enumerate(subline) if x=="-"]
                            val = sep[-1]
                            x,y = float(subline[:val]), float(subline[val:])
                            z = float(split_line[6])
                        
                        if bi==bimatch:   # if base index matches, add atom to the base dictionary
                            base_dict["resname"] = resname
                            base_dict["res_number"] = bi
                            keys = base_dict.keys()
                            if atom in keys:
                                base_dict[atom] = np.array([x,y,z])
          
                        else:                 
                            base1 = base(base_dict)
                            strand1.add_base(base1) # add base to the configuration class
                            base_dict = {}          # reset the base dictionary
                            base_dict = dict.fromkeys(list1)
                            base_dict["resname"] = resname
                            base_dict["res_number"] = bi
                            keys = base_dict.keys()
                            if atom in keys:
                                base_dict[atom] = np.array([x,y,z])
                            bimatch = bi # increments the base index by one
                    elif word=="TER":
                        strand1.add_base(base1) # add base to the configuration class
                        config.add_strand(strand1)             
                        strand1 = strand(self.circular)
                        bimatch += 1 # this just means next strand doesn't pick up last residue of last strand
                    elif word=="ENDMDL":
                        0
                    else:
                        pass
        else:
            step = 1
            skip = False
            config = configuration(step)
            config.load_strand_list(self.strand_list)
            with open(self.file_name,"r") as file:
                long_line = []
                strand_i = 0
                res_i = 0
                res = config.strand_list[0].res[0]
                base_dict = dict.fromkeys(list1)
                atoms, natoms = get_atoms(res)
                base_dict["resname"] = res
                first_line = file.readline()
                
                if len(re.findall("-?\d*\.\d{3}", first_line)) > 0: # skip first line if there are no coordinates (eg Generated by Cpptraj)
                    file.seek(0)
                
                for line in file:
                    if skip == True:
                        skip = False
                    else:
                        split_line = line.split()
                        try:
                            split_line = [float(val) for val in split_line]
                        except ValueError:
                            split_line = re.findall("-?\d*\.\d{3}", line)
                            split_line = [float(val) for val in split_line]
                        long_line += split_line
                        if len(long_line) >= 3 * natoms:                             
                            j = 0
                            for atom in atoms:
                                x,y,z = long_line[j],long_line[j+1],long_line[j+2]
                                j += 3
                                base_dict[atom] = np.array([x,y,z])
                                
                            base1 = base(base_dict)
                            config.strand_list[strand_i].add_base(base1)
                            res_i += 1
                            remainder = len(long_line) - 3 * natoms
                            long_line = split_line[-remainder:]
                            
                            if remainder == 0:
                                long_line = []

                            if res_i == config.strand_list[strand_i].nres:
                                res_i = 0
                                strand_i += 1
                                if strand_i == config.nstrands:
                                    for tup in list_of_tups:
                                        self.quant_dict[tup[0]].append(tup[1](config))
                                    
                                    self.last_conf = config
                                    print("Step:", step, end='\r', flush=True)
                                    # print("Step:", step)
                                    strand_i = 0                             
                                    step += 1                 
                                    config = configuration(step, direction=self.direction)
                                    config.load_strand_list(self.strand_list)
                                    skip = True
                                    if step > max_steps:
                                        break
                         
                            base_dict = dict.fromkeys(list1)
                            res = config.strand_list[strand_i].res[res_i]
                            atoms, natoms = get_atoms(res)
                            base_dict["resname"] = res
                                                 

class strand:
    def __init__(self,circular,res):
        self.base_list = [] # list of base classes
        self.res = res #residue list
        self.nres = len(self.res)
        self.circular = circular
        
    def add_base(self, base):
        self.base_list.append(base)
        
    def reverse_strand(self):
        self.base_list.reverse()
        self.res.reverse()

        
    def get_bb_list(self):
        self.bb_list = [base.bb for base in self.base_list]
        return self.bb_list
            
    def get_atom_list(self,atom):
        """
        gets the list of one particular atom as a function of base-pair index
        
        """
        atom_list = [base.get_atom(atom) for base in self.base_list]    
        return atom_list
            
    def get_base_centres(self):
        self.base_centres = np.array([base.get_base_centre() for base in self.base_list])
        return self.base_centres

        
        

class base:
    def __init__(self, base_dict):
        self.base_dict = base_dict
        self.resname = self.base_dict["resname"]
        self.bb = self.get_bb()

        
    def get_base_centre(self):
        n1,c6,c5,c4,n3,c2 = self.base_dict["N1"],self.base_dict["C6"],self.base_dict["C5"],self.base_dict["C4"],self.base_dict["N3"],self.base_dict["C2"]
        self.base_centre = (n1+c6+c5+c4+n3+c2)/6
        return self.base_centre
        
    def get_atom(self,atom):
        return self.base_dict[atom]
    
    def get_bb(self):
        if len(list(self.resname)) > 2:
            if self.resname[2] == "3":
                self.bb = self.base_dict["HO3'"]
            else:
                self.bb = self.base_dict["HO5'"]
                
        else:
            self.bb = self.base_dict["P"]
        return self.bb
    
    def natoms(self):
        natoms = 0
        for val in self.base_dict.values():
            if isinstance(val, np.ndarray):
                natoms += 1
        return natoms
            



class configuration:
    '''
    decribes a single time-step configuration
    '''
    def __init__(self,step,direction="35"):
        self.step = step
        self.nstrands = 0
        self.direction = direction
        
    def add_strands(self):
        self.strand_list.append(strand)
        self.nstrands +=1
        
    def load_strand_list(self, strand_list):
        self.strand_list = copy.deepcopy(strand_list)
        self.nstrands = len(self.strand_list)
        
    def clear(self):
        self.strand_list = []
      


    def get_duplex_centres(self, helical=True,alpha=0):
        '''
        Returns
        -------
        centres, a array of the midpoints of a duplex

        '''
        strand1 = self.strand_list[0]
        strand2 = self.strand_list[1]
        strand1.get_base_centres()
        strand2.get_base_centres()
        centres = np.array([(strand1.base_centres[i] + strand2.base_centres[-i-1]) /2 for i in range(strand1.nres)])
        self.centres = centres
        return self.centres

        
