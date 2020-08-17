#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:32:52 2020

@author: michaelselby
"""
import sys 

top_file = sys.argv[1]
top = open(top_file, "r").readlines()


# converts an oxdna topology file into its pure sequence and the length of that sequence, e.g. agct and 4
def top_2_seq(top_file):
    top = open(top_file, "r").readlines()
    n_bases = float(top[0].split()[0])
    n_strands = float(top[0].split()[1])
    sequences = ""
    for line in top:
        base = line.split()[1]
        base = base.lower()
        sequences = str(sequences)+str(base)
    sequences = sequences[1:] 
    sequence_length = int(n_bases/n_strands)
    sequence = sequences[:sequence_length]
    return sequence, sequence_length




# write the linear pdb file
sequence, sequence_length = top_2_seq(top_file)
seq_file = "seq"+str(sequence_length)+".dat"
seq_file = open(seq_file,"w")
seq_file.write(sequence)
seq_file.close()










