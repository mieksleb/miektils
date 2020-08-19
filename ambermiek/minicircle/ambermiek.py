#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:32:52 2020

@author: michaelselby
"""
import sys
import re

file_name = sys.argv[1]
sequence = open(file_name,"r").readlines()
sequence  = str(sequence)
sequence = re.search('([agct]+)',sequence).group()
sequence_length = len(sequence)


sequence = "g" + sequence + "c"  # we add a g and c to the sequence as these are destoyed in the circularisation
name = "dv"+str(sequence_length)
linear_pdb_file_name = "dv"+str(sequence_length)+"linear"
file = open(name+".nab", "w")
file.write('molecule m;'+'\n'+'m = bdna( \"' + str(sequence)+'\" );'+'\n'+ 'putpdb( \"'+str(linear_pdb_file_name) +'.pdb\", m, \"-wwpdb -nocid -tr\");' )
file.close()

print(sequence_length)





