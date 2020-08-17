#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:16:29 2020

@author: michaelselby

This scipt takes in a base number bpi, and converts a pdb to one with those bases replaced with 
a TT residue for an all atom simulation in amber
"""

import sys

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]
bp = sys.argv[3]


new_lines=[]
with open(input_file_name,"r") as input_file: 
    with open(output_file_name,"w") as output_file:
        for lineno, line in enumerate(input_file):
            split_line = line.split()
            if len(split_line)>1:
                if (int(split_line[4]) == int(bp) or int(split_line[4]) == int(bp)+1):
                    split_line[3] = 'TT'
                    type_line = split_line[2]
                    m=len(type_line)
                    output_file.write("{:>4}".format(split_line[0]))
                    output_file.write("{:>7}".format(split_line[1]))
                    if m>2:
                        output_file.write("{:>5}".format(split_line[2]))
                    else:
                        output_file.write("{:<2}".format(""))
                        split_line[2]+"  "
                        output_file.write("{:<3}".format(split_line[2]))
		    output_file.write("{:>3}".format(split_line[3]))
                    output_file.write("{:>7}".format(split_line[4]))
                    output_file.write("{:>12}".format(split_line[5]))
                    output_file.write("{:>8}".format(split_line[6]))
                    output_file.write("{:>8}".format(split_line[7]))
                    output_file.write("{:>6}".format(split_line[8]))
                    output_file.write("{:>6}".format(split_line[9]))
                    output_file.write("\n") 
                else:
                    output_file.write(line)
