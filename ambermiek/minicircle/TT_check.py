#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:45:00 2020

This script checks whether or not the choice of TT position is correct, 
if it is incorrect it will return a list of possible choices

@author: michaelselby
"""

import sys


seq_file = sys.argv[1]
TT = int(sys.argv[2])

seq_string = open(seq_file, "r").readlines()

seq = list(seq_string[0])


TT_pairs = []
for i in range(0,len(seq)):
    if (seq[i]=='t' and seq[i+1]=='t'):
        TT_pairs.append(i+1)

if (seq[TT]=='t' and seq[TT+1]=='t'):
    print(1)
else:
    print(0)
    print("This is not a correct choice of TT pair, please try one of the following:")
    print(TT_pairs)
    