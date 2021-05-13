#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:46:59 2019

@author: michaelselby
"""
import sys
from random import choice

N = int(sys.argv[1]) #length of sequence

def sequence(N):
    
       DNA="CIRCULAR DOUBLE "
       for count in range(N):
          DNA+=choice("CGTA")
       return DNA
   
print(sequence(N))


f = open("sequence_"+str(N), "w+")
f.write(sequence(N))
f.close()

