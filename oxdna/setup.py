#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:49:29 2020

@author: michaelselby
"""
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("cython_writhe.pyx")
)