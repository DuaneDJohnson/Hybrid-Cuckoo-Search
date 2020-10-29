#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 10:49:37 2018

@author: rahulsn
"""
import numpy as np

    
    

def standardize(nests, Ub, Lb):
    mu_scale = 0.5*(Ub[0, 1] + Lb[0, 1])
    std_scale = 0.5*(Ub[0, 1] - Lb[0, 1])
    standard_nests = (nests - mu_scale)/std_scale
    return standard_nests    
def unstandardize(nests, Ub, Lb):
    mu_scale = 0.5*(Ub[0, 1] + Lb[0, 1])
    std_scale = 0.5*(Ub[0, 1] - Lb[0, 1])
    new_nests = nests*std_scale + mu_scale
    return new_nests

def simplebounds(x, Ub, Lb):
    x = np.where(x < Lb, Lb, x)
    x = np.where(x > Ub, Ub, x)
    return x
