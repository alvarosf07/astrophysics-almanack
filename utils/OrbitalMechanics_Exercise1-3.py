#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 18:50:32 2019

@author: alvarosanchezfernandez
"""

import numpy as np
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function

rt = 6378
mu = 398600

#Resolución del sistema para calcular e y h
r1 = 2596 + rt
theta_1 = 75
theta_1_rad = theta_1*np.pi/180
r2 = 3369 + rt
theta_2 = 130
theta_2_rad = theta_2*np.pi/180

e = (r1-r2)/(r2*np.cos(theta_2_rad)-r1*np.cos(theta_1_rad))
h = np.sqrt((r1*mu)+(r1*mu*e*np.cos(theta_1_rad)))
p = (h**2)/(mu)

print ("e = ",e)
print ("h = ",h)

# ----------------------------------------------

rp = p/(1+e*np.cos(0))
zp = rp - rt
ra = p/(1+e*np.cos(np.pi))
za = ra - rt

a = (ra+rp)/2
T = (2*np.pi*a**(3/2))/(mu**(1/2))

print ("zp = ",zp)
print ("za = ",za)
print ("a = ",a)
print ("T = ",T/3600)

