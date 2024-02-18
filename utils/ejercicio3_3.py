#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 13:22:44 2019

r = (40000,-80790,-100000,200)

@author: alvarosanchezfernandez
"""

import numpy as np # Librería numérica
import matplotlib.pyplot as plt # Librería para representación

#import practica2_1 as p2_1
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function
"""
ST_deg = 110
lat_deg = -40
rt = 6378.13655
h = 0
 
x1_ECI = -2032
y1_ECI = 4591
z1_ECI = -4544
"""
ST_deg = 200
lat_deg = -25
rt = 6378
h = 0.3

x1_ECI = 40000
y1_ECI = -80790
z1_ECI = -100000

R1_ECI = np.array ([x1_ECI,y1_ECI,z1_ECI])
print ("R1_ECI = ",R1_ECI)

x2_ECI = function.ECI (h,ST_deg,lat_deg) [0]
y2_ECI = function.ECI (h,ST_deg,lat_deg) [1]
z2_ECI = function.ECI (h,ST_deg,lat_deg) [2]
R2_ECI = np.array ([x2_ECI,y2_ECI,z2_ECI])
print ("R2_ECI = ",R2_ECI)

x3_ECI = x1_ECI - x2_ECI
y3_ECI = y1_ECI - y2_ECI
z3_ECI = z1_ECI - z2_ECI
R3_ECI = np.array ([x3_ECI,y3_ECI,z3_ECI])
print ("R3_ECI = ",R3_ECI)

R3_ENU = function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)
x3_ENU = function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)[0]
y3_ENU = function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)[1]
z3_ENU = function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)[2]
print ("R3_ENU = ",R3_ENU)

r3_ENU = function.norm(R3_ENU)
print ("r3_ENU_module = ",r3_ENU)
x3_ENU_unit = (function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)[0])/r3_ENU
y3_ENU_unit = (function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)[1])/r3_ENU
z3_ENU_unit = (function.ECI_to_ENU (x3_ECI,y3_ECI,z3_ECI,ST_deg,lat_deg)[2])/r3_ENU
R3_ENU_unit = np.array ([x3_ENU_unit,y3_ENU_unit,z3_ENU_unit])
print ("R3_ENU_unit = ",R3_ENU_unit)

elev_rad = np.arcsin(z3_ENU_unit)
elev_deg = elev_rad*180/np.pi
print ("a = ",elev_deg)
az_rad = np.arccos(y3_ENU_unit/np.cos(elev_rad))
az_deg = az_rad*180/np.pi
print ("A = ",az_deg)

# Ejericio 3_3 - Pregunta 1
Sol1 = x2_ECI
print ("\n\nPregunta 1: %.5f " %Sol1)

# Ejericio 3_3 - Pregunta 2
Sol2 = z3_ECI
print ("\nPregunta 2: %.5f " %Sol2)

# Ejericio 3_3 - Pregunta 3
Sol3 = - np.sin(ST_deg*np.pi/180)
print ("\nPregunta 3: %.5f " %Sol3)

# Ejericio 3_3 - Pregunta 4
Sol4 = z3_ENU
print ("\nPregunta 4: %.5f " %Sol4)

# Ejericio 3_3 - Pregunta 5
Sol5 = z3_ENU_unit
print ("\nPregunta 5: %.8f " %Sol5)

# Ejericio 3_3 - Pregunta 6
Sol6 = elev_deg
print ("\nPregunta 6: %.5f " %Sol6)

# Ejericio 3_3 - Pregunta 7
Sol7 = az_deg
print ("\nPregunta 7: %.5f " %Sol7)