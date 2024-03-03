#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 19:18:07 2019

@author: alvarosanchezfernandez
"""

import numpy as np
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function

rt = 6378
mu = 398600

r_geo = ((24*3600*np.sqrt(mu))/(2*np.pi))**(2/3)
v_geo = np.sqrt(mu/r_geo)

z_geo = r_geo - rt

print ("z_geo = ",z_geo)
print ("v_geo = ",v_geo)

# ----------------------------------------------

lat = np.arccos(rt/r_geo)
lat_deg = lat*180/np.pi
sup_terr = 4*np.pi*rt**2
sup_vis = 2*np.pi*(rt**2)*(1-np.cos(lat))
sup_percent = (sup_vis/sup_terr)*100

print ("lat = ",lat)
print ("% Sup. Terrestre Vista = ",sup_percent)
