#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 17:22:31 2019

@author: alvarosanchezfernandez
"""
import numpy as np

mu = 132712e+06
r1 = 149.6e+06
r2 = 108.2e+06

# A) Calcular delta-v

# Órbita 1
e1 = 0
h1 = np.sqrt(r1*mu*(1+e1*np.cos(0)))
v1 = h1/r1
T1 = (2*np.pi*r1**(3/2))/(mu**(1/2))

# Órbita 3
ra = r1
rp = r2
e3 = (ra-rp)/(ra+rp)
h3 = np.sqrt (rp*mu*(1+e3*np.cos(0)))
v3a = h3/ra
v3b = h3/rp

# Órbita 2
e2 = 0
h2 = np.sqrt(r2*mu*(1+e1*np.cos(0)))
v2 = h2/r2
T2 = (2*np.pi*r2**(3/2))/(mu**(1/2))

delta_v13 = v1 - v3a
delta_v32 = v3b - v2
delta_v = np.abs(delta_v13) + np.abs(delta_v32)

# B) Tiempo de transferencia A-B
a = (r1+r2)/2
T3 = (2*np.pi*a**(3/2))/(mu**(1/2))
t_ab = T3/2
t_ab_days = (T3/2)/(60*60*24)

# C) Anomalía de Venus respecto a la Tierra en el momento de la salida
theta2_f = 180
theta2_0 = theta2_f - (360*(t_ab/T2))

theta1_0 = 0
theta1_f = theta1_0 + (360*(t_ab/T1))

theta21_0 = theta2_0 - theta1_0
theta21_f = theta2_f - theta1_f


print ("Delta_v = ",delta_v)
print ("t = ",t_ab_days)
print ("theta = ",theta21_0+360)


