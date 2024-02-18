#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 18:05:27 2019
ejercicio_1.5
@author: alvarosanchezfernandez
"""

import numpy as np

MU=398600
RT=6378

# Se pide por pantalla la introducción del vector de estado de la posición inicial del satélite, así como del incremento de anomalía verdadera
print("Introduzca las componentes del vector de estado inicial (r0=rx0 i + ry0 j km, v0=vx0 i + vy0 j km/s):")
#rx0= float(input("rx0?"))
#ry0= float(input("ry0?"))
#vx0= float(input("vx0?"))
#vy0= float(input("vy0?"))

rx0= 8182.4
ry0= -6865.9
vx0= 0.47572
vy0= 8.8116

#delta_theta= float(input("Valor del incremento de anomalía verdadera (grados)?"))
delta_theta=120

# Algoritmo de Cálculo del Vector de Estado - Método Coeficientes de Lagrange

    # 1/ Cálculo de las magnitudes de r0 y v0
r0=np.array([rx0,ry0,0])
v0=np.array([vx0,vy0,0])

mod_r0=np.sqrt((rx0**2)+(ry0**2))
mod_v0=np.sqrt((vx0**2)+(vy0**2))

    # 2/ Cálculo de la componente radial de v0
vr0=np.dot(r0,v0)/mod_r0

    # 3/ Cálculo del momento angular h
h=mod_r0*np.sqrt(mod_v0**2-vr0**2)

    # 4/ Cálculo de la magnitud del vector de posición
mod_r=(h**2/MU)*(1/(1+((((h**2)/(MU*mod_r0))-1)*np.cos(delta_theta*np.pi/180))-((h*vr0*np.sin(delta_theta*np.pi/180))/MU)))

    # 5/ Cálculo de los coeficientes de Lagrange
f=1-(((MU*mod_r)/(h**2))*(1-np.cos(delta_theta*np.pi/180)))
g=((mod_r*mod_r0)/(h))*np.sin(delta_theta*np.pi/180)
df=(MU/h)*((1-np.cos(delta_theta*np.pi/180))/(np.sin(delta_theta*np.pi/180)))*(((MU*(1-np.cos(delta_theta*np.pi/180)))/(h**2))-(1/mod_r0)-(1/mod_r))
dg=1-((MU*mod_r0)/(h**2))*(1-np.cos(delta_theta*np.pi/180))

    # 6/ Cálculo final del vector de estado
r=np.array([f*r0[0]+g*v0[0],f*r0[1]+g*v0[1],0])
v=np.array([df*r0[0]+dg*v0[0],df*r0[1]+dg*v0[1],0])

# Cálculo de la anomalía verdadera en t0
theta0=np.arctan((((h**2/MU)-(mod_r)+(((mod_r0-(h**2/MU))*np.cos(delta_theta*np.pi/180))/(mod_r0/mod_r)))*(mod_r0/mod_r))/(((h**2/MU)-mod_r0)*np.sin(delta_theta*np.pi/180)))

# Cálculo de la excentricidad
e=((h**2/MU)-mod_r0)/(mod_r0*np.cos(theta0))


# Impresión por pantalla de los resultados del problema
print("\nResultados:\n")
print("Magnitud de r0={:.2f} km".format(mod_r0))
print("Magnitud de v0={:.2f} km/s".format(mod_v0))
print("vr0={:.2f} km/s".format(vr0))
print("h={:.2f} km^2/s".format(h))
print("r={:.2f} km".format(mod_r))
print("r={:.2f} i".format(r[0]),"+{:.2f} j km".format(r[1]))
print("v={:.2f} i".format(v[0]),"+{:.2f} j km/s".format(v[1]))
print("e={:.2f}".format(e))
print("Anomalía verdadera en t0={:.2f} grados".format(360-(theta0*180/np.pi)))