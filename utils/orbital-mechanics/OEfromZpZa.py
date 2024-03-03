#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:48:53 2019
Programa: ejercicio1_2.py (SVfromSV0rel) 
Descripción: Cálculo de elementos orbitales a partir de los valores de actitud del perigeo y apogeo
Autor: Álvaro Sánchez Fernández
Fecha: 06/03/2019
"""
import numpy as np

MU=398600
RT=6378

zp= float(input("zp?"))
za= float(input("za?"))

rp=zp+RT
ra=za+RT

a=(rp+ra)/2
e=(ra-rp)/(ra+rp) #e=1-(rp/a)
h=(ra*MU*(1-e))**(1/2)

vacp=(MU*(1+e*np.cos(0)))/(h)
vrp=(MU*e*np.sin(0))/(h)

vaca=(MU*(1+e*np.cos(np.pi)))/(h)
vra=(MU*e*np.sin(np.pi))/(h)

T=(2*np.pi*(a**(3/2)))/(((MU)**(1/2))*3600)

print("e={:.2f}".format(e))
print("h={:.2f}km^2/s".format(h))
print("vp={:.2f}km/s".format(vacp))
print("va={:.2f}km/s".format(vaca))
print("a={:.2f}km".format(a))
print("T={:.3f}h".format(T))







