#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:19:15 2019
Ejercicio: 3.2
@author: alvarosanchezfernandez
"""
import numpy as np # Librería numérica
import matplotlib.pyplot as plt # Librería para representación

#import practica2_1 as p2_1
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function

# Ejericio 3_2 - Pregunta 1
Sol1=function.JD_from_date (21,7,1969,0,0,0)
print ("\n\nPregunta 1: %.5f " %Sol1)

# Ejericio 3_2 - Pregunta 2
Sol2=function.JD_from_date (21,7,1969,2,56,15)
print ("\nPregunta 2: %.5f " %Sol2)

# Ejericio 3_2 - Pregunta 3
Sol3=function.JD_from_date (15,5,2019,17,30,0)
print ("\nPregunta 3: %.5f " %Sol3)

# Ejericio 3_2 - Pregunta 4
Sol4=function.JD_from_date (15,5,2019,17,30,0)-function.JD_from_date (21,7,1969,2,56,15)
print ("\nPregunta 4: %.5f " %Sol4)

# Ejericio 3_2 - Pregunta 5
Sol5=function.ST_deg_from_date_and_longitude (15,5,2019,15,30,0,-5,-34,-18.228)
print ("\nPregunta 5: %.5f " %Sol5)

# Ejericio 3_2 - Pregunta 6
Sol6=function.ST_hour_from_date_and_longitude (15,5,2019,15,30,0,-5,-34,-18.228)
print ("\nPregunta 6: %.5f " %Sol6)

# Ejericio 3_2 - Pregunta 7
Sol7=function.ST_deg_from_date_and_longitude (15,5,2019,8,30,0,-122.33,0,0)
print ("\nPregunta 7: %.5f " %Sol7)

