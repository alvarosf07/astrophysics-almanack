#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Programa: Práctica 3 - Cálculo de Fuerzas Aerodinámicas Sobre Una Aeronave
Fecha Creación: Wed Apr 24 18:10:30 2019
Fecha Última Modificación: Wed Apr 24 18:10:30 2019
@author: alvarosanchezfernandez

Variables:
    z --> cota geométrica (ft)
    h --> cota geopotencial (ft)
    v --> velocidad (kt)
    beta --> ángulo desplazamiento (rad)
    delta_r --> deflexión del timón de dirección (rad) definido entre [1.0,+1.0]
    
Para hacer con arrays en vez de con listas:
    
    np.genfromtxt(file,dtype:front,delimiter=";")
    def interpolar (tabla31,x,y)
        r=tabla31[:,0]
        s=tabla31[:,1]
        t=tabla31[:,2]
        
        puntos=np.stack((r,s),axis=-1)

"""

import numpy as np # Librería numérica
import csv
import matplotlib.pyplot as plt # Librería para representación
from scipy.interpolate import griddata

#import practica2_1 as p2_1
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function


def interpolar_curva(tabla3d,x,y):
    col_x=[float(i) for i in list(zip(*tabla3d))[0]] #primera columna
    col_y=[float(i) for i in list(zip(*tabla3d))[1]] #segunda columna
    col_z=[float(i) for i in list(zip(*tabla3d))[2]] #tercera columna
    puntos =[list(i) for i in zip(col_x,col_y)]
    z=griddata(puntos,col_z,(x,y),method="linear")
    return z



def aerodynamic_forces_boeing747 (z,v,beta,delta_r):
    #Constantes
    rho=function.densidad_aire(z*0.3048)
    S=5500*(0.3048**2)
    M=function.mach(z*0.3048,v*0.5145)
    
    #Obtención de los valores necesarios (Cl,alpha0,Cy_delta_r,Cy_beta,Cd) a partir de las gráficas
    #Cl
    Cl_Mach="cl_mach.csv"
    with open(Cl_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cl_Mach = list(reader)
        tabla_reales_Cl_Mach =[] # Tabla vacía
        tabla_reales_Cl_Mach =[[float(item) for item in fila] for fila in tabla_Cl_Mach]#convertimos los string en tipo float para poder operar
    
    Cl=interpolar_curva(tabla_reales_Cl_Mach,z,M)
    
    #Cy_delta_r
    Cy_delta_r="cy_delta_r.csv"
    with open(Cy_delta_r,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cy_delta_r = list(reader)
        tabla_reales_Cy_delta_r =[] # Tabla vacía
        tabla_reales_Cy_delta_r =[[float(item) for item in fila] for fila in tabla_Cy_delta_r]#convertimos los string en tipo float para poder operar
    
    Cy_delta_r=interpolar_curva(tabla_reales_Cy_delta_r,z,M) 
    
    #Cy_beta
    Cy_beta="cy_beta.csv"
    with open(Cy_beta,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cy_beta = list(reader)
        tabla_reales_Cy_beta =[] # Tabla vacía
        tabla_reales_Cy_beta =[[float(item) for item in fila] for fila in tabla_Cy_beta]#convertimos los string en tipo float para poder operar
    
    Cy_beta=interpolar_curva(tabla_reales_Cy_beta,z,M) 
    
    #Cy
    Cy=Cy_delta_r*delta_r+Cy_beta*beta
    
    #Cd
    Cd_Mach="cd_mach.csv"
    with open(Cd_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cd_Mach = list(reader)
        tabla_reales_Cd_Mach =[] # Tabla vacía
        tabla_reales_Cd_Mach =[[float(item) for item in fila] for fila in tabla_Cd_Mach]#convertimos los string en tipo float para poder operar
    
    Cd=interpolar_curva(tabla_reales_Cd_Mach,z,M)    
    
    #Cálculo de las fuerzas aerodinámicas
    lift=(rho*((v*0.5145)**2)*S*Cl)/2
    sideforce=(rho*((v*0.5145)**2)*S*Cy)/2
    drag=(rho*((v*0.5145)**2)*S*Cd)/2
    aerodynamic_forces=np.array([lift, sideforce, drag, Cl])
    
    return aerodynamic_forces

def aoa_boeing747 (z,v):
    #Constantes
    M=function.mach(z*0.3048,v*0.5145)
    #Obtención de los valores necesarios (Cl,alpha0,Cl_alpha0) a partir de las gráficas 
    #Cl
    Cl_Mach="cl_mach.csv"
    with open(Cl_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cl_Mach = list(reader)
        tabla_reales_Cl_Mach =[] # Tabla vacía
        tabla_reales_Cl_Mach =[[float(item) for item in fila] for fila in tabla_Cl_Mach]#convertimos los string en tipo float para poder operar
    
    Cl=interpolar_curva(tabla_reales_Cl_Mach,z,M)
    
    #alpha0
    alpha0_Mach="alpha0_mach.csv"
    with open(alpha0_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_alpha0_Mach = list(reader)
        tabla_reales_alpha0_Mach =[] # Tabla vacía
        tabla_reales_alpha0_Mach =[[float(item) for item in fila] for fila in tabla_alpha0_Mach]#convertimos los string en tipo float para poder operar
    
    alpha0=interpolar_curva(tabla_reales_alpha0_Mach,z,M)    
    
    #Cl_alpha
    Cl_alpha_Mach="cl_alpha.csv"
    with open(Cl_alpha_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cl_alpha_Mach = list(reader)
        tabla_reales_Cl_alpha_Mach =[] # Tabla vacía
        tabla_reales_Cl_alpha_Mach =[[float(item) for item in fila] for fila in tabla_Cl_alpha_Mach]#convertimos los string en tipo float para poder operar
    
    Cl_alpha=interpolar_curva(tabla_reales_Cl_alpha_Mach,z,M)
    
    #Cálculo del ángulo de ataque (alpha)
    alpha=alpha0+(Cl/Cl_alpha)
    
    return alpha

def main():
    z=20000 #feet
    beta=0.1
    delta_r=-0.05
    v0=184 #Knots
    
    lift=np.array([aerodynamic_forces_boeing747 (z,v0,beta,delta_r)[0]])
    sideforce=np.array([aerodynamic_forces_boeing747 (z,v0,beta,delta_r)[1]])
    drag=np.array([aerodynamic_forces_boeing747 (z,v0,beta,delta_r)[2]])
    cl=np.array([aerodynamic_forces_boeing747 (z,v0,beta,delta_r)[3]])
    alpha=np.array([aoa_boeing747 (z,v0)])
    v=np.array([v0])
    
    while v[-1]<=585: #Bucle que nos permitirá realizar la integración e ir guardando los valores de las posiciones
        #print("\nBucle\n")
        #rho=function.densidad_aire(f[-2])
        h=1 #Paso de la velocidad del bucle
        v=np.append(v,[v[-1]+h])
        
        lift_i=aerodynamic_forces_boeing747 (z,v[-1],beta,delta_r)[0]
        sideforce_i=aerodynamic_forces_boeing747 (z,v[-1],beta,delta_r)[1]
        drag_i=aerodynamic_forces_boeing747 (z,v[-1],beta,delta_r)[2]
        cl_i=aerodynamic_forces_boeing747 (z,v[-1],beta,delta_r)[3]
        alpha_i=aoa_boeing747 (z,v[-1])
        
        lift=np.append(lift,lift_i)
        #print("f=",f)
        sideforce=np.append(sideforce,sideforce_i)
        #print("df=",df)
        drag=np.append(drag,drag_i)
        cl=np.append(cl,cl_i)
        #print("z=",z)
        alpha=np.append(alpha,alpha_i)
        #print("v=",v)
        #t=np.append(t,t[-1]+delta)
        #print("t=",t)
        
    fig1 = plt.figure() # Crea una figura
    ax1 = fig1.add_subplot(2,2,1) # ax1 -> (2 filas, 2 columna, gráfico 1)  ax = fig.add_subplot(2, 3, i)
    ax1.set_xlabel("Speed (knots)") # Asigna título al eje x del gráfico ax1
    ax1.set_ylabel("Lift (N)") # Asigna título al eje y del gráfico ax1
    ax1.plot(v,lift,"y") # Representa la curva x-y en color azul
    
    #fig1 = plt.figure() # Crea una figura
    ax2 = fig1.add_subplot(2,2,2) # ax1 -> (2 filas, 2 columna, gráfico 1) 
    ax2.set_xlabel("Speed (knots)") # Asigna título al eje x del gráfico ax1
    ax2.set_ylabel("Side Force") # Asigna título al eje y del gráfico ax1
    ax2.plot(v,sideforce,"b") # Representa la curva x-y en color azul
    
    #fig1 = plt.figure() # Crea una figura
    ax3 = fig1.add_subplot(2,2,3) # ax1 -> (2 filas, 2 columna, gráfico 1) 
    ax3.set_xlabel("Speed (knots)") # Asigna título al eje x del gráfico ax1
    ax3.set_ylabel("Drag") # Asigna título al eje y del gráfico ax1
    ax3.plot(v,drag,"g") # Representa la curva x-y en color azul
    
    #fig1 = plt.figure() # Crea una figura
    ax4 = fig1.add_subplot(2,2,4) # ax1 -> (2 filas, 2 columna, gráfico 1) 
    ax4.set_xlabel("Speed (knots)") # Asigna título al eje x del gráfico ax1
    ax4.set_ylabel("Alpha") # Asigna título al eje y del gráfico ax1
    ax4.plot(v,alpha,"r") # Representa la curva x-y en color azul
    
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.tight_layout()
    
    fig2 = plt.figure()
    ax5 = fig2.add_subplot(211) # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax5.set_xlabel("Alpha") # Asigna título al eje x del gráfico ax1
    ax5.set_ylabel("Cl") # Asigna título al eje y del gráfico ax1
    ax5.plot(alpha,cl,"r") # Representa la curva x-y en color azul
    
    fig1.show()    
    fig2.show()

    
if __name__ == "__main__": main()


    

