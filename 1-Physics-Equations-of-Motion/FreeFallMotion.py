#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 13:33:38 2019
Título: Simulación del vuelo de un paracaidista
Parámetros de configuración del programa
z0: altitud del salto
z1: altitud a la que se abre el paracaidas
s: superficie de contacto con el aire
m: masa del paracaidista
Cd0: Coeficiente de resistencia con el paracaidas cerrado 
Cd1: coeficiente de resistencia con el paracaidas abierto 
Utiliza unidades del Sistema Internacional
@author: alvarosanchezfernandez
"""

import numpy as np # Librería numérica
import matplotlib.pyplot as plt # Librería para representación

#import practica2_1 as p2_1
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function


def integra (f,df,delta):
    fi=np.array([f[-2]+delta*df[-2],f[-1]+delta*df[-1]])
    return fi

def cota_geopotencial (z):
    h=(6356766*z)/((6356766+z))
    return h
    
def temperatura_aire (z):
    h=cota_geopotencial (z)
    temp=288.15-0.0065*h
    return temp

def presion_aire (z):
    temp=temperatura_aire (z)
    p=101325*np.e**(((9.80665*28.9644)/(8314.32*(-0.0065)))*np.log(288.15/temp))
    return p

def densidad_aire (z):
    temp=temperatura_aire (z)
    p=presion_aire (z)
    d=(p/temp)*(28.9644/8314.32)
    return d

def main(): # Código de la función principal
    z0=3500
    z1=700
    v0=0
    s=1
    m=75
    Cd0=0.6
    Cd1=1.7
    
    G=9.81
    delta=0.1
    rho=function.densidad_aire(z0)
    
    f=np.array([z0,v0])
    print("f=",f)
    df=np.array([v0,((rho*s*Cd0*v0**2)/(2*m))-G])
    print("df=",df)
    z=np.array([f[-2]])
    v=np.array([f[-1]])
    t=np.array([0])
    
    while f[-2]>=0: #Bucle que nos permitirá realizar la integración e ir guardando los valores de las posiciones
        #print("\nBucle\n")
        rho=function.densidad_aire(f[-2])
        if f[-2]>=z1:
            Cd=Cd0
        else:
            Cd=Cd1
        
        fi= integra(f,df,delta)
        #print("fi=",fi)
        dfi=np.array([fi[-1],((rho*s*Cd*fi[-1]**2)/(2*m))-G])
        #print("dfi=",dfi)
        
        f=np.append(f,fi)
        #print("f=",f)
        df=np.append(df,dfi)
        #print("df=",df)
        z=np.append(z,[f[-2]])
        #print("z=",z)
        v=np.append(v,[f[-1]])
        #print("v=",v)
        t=np.append(t,t[-1]+delta)
        #print("t=",t)
        
    # Interpolación de la solución correcta, al encontrarnos con valores de la altura negativos en el último paso del bucle. 
    # Los vectores finales de z,v y t quedarán corregidos, eliminando el ultimo valor correspondiente a una altura negativa, 
    # e introduciendo el valor correspondiente con la altura final de zf=0
    zf=0
    vf=(v[-2])-(v[-2]-v[-1])*((z[-2]-zf)/(z[-2]-z[-1]))
    tf=(t[-2])-(t[-2]-t[-1])*((z[-2]-zf)/(z[-2]-z[-1]))

    z=np.append(z[0:-1],zf)
    v=np.append(v[0:-1],vf)
    t=np.append(t[0:-1],tf)
      
    #Impresión por pantalla de los resultados finales
    print ("\n\n ------  RESULTADOS DEL SALTO EN PARACAÍDAS:  ------ \n")
    print ("\n\nVelocidad de llegada al suelo: {:.2f} m/s".format(f[-1]))
    print ("\nTiempo total: {:.2f} s".format(t[-1]))
    print ("\n\nGráfico de la Trayectoria:\n")
    
    fig1 = plt.figure() # Crea una figura
    ax1 = fig1.add_subplot(211) # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax1.set_xlabel("Tiempo (s)") # Asigna título al eje x del gráfico ax1
    ax1.set_ylabel("Altitud (m)") # Asigna título al eje y del gráfico ax1
    ax1.plot(t,z,"b") # Representa la curva x-y en color azul
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(211) # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax2.set_xlabel("Tiempo (s)") # Asigna título al eje x del gráfico ax1
    ax2.set_ylabel("Velocidad (m/s)") # Asigna título al eje y del gráfico ax1
    ax2.plot(t,v,"r") # Representa la curva x-y en color azul
    
    fig1.show()
    fig2.show()
    
if __name__ == "__main__": 
    #main()
    sys.exit(main())
