#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 23:23:07 2019
Title: Practica2
Description: Atmospheric Model definition. Simulation of a sky-diver fall in ISA Atmospheric Conditions
@author: alvarosanchezfernandez
"""

import numpy as np # Librería numérica
import matplotlib.pyplot as plt # Librería para representación

def integra(fi0,dfi0,h): # Código para calcular el resultado de la integración
    qi=dfi0
    fi=fi0
    fi.extend ([fi0[-4]+(qi[-4]*h),fi0[-3]+(qi[-3]*h),fi0[-2]+(qi[-2]*h),fi0[-1]+(qi[-1]*h)])
    #fi=fi0+(dfi0*h)
    return fi
    
def main(): # Código de la función principal
    R=6370
    G=6.67259e-20
    M=5.976e24
    MU=G*M
    
    #r0=float(input("Introduzca la distancia inicial r0 (km):"))
    #phi0=float(input("Introduzca el ángulo inicial phi0 (º):"))
    #vr0=float(input("Introduzca la velocidad normal inicial Vr0 (m/s):"))
    #vphi0=float(input("Introduzca la velocidad tangencial inicial Vphi0 (m/s):"))
    #h=float(input("Introduzca el periodo de integración (s):"))
    #t=float(input("Introduzca tiempo total de la simulación (s):"))
    
    r0=4*R
    phi0=0
    vr0=4
    vphi0=-3
    h=1
    t_fin=100*60*60
    
    f0=[r0,phi0,vr0,vphi0]
    #print("f0=",f0)
    df0=[vr0,vphi0/r0,((vphi0**2)/r0)-((MU)/(r0**2)),-(vr0*vphi0)/(r0)]
    #print("df0=",df0)
    
    r=[r0]
    rx=[r[0]*np.cos(phi0)]
    ry=[r[0]*np.sin(phi0)]
    phi=[phi0]
    t=[0]
    
    fi0=f0
    dfi0=df0
    ti=0
    
    while t_fin>ti:
        fi=integra(fi0,dfi0,h)
        #print ("fi=",fi)
        dfi=[fi[-2],(fi[-1]/fi[-4]),((fi[-1]**2)/fi[-4])-((MU)/(fi[-4]**2)),-(fi[-2]*fi[-1])/(fi[-4])]
        #print ("dfi=",dfi)
        r.extend([fi[-4]])
        rx.extend([fi[-4]*np.cos(fi[-3])])
        ry.extend([fi[-4]*np.sin(fi[-3])])
        phi.extend([fi[-3]*180/np.pi])
        t.extend([t[-1]+h])
        #print ("rx=",rx)
        #print ("ry=",ry)
        fi0=fi
        dfi0=dfi
        ti=ti+h
        if  fi[-4]<R:
            break
        
    tarray=np.array(t)
    th=tarray/3600
    
    print ("\n")
    print ("Valores de las Variables de Estado al final de la simulación:")  
    print (" R final: {:.2f} (km)".format(r[-1]))
    print (" Ángulo phi final: {:.2f} (deg)".format(phi[-1]))
    print (" Velocidad normal final Vr: {:.4f} (km/s)".format(fi[-2]))
    print (" Velocidad tangencial final Vphi: {:.4f} (km/s)".format(fi[-1]))
    print (" Tiempo Total de Simulación:",t[-1],"s")
    print ("\n")
    print ("Gráficos de la Trayectoria:\n")
    
    fig1 = plt.figure() # Crea una figura

    ax1 = fig1.add_subplot(111) # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax1.set_xlabel("t (h)") # Asigna título al eje x del gráfico ax1
    ax1.set_ylabel("R (km)") # Asigna título al eje y del gráfico ax1
    ax1.plot(th,r,"b") # Representa la curva x-y en color azul
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111) # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax2.set_xlabel("Posición en X (km)") # Asigna título al eje x del gráfico ax1
    ax2.set_ylabel("Posiciónn en Y (km)") # Asigna título al eje y del gráfico ax1
    ax2.plot(rx,ry,"r") # Representa la curva x-y en color azul
    
    fig3 = plt.figure() # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax3 = fig3.add_subplot(111)
    ax3.set_xlabel("Phi (deg)") # Asigna título al eje x del gráfico ax1
    ax3.set_ylabel("R (km)") # Asigna título al eje y del gráfico ax1
    ax3.plot(phi,r,"g") # Representa la curva x-y en color azul
    #fig.tight_layout()

    fig1.show()
    fig2.show()
    fig3.show()
    
    
if __name__ == "__main__": main()
