#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:48:53 2019
Programa: practica1.py 
Descripción: Cálculo de la geometría de un tiro parabólico
Autor: Álvaro Sánchez Fernández
Fecha: 06/03/2019
"""

import numpy as np # Librería numérica
import matplotlib.pyplot as plt # Librería para representación

def integra(fi0,dfi0,h): # Código para calcular el resultado de la integración
    qi=dfi0
    fi=fi0
    fi.extend ([fi0[-4]+(qi[-4]*h),fi0[-3]+(qi[-3]*h),fi0[-2]+(qi[-2]*h),fi0[-1]+(qi[-1]*h)])
    return fi
    
def main(): # Código de la función principal
    rx0=0
    ry0=float(input("Introduzca la altura del lanzamiento (m): "))
    v0=float(input("Introduzca la velocidad inicial de lanzamiento (m/s): "))
    alfa0=float(input("Introduzca el ángulo inicial de lanzamiento (º): "))
 
    vx0=v0*np.cos(alfa0*np.pi/180)
    vy0=v0*np.sin(alfa0*np.pi/180)
    G=-9.81
    h=0.01
    
    f0=[rx0,ry0,vx0,vy0]
    #print("f0=",f0)
    df0=[vx0,vy0,0,G]
    #print("df0=",df0)
    
    rx=[rx0]
    ry=[ry0]
    t=[0]
    
    fi0=f0
    dfi0=df0
    ryi=ry0
    
    while ryi>=0: #Bucle que nos permitirá realizar la integración e ir guardando los valores de las posiciones
        fi=integra(fi0,dfi0,h)
        dfi=[fi[-2],fi[-1],0,G]
        
        rx.extend([fi[-4]])
        ry.extend([fi[-3]])
        t.extend([t[-1]+h])

        fi0=fi
        dfi0=dfi
        ryi=fi[-3]
        
    if fi[-3]<0: #Mediante este condicional, interpolaremos la solución correcta en el caso de que el paso de integración nos devuelva posiciones del eje "y" negativas
        fi[-4]=fi[-8]-(((fi[-7])/(fi[-7]-fi[-3]))*(fi[-8]-fi[-4]))
        fi[-3]=0
        
    print ("\n\n ------  RESULTADOS DEL LANZAMIENTO:  ------ \n\n")
    print ("\n\nPosición inicial: [",rx0,",",ry0,"] (m)")
    print ("\nPosición final: [ {:.2f} ,".format(fi[-4]),"{:}] (m)".format(fi[-3]))
    print ("\nTiempo total: {:.2f} s".format(t[-1]))
    print ("\n\nGráfico de la Trayectoria:\n")
    
    fig = plt.figure() # Crea una figura

    ax1 = fig.add_subplot(111) # ax1 -> (2 filas, 1 columna, gráfico 1) 
    ax1.set_xlabel("Alcance (m)") # Asigna título al eje x del gráfico ax1
    ax1.set_ylabel("Cota (m)") # Asigna título al eje y del gráfico ax1
    ax1.plot(rx,ry,"b") # Representa la curva x-y en color azul

    fig.show()
    
if __name__ == "__main__": main()
