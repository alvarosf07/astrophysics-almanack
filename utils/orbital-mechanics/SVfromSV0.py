# -*- coding: utf-8 -*-

"""
Programa: ejercicio1.py (SVfromSV0rel) 
Descripción: Cálculo de Vector de Estado (State Vector, SV (t)) en un Tiempo (t) a partir del vector de estado (SV0) en un tiempo incial (t0)
Autor: Álvaro Sánchez Fernández
Fecha: 05/03/2019
"""

"""
def f(x):
    return x+1

y=f(4.597)
print("nNumero con dos decimales: %, 2f, %3f" %(y,y))
"""

import numpy as np

#print ("\n\nCOMIENZO DEL PROGRAMA\n\n")

m1=5.974e24
m2=1000
G=6.6742e-11
mu=G*(m1+m2)
Rt=6378000

r0=np.array((9700000,0,4200000))
v0=np.array((0,6000,3000))
tiempo=4*60*60
delta_t=1


# 1) Definición de todas los vectores de estado en el instante inicial
f0=[r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]]
mod_r0=((f0[0]**2)+(f0[1]**2)+(f0[2]**2))**(1/2)
mod_v0=((f0[3]**2)+(f0[4]**2)+(f0[5]**2))**(1/2)
#print("f0 = ",f0)

# 2) Definición de todas los vectores de estado en el instante inicial
df0=[f0[3],f0[4],f0[5],-(mu*f0[0])/(mod_r0**3),-(mu*f0[1])/(mod_r0**3),-(mu*f0[2])/(mod_r0**3)]
#print("df0 = ",df0)

# 3) Algoritmo Integración de Runge Kutta para la obtención del vector de estado en el instante ti

# Función Integrar
def integra (fi0,dfi0,h):
    #print("\nFunción integra")
    qi=np.array(dfi0)
    #print("qi = ", qi)
    fi=np.array(fi0+(qi*h))
    return fi
    
fi0=np.array(f0)
dfi0=np.array(df0)
h=delta_t
    
h1=mod_r0
v1=mod_v0
h2=mod_r0
v2=mod_v0
    
for i in range (4*60*60+1):
    #print("\n\nEntrada en For", i,"\n")
    #print("fi0 = ", fi0)
    #print("dfi0 = ", dfi0)
    fi=integra (fi0,dfi0,h)
    #print("fi = ", fi)
    mod_ri=((fi[0]**2)+(fi[1]**2)+(fi[2]**2))**(1/2)
    mod_vi=((fi[3]**2)+(fi[4]**2)+(fi[5]**2))**(1/2)
    dfi=np.array((fi[3],fi[4],fi[5],-(mu*fi[0])/(mod_ri**3),-(mu*fi[1])/(mod_ri**3),-(mu*fi[2])/(mod_ri**3)))
    #print ("dfi = ", dfi)
    
    #print ("\nModulo R mínima: ", h1)
    #print ("Modulo R máxima: ", h2)
    #print ("Modulo Ri: ", mod_ri)
    #print ("Modulo vi: ", mod_vi) 
    
    if mod_ri < h1:
        h1=mod_ri
        v1=mod_vi
        
    if mod_ri > h2:
        h2=mod_ri
        v2=mod_vi
        
    fi0=fi
    dfi0=dfi

#h1= round ((h1-Rt)/1000,1)
#v1= round (v1/1000,1)
#h2= round ((h2-Rt)/1000,1)
#v2= round (v2/1000,1)

print("h1(km)={:.1f}".format((h1-Rt)/1000))
print("v1(km/s)={:.1f}".format(v1/1000)) 
print("h2(km)={:.1f}".format((h2-Rt)/1000)) 
print("v2(km/s)={:.1f}".format(v2/1000)) 

    
    
        
    
        
        
        
        
        

    
    





