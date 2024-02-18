#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Programa: Práctica 4 - Simulación de las Ecuaciones de Movmiento de Una Aeronave
Created on Sun May 19 17:54:12 2019

@author: alvarosanchezfernandez
"""

import numpy as np # Librería numérica
import csv
import matplotlib.pyplot as plt # Librería para representación
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

#import practica2_1 as p2_1
import sys
sys.path.insert(0,"../funciones.py")
import funciones as function


def integra (f,df,delta): #siendo f & df arrays
    #fi=np.array([f[-3]+delta*df[-3],f[-2]+delta*df[-2],f[-1]+delta*df[-1]])
    fi=np.array(f+delta*df)
    return fi


def interpolar_curva (tabla3d,x,y):
    col_x = [float(i) for i in list(zip(*tabla3d))[0]] #primera columna
    col_y = [float(i) for i in list(zip(*tabla3d))[1]] #segunda columna
    col_z = [float(i) for i in list(zip(*tabla3d))[2]] #tercera columna
    puntos = [list(i) for i in zip(col_x,col_y)]
    z = griddata(puntos,col_z,(x,y),method="linear")
    return z


"""
def Cl_boeing747 (z,Vc):
    M = function.mach(z,Vc)
    if M>1:
        M=1
    Cl_Mach = "cl_mach.csv"
    with open(Cl_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cl_Mach = list(reader)
        tabla_reales_Cl_Mach = [] # Tabla vacía
        tabla_reales_Cl_Mach = [[float(item) for item in fila] for fila in tabla_Cl_Mach] #convertimos los string en tipo float para poder operar
    Cl = interpolar_curva(tabla_reales_Cl_Mach,z,M)
    return Cl


def Cy_delta_r_boeing747 (z,Vc):
    M=function.mach(z,Vc)
    if M>1:
        M=1
    Cy_delta_r="cy_delta_r.csv"
    with open(Cy_delta_r,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cy_delta_r = list(reader)
        tabla_reales_Cy_delta_r =[] # Tabla vacía
        tabla_reales_Cy_delta_r =[[float(item) for item in fila] for fila in tabla_Cy_delta_r]#convertimos los string en tipo float para poder operar
    Cy_delta_r=interpolar_curva(tabla_reales_Cy_delta_r,z,M)
    return Cy_delta_r


def Cy_beta_boeing747 (z,Vc):
    M=function.mach(z,Vc)
    if M>1:
        M=1
    Cy_beta="cy_beta.csv"
    with open(Cy_beta,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cy_beta = list(reader)
        tabla_reales_Cy_beta =[] # Tabla vacía
        tabla_reales_Cy_beta =[[float(item) for item in fila] for fila in tabla_Cy_beta]#convertimos los string en tipo float para poder operar
    Cy_beta=interpolar_curva(tabla_reales_Cy_beta,z,M)
    return Cy_beta


def Cd_Mach_boeing747 (z,Vc):
    M=function.mach(z,Vc)
    if M>1:
        M=1
    Cd_Mach="cd_mach.csv"
    with open(Cd_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Cd_Mach = list(reader)
        tabla_reales_Cd_Mach =[] # Tabla vacía
        tabla_reales_Cd_Mach =[[float(item) for item in fila] for fila in tabla_Cd_Mach]#convertimos los string en tipo float para poder operar
    Cd=interpolar_curva(tabla_reales_Cd_Mach,z,M)
    return Cd

"""

def thrust_Mach_boeing747 (z,Vc):
    M=function.mach(z,Vc)
    #print("Mach ",M) 
    Thrust_Mach="empuje_mach.csv"
    with open(Thrust_Mach,"r") as f:
        reader = csv.reader(f,delimiter=";")
        tabla_Thrust_Mach = list(reader)
        tabla_reales_Thrust_Mach =[] # Tabla vacía
        tabla_reales_Thrust_Mach =[[float(item) for item in fila] for fila in tabla_Thrust_Mach]#convertimos los string en tipo float para poder operar
    T = interpolar_curva(tabla_reales_Thrust_Mach,25000,M)
    return T



def main():
    
# 1/ Parámetros/Constantes Aeronave
    
    # Embergadura - Wingspan
    b = 195.68*0.3048
    
    # Cuerda Alar de Referencia - Chord
    c = 27.31*0.3048
    
    # Posición relativa del centro de gravedad con respecto a su cuerda alar
    cg = 0.25
    
    # Superficie Alar de Referencia
    S = 5500*(0.3048**2)
    
    # Masa
    m = 288756
    g = 9.8
    
    # Momentos de Inercia
    Ix = 18200000*1.35582
    Iy = 33100000*1.35582
    Iz = 49700000*1.35582
    Ixz = 970000*1.35582
    
    # Definición de la iteración 0: Valores Iniciales
    n=0
    
    # Radio Terrestre (m)
    R = 6378000
    
    
# 2/ Definición del Vector Entradas de Control
    # delta_e : deflexión del timón de profundidad ([-1.0,1.0])
    # delta_a : deflexión de alerones ([-1.0,1.0])
    # delta_r : deflexión del timón de dirección ([-1.0,1-0])
    # delta_t : posición de la palanca de gases ([0,1.0]
    delta_e = 0
    delta_a = 0
    delta_r = 0
    delta_t = 0
    #entry_array = np.asarray ([delta_e,delta_a,delta_r,delta_t])


# 3/ Valores Iniciales
    
    # Velocidad Lineal (m/s)
    u = 203
    v = 0
    w = 0
    aerodynamic_velocities_array = np.array ([u,v,w])
    
    u_dot = 0
    v_dot = 0
    w_dot = 0
    
    # Velocidad Angular
    p = 0
    q = 0
    r = 0
    angular_velocities_array = np.array ([p,q,r])
    
    p_dot = 0
    q_dot = 0
    r_dot = 0
    
    # Posición (Latitud [lambda], Longitud [mu]], Cota [z]) (rad) & (meters)
    lambd = 42.600029*np.pi/180
    mu = -5.5703201*np.pi/180
    z = 8000
    position_array = np.array ([lambd,mu,z])
    
    
    # Actitud (Pitch [theta], Roll [phi]], Yaw[psi]) (rad)
    theta = 0
    phi = 0
    psi = 0
    
    # Cuaternión
    e0 = np.cos(psi/2)*np.cos(theta/2)*np.cos(phi/2)+np.sin(psi/2)*np.sin(theta/2)*np.sin(phi/2)
    e1 = np.cos(psi/2)*np.cos(theta/2)*np.sin(phi/2)-np.sin(psi/2)*np.sin(theta/2)*np.cos(phi/2)
    e2 = np.cos(psi/2)*np.sin(theta/2)*np.cos(phi/2)+np.sin(psi/2)*np.cos(theta/2)*np.sin(phi/2)
    e3 = np.sin(psi/2)*np.cos(theta/2)*np.cos(phi/2)-np.cos(psi/2)*np.sin(theta/2)*np.sin(phi/2)
    norm_e = np.sqrt(e0**2+e1**2+e2**2+e3**2)
    e0 = e0/norm_e
    e1 = e1/norm_e
    e2 = e2/norm_e
    e3 = e3/norm_e
    e = np.array([e0,e1,e2,e3])
    etta = 1-(e0**2+e1**2+e2**2+e3**2) 
    
    a11 = e0**2+e1**2-e2**2-e3**2
    a12 = 2*(e1*e2-e0*e3)
    a13 = 2*(e0*e2+e1*e3)
    a21 = 2*(e1*e2+e0*e3)
    a22 = e0**2-e1**2+e2**2-e3**2
    a23 = 2*(e2*e3-e0*e1)
    a31 = 2*(e1*e3-e0*e2)
    a32 = 2*(e2*e3+e0*e1)
    a33 = e0**2-e1**2-e2**2+e3**2
    #DCM = np.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])


# 5/ Definición del Vector Tiempo
    t0 = 0
    t = t0
    tf = 16
    t_step = 0.01
    time_array = np.arange (t0,tf+t_step,t_step)
    #print (time_array)

        
# 4/ Definición del Vector Estado
    state_array = np.asarray ([t,u,v,w,p,q,r,lambd,mu,z,theta,phi,psi,e0,e1,e2,e3])
    
    
# 6/ Definición de la Matriz de Estado - Temporal
    # Crea una matriz de ceros:
    #   Número de filas: número de elementos de tiempo -> len(time_array)
    #   Número de columnas: número de elementos de estado -> len(state_array)
    state_matrix = np.zeros ((len(time_array),len(state_array)))    
    
    
# 7/ Bucle de Cálculo de las Ecuaciones de Movimiento
    
    for n in range(0,len(time_array)):
        
    # 7.0/ Guardamos en la matriz el vector de estado correspondiente a la fila n
        
        state_matrix [n,:]= state_array
        
    # 7.1/ Actualizar los valores de las entradas de actuación sobre la aeronave
        if time_array[n] < 10.0:
            delta_e = 0
            delta_a = 0
            delta_r = 0
            delta_t = 0.0
            #entry_array = np.asarray ([delta_e,delta_a,delta_r,delta_t])

        if 10 < time_array[n] < 20.0:
            delta_e = 0.0
            delta_a = 0.0
            delta_r = 0.0
            delta_t = 0.0
            #entry_array = np.asarray ([delta_e,delta_a,delta_r,delta_t])
            
        if 20.0 < time_array[n] < 30.0:
            delta_e = 0.0
            delta_a = 0.0
            delta_r = 0.0
            delta_t = 0.0
            #entry_array = np.asarray ([delta_e,delta_a,delta_r,delta_t])
        
        #print ("\nIteration ",n," completed: t = {:.1f}".format(t))
        n = n+1
        t = t+t_step
            
    # 7.2/ Cálculo de la velocidad del aire, ángulo de ataque y ángulo de deslizamiento lateral
        rho = function.densidad_aire(z)
        #print ("u = ",u,"v = ",v,"w = ",w)
        Vc = np.sqrt(u**2+v**2+w**2)
        alpha = np.arctan(w/u)
        alpha0 = 8.5*np.pi/180
        alpha_w = alpha + alpha0
        alpha_dot = (u*w_dot-w*u_dot)/(np.sqrt(u**2+w**2))
        beta = np.arctan(v/np.sqrt(u**2+w**2))
        #beta_dot = (v_dot*np.sqrt(u**2+w**2)-v*(u*u_dot+w*w_dot))/(np.sqrt(u**2+w**2)*(u**2+v**2+w**2))
    
    
    # 7.3/ Obtener los valores de los coeficientes aerodinámicos en el estado actual
        Cl = 1.76
        Cy_delta_r = 0.179
        Cy_beta = -1.08
        Cy = Cy_delta_r*delta_r+Cy_beta*beta
        Cd = 0.263
            
        Cm0 = 0
        Cm_alpha = -1.45
        Cm_delta_e = -1.40
        Cm_q = -21.4
        Cm_alpha_dot = -3.3
            
        Cl_beta = -0.281
        Cl_delta_a = 0.0530
        Cl_delta_r = 0
        Cl_p = -0.502
        Cl_r = 0.195
            
        Cn_beta = 0.184
        Cn_delta_a = 0.0083
        Cn_delta_r = -0.112
        Cn_p = -0.222
        Cn_r = -0.36
        
    
    # 7.4/ Cálculo de las Fuerzas Aerodinámicas
        L = (rho*(Vc**2)*S*Cl)/2
        D = (rho*(Vc**2)*S*Cd)/2
        SF = (rho*(Vc**2)*S*Cy)/2
    
    # 7.5/ Cálculo de las Fuerzas y Momentos Motrices
        Ex = delta_t*thrust_Mach_boeing747 (z,Vc)
        #print("3) Thrust",thrust_Mach_boeing747 (z,Vc))
        Ey = 0
        Ez = 0
        
        El = 0
        Em = 0
        En = 0
    
    # 7.6/ Cálculo de las Fuerzas en Ejes Cuerpo
        Fx = L*np.sin(alpha)-D*np.cos(alpha)-m*g*np.sin(theta)+Ex
        Fy = SF+m*g*np.sin(phi)*np.cos(theta)+Ey
        Fz = -L*np.cos(alpha)-D*np.sin(alpha)+m*g*np.cos(theta)*np.cos(phi)+Ez
        
    # 7.7/ Cálculo de las Aceleraciones en Ejes Cuerpo
        u_dot = Fx/m - q*w+r*v
        v_dot = Fy/m - r*u+p*w
        w_dot = Fz/m - p*v+q*u
        aerodynamic_accelerations_array = np.array ([u_dot,v_dot,w_dot])
        
    # 7.8/ Cálculo de las Velocidades en Ejes Cuerpo
        u = integra (aerodynamic_velocities_array,aerodynamic_accelerations_array,t_step) [0]
        v = integra (aerodynamic_velocities_array,aerodynamic_accelerations_array,t_step) [1]
        w = integra (aerodynamic_velocities_array,aerodynamic_accelerations_array,t_step) [2]
        aerodynamic_velocities_array = np.array ([u,v,w])

    # 7.9/ Cálculo de las Velocidades Terrestres
        Vn = u*a11+v*a12+w*a13
        Ve = u*a21+v*a22+w*a23
        Vd = u*a31+v*a32+w*a33
    
    # 7.10/ Cálculo de las Velocidades de Cambio de Posición
        lambd_dot = Vn/(R+z)
        mu_dot = Ve/(np.cos(lambd)*(R+z))
        z_dot = -Vd
        position_change_array = np.array ([lambd_dot,mu_dot,z_dot])

    # 7.11/ Cálculo de las Posiciones
        lambd = integra (position_array,position_change_array,t_step) [0]
        mu = integra (position_array,position_change_array,t_step) [1]
        z = integra (position_array,position_change_array,t_step) [2]
        position_array = np.array ([lambd,mu,z])
        
    # 7.12/ Cálculo de las Velocidades Angulares en Ejes de Estabilidad
        p_stab = p*np.cos(alpha)+r*np.sin(alpha)
        r_stab = r*np.cos(alpha)-p*np.sin(alpha)

    # 7.13/ Cálculo de los Momentos en los Ejes de Estabilidad
        M_stab = (1/2)*(rho*(Vc**2)*S*c*(Cm0+Cm_alpha*alpha_w+Cm_delta_e*delta_e))+((1/4)*rho*Vc*S*(c**2)*(Cm_q*q+Cm_alpha_dot*alpha_dot))
        L_stab = (1/2)*(rho*(Vc**2)*S*b*(Cl_beta*beta+Cl_delta_a*delta_a+Cl_delta_r*delta_r))+((1/4)*rho*Vc*S*(b**2)*(Cl_p*p_stab+Cl_r*r_stab))
        R_stab = (1/2)*(rho*(Vc**2)*S*b*(Cn_beta*beta+Cn_delta_a*delta_a+Cn_delta_r*delta_r))+((1/4)*rho*Vc*S*(b**2)*(Cn_p*p_stab+Cn_r*r_stab))
        #stability_moments=np.array([M_stab, L_stab, R_stab])
        
    # 7.14/ Cálculo de los Momentos en los Ejes Cuerpo
        M = M_stab+L*(cg-0.25)*c*np.cos(alpha)+D*(cg-0.25)*c*np.sin(alpha)+Em
        L = L_stab*np.cos(alpha)-R_stab*np.sin(alpha)+El
        R = R_stab*np.cos(alpha)+L_stab*np.sin(alpha)-SF*(cg-0.25)*c+En
        
    # 7.15/ Cálculo de las Aceleraciones Angulares en Ejes Cuerpo
        p_dot = (L+((Iy-Iz)*q*r)+(Ixz*(r_dot+p*q)))/(Ix)
        q_dot = (M+(Iz-Ix)*p*r+Ixz*(r**2-p**2))/(Iy)
        r_dot = (R+(Ix-Iy)*p*q+Ixz*(p_dot-q*r))/(Iz)
        angular_accelerations_array = np.array ([p_dot,q_dot,r_dot])
        #p_dot0 = p_dot
        
    # 7.16/ Cálculo de las Velocidades Angulares
        p = integra (angular_velocities_array,angular_accelerations_array,t_step) [0]
        q = integra (angular_velocities_array,angular_accelerations_array,t_step) [1]
        r = integra (angular_velocities_array,angular_accelerations_array,t_step) [2]
        angular_velocities_array = np.array ([p,q,r])
        
    # 7.17/ Cálculo Cuaterniones
        e = np.array([e0,e1,e2,e3])
        etta = 1-(e0**2+e1**2+e2**2+e3**2)
        e0_dot = -(1/2)*(e1*p+e2*q+e3*r)+etta*e0
        e1_dot = (1/2)*(e0*p+e2*r-e3*q)+etta*e1
        e2_dot = (1/2)*(e0*q+e3*p-e1*r)+etta*e2
        e3_dot = (1/2)*(e0*r+e1*q-e2*p)+etta*e3
        e_dot = np.array([e0_dot,e1_dot,e2_dot,e3_dot])
        
        e0 = integra (e,e_dot,t_step) [0]
        e1 = integra (e,e_dot,t_step) [1]
        e2 = integra (e,e_dot,t_step) [2]
        e3 = integra (e,e_dot,t_step) [3]
        norm_e = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        
        e0 = e0/norm_e
        e1 = e1/norm_e
        e2 = e2/norm_e
        e3 = e3/norm_e
        #e = np.array([e0,e1,e2,e3])
        #etta = 1-(e0**2+e1**2+e2**2+e3**2)
        
    # 7.17/ Cálculo DCM        
        a11 = e0**2+e1**2-e2**2-e3**2
        a12 = 2*(e1*e2-e0*e3)
        a13 = 2*(e0*e2+e1*e3)
        a21 = 2*(e1*e2+e0*e3)
        a22 = e0**2-e1**2+e2**2-e3**2
        a23 = 2*(e2*e3-e0*e1)
        a31 = 2*(e1*e3-e0*e2)
        a32 = 2*(e2*e3+e0*e1)
        a22 = e0**2-e1**2-e2**2+e3**2
        
        #DCM = np.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
        
        # Actitud (Pitch [theta], Roll [phi]], Yaw[psi]) (rad)
        theta = np.arcsin(-a31)
        phi = np.arctan2(a32,a33)
        psi = np.arctan2(a21,a11)
        
        rat = psi//np.pi
        psi = psi-rat*np.pi
        
    
    # 7.5/ Actualizar el valor del vector de estado
    
        state_array = np.asarray ([t,u,v,w,p,q,r,lambd,mu,z,theta,phi,psi,e0,e1,e2,e3])
    
    
# 8/ Resultados Finales: Gráficos
    
    #print (state_matrix[:,0])
    #print (state_matrix[:,9])
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Altitude (m)")
    ax1.plot(state_matrix [:,0],state_matrix [:,9],"y")
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(3,1,1)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Pitch (rad)")
    ax2.plot(time_array,state_matrix [:,10],"b")
    
    ax3 = fig2.add_subplot(3,1,2)
    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Roll (rad)")
    ax3.plot(time_array,state_matrix [:,11],"r")
    
    ax4 = fig2.add_subplot(3,1,3)
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("Yaw (rad)")
    ax4.plot(time_array,state_matrix [:,12],"g")
    
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.tight_layout()
    
    fig3 = plt.figure ()
    ax5 = fig3.add_subplot (111, projection="3d")
    ax5.plot(state_matrix [:,7],state_matrix [:,8],state_matrix [:,9],label= " Latitude - Longitude - Altitude ")
    ax5.set_xlabel("Latitud (rad)")
    ax5.set_ylabel("Longitud (rad)")
    ax5.set_zlabel("Cota (m)")
    
    plt.tight_layout()
    
    fig1.show()    
    fig2.show()   
    fig3.show()
    
if __name__ == "__main__": main()
