#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:43:33 2019
Funciones:
@author: alvarosanchezfernandez
"""
import numpy as np

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

def velocidad_sonido (z):
    temp=temperatura_aire(z)
    a=(1.4*temp*(8314.32/28.9644))**(1/2)
    return a

def mach (z,v):
    a=velocidad_sonido(z)
    num_mach=v/a
    return num_mach


def norm (V):
    v = (np.dot(V,V))**(1/2)
    return v

def integra (f,df,delta): #siendo f & df arrays
    #fi=np.array([f[-3]+delta*df[-3],f[-2]+delta*df[-2],f[-1]+delta*df[-1]])
    fi=np.array(f+delta*df)
    return fi

def stumpC (z):
    if z > 0:
        C=(1-np.cos(np.sqrt(z)))/z
    elif z < 0:
        C=(np.cosh(np.sqrt(-z))-1)/(-z)
    else:
        C=1/2;
        
    return C

def stumpS (z):
    if z > 0:
        S=(np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3
    elif z < 0:
        S=(np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3
    else:
        S=1/6
        
    return S

def lagrange_coefficients (R0,V0,delta_theta,MU):
    #MU=398600
    
    mod_r0=np.sqrt((R0[0]**2)+(R0[1]**2)+(R0[2]**2))
    mod_v0=np.sqrt((V0[0]**2)+(V0[1]**2)+(V0[2]**2))
    
    # 2/ Cálculo de la componente radial de v0
    vr0=np.dot(R0,V0)/mod_r0
    
    # 3/ Cálculo del momento angular h
    h=mod_r0*np.sqrt(mod_v0**2-vr0**2)
    
    # 4/ Cálculo de la magnitud del vector de posición
    mod_r=(h**2/MU)*(1/(1+((((h**2)/(MU*mod_r0))-1)*np.cos(delta_theta*np.pi/180))-((h*vr0*np.sin(delta_theta*np.pi/180))/MU)))
    
    # 5/ Cálculo de los coeficientes de Lagrange
    f=1-(((MU*mod_r)/(h**2))*(1-np.cos(delta_theta*np.pi/180)))
    g=((mod_r*mod_r0)/(h))*np.sin(delta_theta*np.pi/180)
    df=(MU/h)*((1-np.cos(delta_theta*np.pi/180))/(np.sin(delta_theta*np.pi/180)))*(((MU*(1-np.cos(delta_theta*np.pi/180)))/(h**2))-(1/mod_r0)-(1/mod_r))
    dg=1-(((MU*mod_r0)/(h**2))*(1-np.cos(delta_theta*np.pi/180)))
    
    # 6/ Cálculo de la anomalía verdadera en t0
    theta0=np.arctan((((h**2/MU)-(mod_r)+(((mod_r0-(h**2/MU))*np.cos(delta_theta*np.pi/180))/(mod_r0/mod_r)))*(mod_r0/mod_r))/(((h**2/MU)-mod_r0)*np.sin(delta_theta*np.pi/180)))

    # 7/ Cálculo de la excentricidad
    e=((h**2/MU)-mod_r0)/(mod_r0*np.cos(theta0))
    
    lagrange_coefficients=np.array([f,g,df,dg,theta0,e])
    
    return lagrange_coefficients

def Kepler_Solver_Ellipse (e,Me):
    # Error Tolerance
    error = 1e-06
    
    # 1/ Select the Starting Value for E
    if Me < np.pi:
        E = Me + e/2
    else:
        E = Me - e/2
    #print ("E0 = ",E)
    
    # 2/ Iteration
    ratio = 1
    while np.abs(ratio) > error:
        ratio = (E-e*np.sin(E)-Me)/(1-e*np.cos(E))
        E = E - ratio
        #print ("Eit = ",E)
    
    return E

def Kepler_Solver_Hyperbola (e,Mh):
    # Error Tolerance
    error = 1e-06
    
    # 1/ Select the Starting Value for F
    F = Mh
    
    # 2/ Iteration
    ratio = 1
    while np.abs(ratio) > error:
        ratio = (e*np.sinh(F)-F-Mh)/(e*np.cosh(F)-1)
        F = F - ratio
        #print ("Fit = ",F)
    
    return F



def sv_from_sv0 (R0,V0,delta_theta,MU): #siendo R0 y V0 vectores de la forma: R0=np.array([rx0,ry0,0]) && V0=np.array([vx0,vy0,0]) || delta_theta: degrees || MU_Earth=398600
    """
    MU=398600
    
    mod_r0=np.sqrt((R0[0]**2)+(R0[1]**2)+(R0[2]**2))
    mod_v0=np.sqrt((V0[0]**2)+(V0[1]**2)+(V0[2]**2))

    # 2/ Cálculo de la componente radial de v0
    vr0=np.dot(R0,V0)/mod_r0

    # 3/ Cálculo del momento angular h
    h=mod_r0*np.sqrt(mod_v0**2-vr0**2)

    # 4/ Cálculo de la magnitud del vector de posición
    mod_r=(h**2/MU)*(1/(1+((((h**2)/(MU*mod_r0))-1)*np.cos(delta_theta*np.pi/180))-((h*vr0*np.sin(delta_theta*np.pi/180))/MU)))

    # 5/ Cálculo de los coeficientes de Lagrange
    f=1-(((MU*mod_r)/(h**2))*(1-np.cos(delta_theta*np.pi/180)))
    g=((mod_r*mod_r0)/(h))*np.sin(delta_theta*np.pi/180)
    df=(MU/h)*((1-np.cos(delta_theta*np.pi/180))/(np.sin(delta_theta*np.pi/180)))*(((MU*(1-np.cos(delta_theta*np.pi/180)))/(h**2))-(1/mod_r0)-(1/mod_r))
    dg=1-((MU*mod_r0)/(h**2))*(1-np.cos(delta_theta*np.pi/180))
    
    # Cálculo de la anomalía verdadera en t0
    theta0=np.arctan((((h**2/MU)-(mod_r)+(((mod_r0-(h**2/MU))*np.cos(delta_theta*np.pi/180))/(mod_r0/mod_r)))*(mod_r0/mod_r))/(((h**2/MU)-mod_r0)*np.sin(delta_theta*np.pi/180)))

    # Cálculo de la excentricidad
    e=((h**2/MU)-mod_r0)/(mod_r0*np.cos(theta0))
    """
    # Cálculo de los coeficientes de Lagrante, theta0 y e
    f = lagrange_coefficients (R0,V0,delta_theta,MU)[0]
    g = lagrange_coefficients (R0,V0,delta_theta,MU)[1]
    df = lagrange_coefficients (R0,V0,delta_theta,MU)[2]
    dg = lagrange_coefficients (R0,V0,delta_theta,MU)[3]
    
    theta0 = lagrange_coefficients (R0,V0,delta_theta,MU)[4]
    e = lagrange_coefficients (R0,V0,delta_theta,MU)[5]

    # Cálculo final del vector de estado
    R = np.array([f*R0[0]+g*V0[0],f*R0[1]+g*V0[1],0])
    V = np.array([df*R0[0]+dg*V0[0],df*R0[1]+dg*V0[1],0])
    
    SV = np.array([R[0],R[1],R[2],V[0],V[1],V[2],theta0,e])
    
    return SV


def t_from_theta (e,h,theta,mu):
    
    if e==0:
        t = ((h**3)*theta)/(mu**2)
        
    if 0<e<1:
        E = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(theta*np.pi/360))
        print ("E = ",E)
        Me = E - e*np.sin(E)
        print ("Me = ",Me)
        t = (Me/(2*np.pi))*(((2*np.pi)/(mu**2))*(h/np.sqrt(1-e**2))**3)
        print ("t = ",t)
        
    if e==1:
        t = (h**3/mu**2)*(((np.tan(theta*np.pi/360))/2)+((np.tan(theta*np.pi/360)**3)/6))
        
    if e>1:
        F = 2*np.arctanh(np.sqrt((e-1)/(e+1))*np.tan(theta*np.pi/360))
        #print ("F = ",F)
        Mh = e*np.sinh(F) - F
        #print ("Mh = ",Mh)
        t = (Mh*h**3)/((mu**2)*(((e**2)-1)**(3/2)))
    
    return t


def theta_from_t (e,h,t,mu):
    
    if e==0:
        theta = (t*mu**2)/(h**3)
        
    if 0<e<1:
        Me = (t*2*np.pi)/(((2*np.pi)/(mu**2))*(h/np.sqrt(1-e**2))**3)
        #print ("Me = ",Me)
        E = Kepler_Solver_Ellipse (e,Me)
        theta = 360*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))/np.pi

    if e==1:
        Mp = (t*(mu**2))/(h**3)
        theta = 360*np.arctan(((3*Mp+np.sqrt(((3*Mp)**2)+1))**(1/3))-((3*Mp+np.sqrt(((3*Mp)**2)+1))**(-1/3)))/np.pi
        
    if e>1:
        Mh = ((mu**2)*(((e**2)-1)**(3/2))*t)/(h**3)
        #print ("Mh = ",Mh)
        F = Kepler_Solver_Hyperbola (e,Mh)
        theta = 360*np.arctan(np.sqrt((e+1)/(e-1))*np.tanh(F/2))/np.pi
        
    if theta < 0:
        theta = theta + 360
    
    return theta


def OE_from_SV (R,V,mu):
    
    r = norm (R)
    v = norm (V)
    
    vr = np.dot(R,V)/r
    #print ("vr = ",vr)
    
    H = np.cross (R,V)
    h = norm (H)
    #print ("H = ",H)
    #print ("h = ",h)
    
    i = np.arccos(H[2]/h) * 180/np.pi
    
    K = np.array ([0,0,1])
    N = np.cross (K,H)
    n = norm (N)
    
    if N [1] >=0:
        RA = np.arccos (N[0]/n) * 180/np.pi
    else:
        RA = 360 - (np.arccos (N[0]/n) * 180/np.pi)
        
    E = ((v**2-(mu/r))*R - (r*vr)*V)/mu
    e = norm (E)
    print ("E = ",E)
    print ("e = ",e)
    
    if E [2] >=0:
        w = np.arccos (np.dot(N,E)/(n*e)) * 180/np.pi
    else:
        w = 360 - (np.arccos (np.dot(N,E)/(n*e)) * 180/np.pi)
    
    if vr >=0:
        TA = np.arccos (np.dot(E,R)/(e*r)) * 180/np.pi
    else:
        TA = 360 - (np.arccos (np.dot(E,R)/(e*r)) * 180/np.pi)
    
    OE = np.array([h,e,i,RA,w,TA])
    
    return OE


def SV_from_OE (h,e,i,RA,w,TA,mu):
    
    rx_PER = ((h**2)/(mu*(1+e*np.cos(TA*np.pi/180)))) * np.cos(TA*np.pi/180)
    ry_PER = ((h**2)/(mu*(1+e*np.cos(TA*np.pi/180)))) * np.sin(TA*np.pi/180)
    rz_PER = ((h**2)/(mu*(1+e*np.cos(TA*np.pi/180)))) * 0
    R_PER = np.array ([rx_PER,ry_PER,rz_PER])
    print ("R_PER = ",R_PER)
    
    vx_PER = (mu/h) * (-np.sin(TA*np.pi/180))
    vy_PER = (mu/h) * (e + np.cos(TA*np.pi/180))
    vz_PER = (mu/h) * 0
    V_PER = np.array ([vx_PER,vy_PER,vz_PER])
    print ("V_PER = ",V_PER)
    
    rx_ECI = PERIFOCAL_to_ECI (rx_PER,ry_PER,rz_PER,i,RA,w) [0]
    ry_ECI = PERIFOCAL_to_ECI (rx_PER,ry_PER,rz_PER,i,RA,w) [1]
    rz_ECI = PERIFOCAL_to_ECI (rx_PER,ry_PER,rz_PER,i,RA,w) [2]
    R_ECI = np.array([rx_ECI,ry_ECI,rz_ECI])
    print ("R_ECI = ",R_ECI)
    
    vx_ECI = PERIFOCAL_to_ECI (vx_PER,vy_PER,vz_PER,i,RA,w) [0]
    vy_ECI = PERIFOCAL_to_ECI (vx_PER,vy_PER,vz_PER,i,RA,w) [1]
    vz_ECI = PERIFOCAL_to_ECI (vx_PER,vy_PER,vz_PER,i,RA,w) [2]
    V_ECI = np.array([vx_ECI,vy_ECI,vz_ECI])
    print ("V_ECI = ",V_ECI)
    
    SV = np.array([R_ECI [0], R_ECI [1], R_ECI [2], V_ECI [0], V_ECI [1], V_ECI [2]])
    
    return SV

def SV_from_Observation (z,ST_deg,lat_deg,r3_TCI,az_deg,elev_deg,r3_TCI_dot,az_deg_dot,elev_deg_dot):
    
    x2_ECI = ECI (z,ST_deg,lat_deg) [0]
    y2_ECI = ECI (z,ST_deg,lat_deg) [1]
    z2_ECI = ECI (z,ST_deg,lat_deg) [2]
    R2_ECI = np.array([x2_ECI,y2_ECI,z2_ECI])
    #print ("x2_ECI = ",x2_ECI) 
    #print ("y2_ECI = ",y2_ECI)
    #print ("z2_ECI = ",z2_ECI)
    
    
    """
    x3_ENU = ENU (r_sigma,az_deg,elev_deg) [0]
    y3_ENU = ENU (r_sigma,az_deg,elev_deg) [1]
    z3_ENU = ENU (r_sigma,az_deg,elev_deg) [2]
    R3_ENU = np.array([x3_ENU,y3_ENU,z3_ENU])
    """
    
    x3_TCI = ENU_to_TCI (r3_TCI,ST_deg,lat_deg,az_deg,elev_deg) [0]
    y3_TCI = ENU_to_TCI (r3_TCI,ST_deg,lat_deg,az_deg,elev_deg) [1]
    z3_TCI = ENU_to_TCI (r3_TCI,ST_deg,lat_deg,az_deg,elev_deg) [2]
    R3_TCI = np.array([x3_TCI,y3_TCI,z3_TCI])
    #print ("x3_TCI = ",x3_TCI) 
    #print ("y3_TCI = ",y3_TCI)
    #print ("z3_TCI = ",z3_TCI)
    
    x1_ECI = x2_ECI + x3_TCI
    y1_ECI = y2_ECI + y3_TCI
    z1_ECI = z2_ECI + z3_TCI
    R1_ECI = np.array([x1_ECI,y1_ECI,z1_ECI])
    
    R2_ECI_dot = np.cross ([0,0,72.9217e-06],R2_ECI)
    #print ("R2_ECI_dot = ",R2_ECI_dot)
    dec = np.arcsin(np.cos(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180)+np.sin(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180))
    #print ("dec = ",dec*180/np.pi)
    dec_dot = (-az_deg_dot*(np.pi/180)*np.cos(lat_deg*np.pi/180)*np.sin(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180)+elev_deg_dot*(np.pi/180)*(np.sin(lat_deg*np.pi/180)*np.cos(elev_deg*np.pi/180)-np.cos(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)))/(np.cos(dec))
    print ("dec_dot = ",dec_dot)
    
    if 0 < az_deg < 180 :
        hh = 360 - np.arccos((np.cos(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)-np.sin(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180))/(np.cos(dec)))*180/np.pi
    else:
        hh = np.arccos((np.cos(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)-np.sin(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180))/(np.cos(dec)))*180/np.pi
        
    ra = (ST_deg - hh)*np.pi/180
    #print ("Ra = ",ra*180/np.pi)
    ra_dot = 72.9217e-06 + (az_deg_dot*(np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180)-elev_deg_dot*(np.pi/180)*np.sin(az_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)+dec_dot*np.sin(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180)*np.tan(dec))/(np.cos(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)-np.sin(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180))
    print ("ra_dot = ",ra_dot)
    R3_ECI_dot  = np.array ([-ra_dot*np.sin(ra)*np.cos(dec)-dec_dot*np.cos(ra)*np.sin(dec),+ra_dot*np.cos(ra)*np.cos(dec)-dec_dot*np.sin(ra)*np.sin(dec),dec_dot*np.cos(dec)])
    
    
    V1_ECI = R2_ECI_dot + r3_TCI_dot*R3_TCI + r3_TCI*R3_ECI_dot    
    
    SV = np.array([R1_ECI [0], R1_ECI [1], R1_ECI [2], V1_ECI [0], V1_ECI [1], V1_ECI [2]])
    
    return SV 


def gibbs (R1,R2,R3,mu):
    """
    %
    % This function uses the Gibbs method of orbit determination
    % to compute the velocity corresponding to the second of
    % three supplied position vectors.
    %
    % mu            - gravitational parameter (kmˆ3/sˆ2)
    % R1, R2, R3    - three coplanar geocentric position vectors
    % (km)
    % r1, r2, r3    - the magnitudes of R1, R2 and R3 (km)
    % c12, c23, c31 - three independent cross products among R1, R2 and R3
    % N, D, S       - vectors formed from R1, R2 and R3 during the Gibbs’ procedure
    % tol           - tolerance for determining if R1, R2 and R3 are coplanar
    % ierr          - = 0 if R1, R2, R3 are found to be coplanar = 1 otherwise
    % V2            - the velocity corresponding to R2 (km/s)
    %
    % User M-functions required: none
    % -----------------------------------------------------------
    """
    #MU=398600 #Tierra
    
    tol = 1e-4
    ierr = 0
    
    # 1/ Calculate r1 r2 and r3
    r1 = norm (R1)
    r2 = norm (R2)
    r3 = norm (R3)
    
    # 2/ Prograde or Retrograde
    CROSS_R1_R2 = np.cross (R1,R2)
    CROSS_R2_R3 = np.cross (R2,R3)
    CROSS_R3_R1 = np.cross (R3,R1)
    
    # 3/ Check that R1, R2 and R3 are coplanar; if not set error flag
    if abs(np.dot(R1,CROSS_R2_R3)/r1/norm(CROSS_R2_R3)) > tol:
        ierr = 1
        print("\n\n Error ",ierr,": No Coplanar Vectors")
    
    # 4/ Calculation of N
    N = r1*CROSS_R2_R3 + r2*CROSS_R3_R1 + r3*CROSS_R1_R2
    
    # 5/ Calculation of N
    D = CROSS_R1_R2 + CROSS_R2_R3 + CROSS_R3_R1
    
    # 6/ Calculation of N
    S = R1*(r2 - r3) + R2*(r3 - r1) + R3*(r1 - r2)
    
    # 7/ Calculation of V2
    V2 = np.sqrt(mu/norm(N)/norm(D))*(np.cross(D,R2)/r2 + S)
    
    print ("\nR2 = ",R2)
    print ("\nV2 = ",V2)
    SV2=np.array([R2[0],R2[1],R2[2],V2[0],V2[1],V2[2]])
    
    return SV2
    

def lambert_y (z,r1,r2,A):
    y=r1+r2+((A*(z*stumpS(z)-1))/(np.sqrt(stumpC(z))))
    return y

def lambert_F (z,t,r1,r2,A,mu):
    F = ((lambert_y(z,r1,r2,A)/stumpC(z))**1.5)*stumpS(z)+A*np.sqrt(lambert_y(z,r1,r2,A))-np.sqrt(mu)*t
    return F

def lambert_dF (z,r1,r2,A):
    if z==0:
        dF = (np.sqrt(2)/40)*((lambert_y(0,r1,r2,A))**(1.5))+(A/8)*(np.sqrt(lambert_y(0,r1,r2,A))+A*np.sqrt(1/(2*lambert_y(0,r1,r2,A))))
    else:
        #dF = ((lambert_y(z,r1,r2,A)/stumpC(z))**1.5)*((1/2*z)*(stumpC(z)-((3*stumpS(z))/2*stumpC(z)))+((3*stumpS(z)**2)/(4*stumpC(z))))+(A/8)*(3*((stumpS(z)/(stumpC(z)))*np.sqrt(lambert_y(z,r1,r2,A)))+A*np.sqrt(stumpC(z)/lambert_y(z,r1,r2,A))) #Calcula dF
        dF=(lambert_y(z,r1,r2,A)/stumpC(z))**1.5*(1/2/z*(stumpC(z)-3*stumpS(z)/2/stumpC(z))+3*stumpS(z)**2/4/stumpC(z))+A/8*(3*stumpS(z)/stumpC(z)*np.sqrt(lambert_y(z,r1,r2,A)) + A*np.sqrt(stumpC(z)/lambert_y(z,r1,r2,A)))
    return dF

def lambert (R1,R2,t,orbit_kind,mu):
    """
    %
    % This function solves Lambert’s problem.
    %
    % MU            - gravitational parameter (kmˆ3/sˆ2)
    % R1, R2        - initial and final position vectors (km)
    % r1, r2        - magnitudes of R1 and R2 (km)
    % t             - the time of flight from R1 to R2
    %                 (a constant) (s)
    % V1, V2        - initial and final velocity vectors (km/s)
    % c12           - cross product of R1 into R2
    % theta         - angle between R1 and R2
    % orbit_kind    - 'pro' if the orbit is prograde
    %                 'retro' if the orbit is retrograde
    % A             - a constant given by Equation 5.35
    % z             - alpha*xˆ2, where alpha is the reciprocal of the
    %                 semimajor axis and x is the universal anomaly
    % y(z)          - a function of z given by Equation 5.38
    % F(z,t)        - a function of the variable z and constant t,
    %                 given by Equation 5.40
    % dFdz(z)       - the derivative of F(z,t), given by
    %                 Equation 5.43
    % ratio         - F/dFdz
    % tol           - tolerance on precision of convergence
    % nmax          - maximum number of iterations of Newton’s
    %                 procedure
    % f, g          - Lagrange coefficients
    % gdot          - time derivative of g
    % C(z),S(z)     - Stumpff functions
    % dum           - a dummy variable
    %
    % User M-functions required: stumpC and stumpS
    % -----------------------------------------------------------
    """
    #MU=398600 #Tierra
    
    # 1/ Calculate r1 and r2
    r1=np.sqrt((R1[0]**2)+(R1[1]**2)+(R1[2]**2))
    r2=np.sqrt((R2[0]**2)+(R2[1]**2)+(R2[2]**2))
    
    # 2/ Prograde or Retrograde  and calculation of theta. IMPORTANT: THETA IS ACTUALLY delta_theta (not theta measured from periapsis)
    cross_vector_R1_R2=np.cross(R1,R2)
    cross_R1_R2=np.dot(cross_vector_R1_R2,cross_vector_R1_R2)
    theta=np.arccos(np.dot(R1,R2)/(r1*r2))
    
    if orbit_kind == "pro":
        if cross_R1_R2 < 0:
            theta = 2*np.pi - theta
    elif orbit_kind == "retro":
        if cross_R1_R2 >= 0:
            theta = 2*np.pi - theta
    else:
        orbit_kind = 'pro'
        print('\n ** Prograde trajectory assumed **.\n')
        if cross_R1_R2 < 0:
            theta = 2*np.pi - theta
    """
    r1=273378
    r2=146378
    theta=5*np.pi/180
    """
    
    # 3/ Calculate A
    A=np.sin(theta)*np.sqrt(r1*r2/(1-np.cos(theta)))
    
    # 4/ Calculate Z by iteration
    z=0;
    """
    while lambert_F(z,t,r1,r2,A,mu) < 0:
        z=z+0.1
    """
    
    tolerance=1e-2
    nmax=10
    ratio=1
    n=0
    
    while (np.abs(ratio) > tolerance) & (n <= nmax):
        
        ratio = lambert_F(z,t,r1,r2,A,mu)/lambert_dF(z,r1,r2,A)
        
        print("\nF iteration ",n,"-->",lambert_F(z,t,r1,r2,A,mu))
        print("dF iteration ",n,"-->",lambert_dF(z,r1,r2,A))
        print("ratio iteration ",n,"-->",ratio)
        print("z iteration ",n,"-->",z)
        
        n=n+1
        z=z-ratio
    
    if n >= nmax:
        print('\n\n ** Number of iterations exceeds')
    
    # 5/ Calculate f,g and gdot
    f=1-(lambert_y(z,r1,r2,A)/r1)
    g=A*(lambert_y(z,r1,r2,A)/mu)**(1/2)
    gdot=1-(lambert_y(z,r1,r2,A)/r2)
    
    # 6/ Calculate v1 and v2
    V1=1/g*(R2-f*R1)
    V2=1/g*(gdot*R2-R1)
    print ("\nR1 = ",R1)
    print ("\nV1 = ",V1)
    print ("\nR2 = ",R2)
    print ("\nV2 = ",V2)
    
    #SV1=np.array([R1[0],R1[1],R1[2],V1[0],V1[1],V1[2]])
    SV2=np.array([R2[0],R2[1],R2[2],V2[0],V2[1],V2[2]])
    
    return SV2


def Gauss (z_A,ST_deg_A,lat_deg_A,t_A,RA_top_A,DEC_top_A,z_B,ST_deg_B,lat_deg_B,t_B,RA_top_B,DEC_top_B,z_C,ST_deg_C,lat_deg_C,t_C,RA_top_C,DEC_top_C,mu):
    """
    mu = 398600

    z_A = 1
    ST_deg_A = 44.506
    lat_deg_A = 40
    t_A = 0
    RA_top_A = 43.537
    DEC_top_A = -8.7833
    
    z_B = 1
    ST_deg_B = 45
    lat_deg_B = 40
    t_B = 118.1
    RA_top_B = 54.42
    DEC_top_B = -12.074
    
    z_C = 1
    ST_deg_C = 45.499
    lat_deg_C = 40
    t_C = 237.58
    RA_top_C = 64.318
    DEC_top_C = -15.105
    """
    delta_t12 = t_A - t_B
    delta_t32 = t_C - t_B
    delta_t31 = t_C - t_A
    #print ("delta_12 = ",delta_t12)
    #print ("delta_32 = ",delta_t32)
    #print ("delta_31 = ",delta_t31)
    
    R2_A_ECI = ECI (z_A,ST_deg_A,lat_deg_A)
    R2_B_ECI = ECI (z_B,ST_deg_B,lat_deg_B)
    R2_C_ECI = ECI (z_C,ST_deg_C,lat_deg_C)
    #print ("R2_A_ECI = ",R2_A_ECI) 
    #print ("R2_B_ECI = ",R2_B_ECI)
    #print ("R2_C_ECI = ",R2_C_ECI)
    
    R3_A_TCI = TCI (z_A,RA_top_A,DEC_top_A)
    R3_B_TCI = TCI (z_B,RA_top_B,DEC_top_B)
    R3_C_TCI = TCI (z_C,RA_top_C,DEC_top_C)
    #print ("R3_A_TCI = ",R3_A_TCI) 
    #print ("R3_B_TCI = ",R3_B_TCI)
    #print ("R3_C_TCI = ",R3_C_TCI)
    
    CROSS_BC = np.cross (R3_B_TCI,R3_C_TCI)
    CROSS_AC = np.cross (R3_A_TCI,R3_C_TCI)
    CROSS_AB = np.cross (R3_A_TCI,R3_B_TCI)
    #print ("CROSS_BC = ",CROSS_BC) 
    #print ("CROSS_AC = ",CROSS_AC)
    #print ("CROSS_AB = ",CROSS_AB)
    
    d0 = np.dot(R3_A_TCI,CROSS_BC)
    d11 = np.dot(R2_A_ECI,CROSS_BC)
    d12 = np.dot(R2_A_ECI,CROSS_AC)
    d13 = np.dot(R2_A_ECI,CROSS_AB)
    d21 = np.dot(R2_B_ECI,CROSS_BC)
    d22 = np.dot(R2_B_ECI,CROSS_AC)
    d23 = np.dot(R2_B_ECI,CROSS_AB)
    d31 = np.dot(R2_C_ECI,CROSS_BC)
    d32 = np.dot(R2_C_ECI,CROSS_AC)
    d33 = np.dot(R2_C_ECI,CROSS_AB)
    #print ("d0 = ",d0) 
    #print ("d11 = ",d11)
    #print ("d12 = ",d12)
    #print ("d13 = ",d13)
    #print ("d21 = ",d21)
    #print ("d22 = ",d22)
    #print ("d23 = ",d23)
    #print ("d31 = ",d31)
    #print ("d32 = ",d32)
    #print ("d33 = ",d33)

    
    aa = (((-d12*delta_t32)/(delta_t31))+d22+((d32*delta_t12)/(delta_t31)))/(d0)
    bb = (((d12*delta_t32)*(delta_t32**2-delta_t31**2)/(delta_t31))+((d32*delta_t12)*(delta_t31**2-delta_t12**2)/(delta_t31)))/(6*d0)
    ee = np.dot (R2_B_ECI,R3_B_TCI)
    r_22 = np.dot (R2_B_ECI,R2_B_ECI)
    #print ("A = ",aa) 
    #print ("B = ",bb) 
    #print ("e = ",ee) 
    #print ("r_22 = ",r_22)

    
    a = - ((aa**2)+(2*aa*ee)+(r_22))
    b = -2*mu*bb*(aa+ee)
    c = -(mu**2)*(bb**2)
    #print ("a = ",a) 
    #print ("b = ",b) 
    #print ("c = ",c) 
    
    # Iteration
    tolerance = 1e-03
    #Pol = np.roots([1,0,-4e+07,0,0,-2.36e+19,0,0,-9.339e30]) 
    Pol = np.roots ([1,0,a,0,0,b,0,0,c])
    Sol_real = Pol.real[np.abs(Pol.imag)<1e-05]
    r2_B = Sol_real.real[Sol_real.real>0][0]
    #r2_B = 9000
    error = r2_B
    
    while (np.abs(error) > tolerance):
        r2_Bi = r2_B - ((r2_B**8+a*r2_B**6+b*r2_B**3+c)/(8*r2_B**7+6*a*r2_B**5+3*b*r2_B**2))
        error = r2_Bi - r2_B
        r2_B = r2_Bi
    
    #print ("POL = ",Pol) 
    #print ("r2_B = ",r2_B)    
    
    r3_A = (((6*(d31*(delta_t12/delta_t32) + d21*(delta_t31/delta_t32))*(r2_B**3) + mu*d31*(delta_t31**2-delta_t12**2)*(delta_t12/delta_t32))/(6*(r2_B**3) + mu*(delta_t31**2-delta_t32**2)))-d11)/d0
    r3_B = aa + ((mu*bb)/(r2_B**3))
    r3_C = (((6*(d13*(delta_t32/delta_t12) - d23*(delta_t31/delta_t12))*(r2_B**3) + mu*d13*(delta_t31**2-delta_t32**2)*(delta_t32/delta_t12))/(6*(r2_B**3) + mu*(delta_t31**2-delta_t32**2)))-d33)/d0
    #print ("r3_A = ",r3_A)
    #print ("r3_B = ",r3_B)
    #print ("r3_C = ",r3_C)

    R1_A_ECI = R2_A_ECI + r3_A*R3_A_TCI
    R1_B_ECI = R2_B_ECI + r3_B*R3_B_TCI
    R1_C_ECI = R2_C_ECI + r3_C*R3_C_TCI
    
    # Lagrange Coefficients
    f1 = 1-(mu*delta_t12**2)/(2*r2_B**3)
    g1 = delta_t12-(mu*delta_t12**2)/(6*r2_B**3)
    f3 = 1-(mu*delta_t32**2)/(2*r2_B**3)
    g3 = delta_t32-(mu*delta_t32**2)/(6*r2_B**3)
    
    V1_B_ECI = (-f3*R1_A_ECI + f1*R1_C_ECI)/(f1*g3 - f3*g1)
    
    R2 = R1_B_ECI
    V2 = V1_B_ECI
    
    SV2=np.array([R2[0],R2[1],R2[2],V2[0],V2[1],V2[2]])
    
    return SV2


def JD_from_date (day,month,year,hour,minute,second):
    """
    day = [1 - 30]
    month = [1 - 12]
    year = [1901 - 2099]
    """
    J0=367*year-((7*(year+((month+9)//(12))))//(4))+((275*month)//(9))+day+1721013.5
    UT=hour+(minute/60)+(second/3600)
    JD=J0+(UT/24)
    
    return JD


def ST_deg_from_date_and_longitude (day,month,year,hour,minute,second,long_deg,long_min,long_sec):
    
    # 1/ Calculate the Julian Day J0 of the Given Date
    J0=367*year-((7*(year+((month+9)//(12))))//(4))+((275*month)//(9))+day+1721013.5
    
    # 2/ Time in Julian Centuries between the given day J0 and the reference day J2000
    T0=(J0-2451545)/(36525)
    
    # 3/ Calculation of the Greenwich Sidereal Time at UT
    ST_Green_0 =(100.4606184)+(36000.77004*T0)+(0.000387933*T0**2)-((2.583e-8)*(T0**3))
    
    if ST_Green_0 > 360:
        Ngreater1=ST_Green_0//360
        ST_Green_0=ST_Green_0-Ngreater1*360
    
    # 4/ Calculation of the Greenwich Sidereal Time at any other Universal Time
    UT=hour+(minute/60)+(second/3600)
    ST_Green=ST_Green_0+(360.98564724*(UT/24))
    
    # 5/ Local Sidereal Time of the Site, is the Greenwich Sidereal Time adding the East Longitude of the Site
    East_Longitude=long_deg+(long_min/60)+(long_sec/3600)
    ST_deg=ST_Green+East_Longitude
    if ST_deg > 360:
        Ngreater2=ST_deg//360
        ST_deg=ST_deg-Ngreater2*360
    
    return ST_deg


def ST_hour_from_date_and_longitude (day,month,year,hour,minute,second,long_deg,long_min,long_sec):
    
    # 1/ Calculate the Julian Day J0 of the Given Date
    J0=367*year-((7*(year+((month+9)//(12))))//(4))+((275*month)//(9))+day+1721013.5
    print ("J0 =",J0)
    
    # 2/ Time in Julian Centuries between the given day J0 and the reference day J2000
    T0=(J0-2451545)/(36525)
    print ("T0 =",T0)
    
    # 3/ Calculation of the Greenwich Sidereal Time at UT
    ST_Green_0 =(100.4606184)+(36000.77004*T0)+(0.000387933*T0**2)-((2.583e-8)*(T0**3))
    print ("ST_Green_0 =",ST_Green_0)
    
    if ST_Green_0 > 360:
        Ngreater1=ST_Green_0//360
        ST_Green_0=ST_Green_0-Ngreater1*360
    
    print ("ST_Green_0 =",ST_Green_0)
    
    # 4/ Calculation of the Greenwich Sidereal Time at any other Universal Time
    UT=hour+(minute/60)+(second/3600)
    ST_Green=ST_Green_0+(360.98564724*(UT/24))
    print ("ST_Green =",ST_Green)
    
    # 5/ Local Sidereal Time of the Site, is the Greenwich Sidereal Time adding the East Longitude of the Site
    East_Longitude=long_deg+(long_min/60)+(long_sec/3600)
    ST_deg=ST_Green+East_Longitude
    print ("ST_deg =",ST_deg)
    if ST_deg > 360:
        Ngreater2=ST_deg//360
        ST_deg=ST_deg-Ngreater2*360
    print ("ST_deg =",ST_deg)
    
    # 6/ Convert Degrees To Hours
    ST_hour=(ST_deg*24)/360.98564724
    
    return ST_hour


def PERIFOCAL_to_ECI (x_PER,y_PER,z_PER,i,RA,w):
    i = i * np.pi/180
    RA = RA * np.pi/180
    w = w * np.pi/180
    
    x_ECI = (np.cos(RA)*np.cos(w)-np.sin(RA)*np.sin(w)*np.cos(i)) * x_PER + (-np.cos(RA)*np.sin(w)-np.sin(RA)*np.cos(i)*np.cos(w)) * y_PER + (np.sin(RA)*np.sin(i)) * z_PER 
    y_ECI = (np.sin(RA)*np.cos(w)+np.cos(RA)*np.cos(i)*np.sin(w)) * x_PER + (-np.sin(RA)*np.sin(w)+np.cos(RA)*np.cos(i)*np.cos(w)) * y_PER + (-np.cos(RA)*np.sin(i)) * z_PER
    z_ECI = (np.sin(i)*np.sin(w)) * x_PER + (np.sin(i)*np.cos(w)) * y_PER + (np.cos(i)) * z_PER
    
    return np.array ([x_ECI,y_ECI,z_ECI])
 
def ECI (h,RA_deg,DEC_deg):
    rt_ec = 6378.13655
    #rt_ec = 6378
    f = 1/298.256421867
    #f = 0.003353
    
    x_ECI = (rt_ec/((1-(2*f-f*f)*np.sin(DEC_deg*np.pi/180)**2)**(1/2))+h) * (np.cos(DEC_deg*np.pi/180)) * (np.cos(RA_deg*np.pi/180))
    y_ECI = (rt_ec/((1-(2*f-f*f)*np.sin(DEC_deg*np.pi/180)**2)**(1/2))+h) * (np.cos(DEC_deg*np.pi/180)) * (np.sin(RA_deg*np.pi/180))
    z_ECI = (rt_ec*(1-f)**2/((1-(2*f-f*f)*np.sin(DEC_deg*np.pi/180)**2)**(1/2))+h) * (np.sin(DEC_deg*np.pi/180))
    #z_ECI = 4078.5
    
    R_ECI = np.array ([x_ECI,y_ECI,z_ECI])
    
    return R_ECI


def ECI_to_ECEF (x_ECI,y_ECI,z_ECI,ST_deg):
    x_ECEF = np.cos(ST_deg*np.pi/180)*x_ECI + np.sin(ST_deg*np.pi/180)*y_ECI
    y_ECEF = -np.sin(ST_deg*np.pi/180)*x_ECI + np.cos(ST_deg*np.pi/180)*y_ECI
    z_ECEF = z_ECI
    
    R_ECEF = np.array ([x_ECEF,y_ECEF,z_ECEF])
    
    return R_ECEF


def ECI_to_ENU (x_ECI,y_ECI,z_ECI,ST_deg,lat_deg):
    x_ENU = - np.sin(ST_deg*np.pi/180)*x_ECI + np.cos(ST_deg*np.pi/180)*y_ECI
    y_ENU = - np.sin(lat_deg*np.pi/180)*np.cos(ST_deg*np.pi/180)*x_ECI - np.sin(lat_deg*np.pi/180)*np.sin(ST_deg*np.pi/180)*y_ECI + np.cos(lat_deg*np.pi/180)*z_ECI
    z_ENU = + np.cos(lat_deg*np.pi/180)*np.cos(ST_deg*np.pi/180)*x_ECI + np.cos(lat_deg*np.pi/180)*np.sin(ST_deg*np.pi/180)*y_ECI + np.sin(lat_deg*np.pi/180)*z_ECI
    
    R_ENU = np.array ([x_ENU,y_ENU,z_ENU])
    
    return R_ENU


def ECEF (rt_h,long_deg,lat_deg):
    x_ECEF = (rt_h) * np.cos(lat_deg*np.pi/180) * (np.cos(long_deg*np.pi/180))
    y_ECEF = (rt_h) * np.cos(lat_deg*np.pi/180) * (np.sin(long_deg*np.pi/180))
    z_ECEF = (rt_h) * np.sin(lat_deg*np.pi/180)
    
    R_ECEF = np.array ([x_ECEF,y_ECEF,z_ECEF])
    
    return R_ECEF

"""
def ECEF_to_ENU (rt_h,x_ECEF,y_ECEF,z_ECEF,long_deg,lat_deg,ST_deg):
    x_ENU = rt_h*np.cos(lat_deg*np.pi/180)*np.cos(long_deg*np.pi/180) + ( - np.sin(ST_deg*np.pi/180)*x_ECEF + np.cos(ST_deg*np.pi/180)*y_ECEF )
    y_ENU = rt_h*np.cos(lat_deg*np.pi/180)*np.sin(long_deg*np.pi/180) + ( - np.sin(lat_deg*np.pi/180)*np.cos(ST_deg*np.pi/180)*x_ECEF - np.sin(lat_deg*np.pi/180)*np.sin(ST_deg*np.pi/180)*y_ECEF + np.cos(lat_deg*np.pi/180)*z_ECEF)
    z_ENU = rt_h*np.sin(lat_deg*np.pi/180) + ( + np.cos(lat_deg*np.pi/180)*np.cos(ST_deg*np.pi/180)*x_ECEF + np.cos(lat_deg*np.pi/180)*np.sin(ST_deg*np.pi/180)*y_ECEF + np.sin(lat_deg*np.pi/180)*z_ECEF)
    
    R_ENU = np.array([x_ENU,y_ENU,z_ENU])
    
    return R_ENU
"""

def TCI (h,RA_top,DEC_top):
    x_TCI = h * np.cos(DEC_top*np.pi/180) * (np.cos(RA_top*np.pi/180))
    y_TCI = h * np.cos(DEC_top*np.pi/180) * (np.sin(RA_top*np.pi/180))
    z_TCI = h * np.sin(DEC_top*np.pi/180)
    
    R_TCI = np.array ([x_TCI,y_TCI,z_TCI])
    
    return R_TCI

def ENU (r_rho,az_deg,elev_deg):
    x_ENU = (r_rho) * np.cos(elev_deg*np.pi/180) * (np.sin(az_deg*np.pi/180))
    y_ENU = (r_rho) * np.cos(elev_deg*np.pi/180) * (np.cos(az_deg*np.pi/180))
    z_ENU = (r_rho) * np.sin(elev_deg*np.pi/180)
    
    R_ECEF = np.array ([x_ENU,y_ENU,z_ENU])
    
    return R_ECEF

def ENU_to_TCI (r_rho,ST_deg,lat_deg,az_deg,elev_deg):
    
    DEC_top = np.arcsin(np.cos(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180)+np.sin(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180))*180/np.pi
    #print ("DEC_top = ",DEC_top)
    
    if 0 < az_deg < 180:
        hh_deg = 360 - np.arccos((np.cos(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)-np.sin(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180))/(np.cos(DEC_top*np.pi/180)))*180/np.pi
    else:
        hh_deg = np.arccos((np.cos(lat_deg*np.pi/180)*np.sin(elev_deg*np.pi/180)-np.sin(lat_deg*np.pi/180)*np.cos(az_deg*np.pi/180)*np.cos(elev_deg*np.pi/180))/(np.cos(DEC_top*np.pi/180)))*180/np.pi
    
    #print ("hh = ",hh_deg)
    
    RA_top = (ST_deg - hh_deg)
    #print ("RA_top = ",RA_top)
    
    x_TCI = TCI (r_rho,RA_top,DEC_top) [0]
    y_TCI = TCI (r_rho,RA_top,DEC_top) [1]
    z_TCI = TCI (r_rho,RA_top,DEC_top) [2]
    
    R_TCI = np.array ([x_TCI,y_TCI,z_TCI])
    
    return R_TCI

def ENU_to_ECI (x_ENU,y_ENU,z_ENU,ST_deg,lat_deg):
    x_ECI = - np.sin(ST_deg*np.pi/180)*x_ENU - np.sin(lat_deg*np.pi/180)*np.cos(ST_deg*np.pi/180)*y_ENU + np.cos(lat_deg*np.pi/180)*np.cos(ST_deg*np.pi/180)*z_ENU
    y_ECI = + np.cos(ST_deg*np.pi/180)*x_ENU - np.sin(lat_deg*np.pi/180)*np.sin(ST_deg*np.pi/180)*y_ENU + np.cos(lat_deg*np.pi/180)*np.sin(ST_deg*np.pi/180)*z_ENU
    z_ECI = + np.cos(lat_deg*np.pi/180)*y_ENU + np.sin(lat_deg*np.pi/180)*z_ENU
    
    R_ECI = np.array ([x_ECI,y_ECI,z_ECI])
    
    return R_ECI
    
    
    
