# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:22:04 2023

@author: Lucie
"""

import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

#---------------------------------------------Data---------------------------------------------
#Radius (km)
Re = 6378.135
Rl = 1738
Rs = 696340
Rma = 3389.5
Rv = 6051
Rme = 2439
Rj = 69911
Rn = 24622
Ru = 25362 
Rst = 58232

#Gravitationnal constant (km3/s2)
mueE = 3.986005*(10**5)
mueS = 13.2*10**10
mueL = 4903 
mueST = 38*10**6

#Rotation speed ()
Je = 1.08263*10**-3

#Altitude (km)
zGEO = 35786
zminLEO = 180
zmaxLEO = 2000

rGEO = zGEO+Re
rminLEO = zminLEO+Re
rmaxLEO = zmaxLEO+Re

#Distances (m)
au = 149597870700

#Time (s)
day = 86400
hour = 3600
minute = 60

#Angles
omgES = 0.986 #deg/day
omgE = 7.292*10**-5 #rad/s
omgE_rad_d = 6.300387486749 #rad/d

#Others
G = 6.67384*10**(-11)
g = 9.80665
gL = 1.62200

#ISS characteristics 
mISS = 450*10**3
zISS = 408
rISS = zISS + Re
rhoISS = 3.8*10**-12 #kh/m3
CdISS = 2.07
areaISS = 1000   #m2
iIss = 51.6 #deg

#-------------------------------------------Converters--------------------------------------------
def rad(nb_deg):
    #Convert deg to rad
    return nb_deg*np.pi/180

def deg(nb_rad):
    #Convert rad to deg
    return nb_rad*180/np.pi

#-------------------------------------------Space Algorithms---------------------------------------

def e(v,z,gamma,R=Re,mue=mueE):
    #eccentricity 
    #/!\ gamma en radian !!
    return np.sqrt(np.sin(gamma)**2+((1-((z+R)*(v**2))/mue)**2)*np.cos(gamma)**2)

def e_r(rp,ra):
    return (ra-rp)/(ra+rp)

def a(v,z,R=Re):
    #semi-major axis
    return 1/(2/(z+R)-(v**2)/mueE)

def a_r(rp,ra):
    return 1/2*(rp+ra)

def a_z(zp,za,R=Re):
    return Re+1/2*(zp+za)

def a_circular(z,R=Re):
    return Re+z

def a_t(t,mue=mueE):
    return (mue*t**2/(4*np.pi**2))**(1/3)

def z(r,R=Re):
    return r-Re

def n(T):
    return 2*np.pi/T

def n_a(a,mue=mueE):
    return np.sqrt(mue/a**3)

def i(omg,e,r,n,R=Re,J=Je):
    return np.arccos(rad(omg/(-3/2*J*n*((R/r)**2)*(1/(1-e**2)**2))))

def i_circular(omg,J=Je):
    return np.arccos(omg/(-3*np.pi*J))

def i_man(vh,i1,i2):
    return 2*vh*np.sin((i2-i1)/2)

def omgt(i,e,r,n,R=Re,J=Je):
    #ATTENTION: i en rad pour le cos et n en degré par jour, donc omg en degré par jour
    return deg(-3/2*J*(n*3600*24)*((R/r)**2)*(1/((1-e**2)**2))*np.cos(i))

def rp(e,a):
    #perigee
    return a*(1-e)

def ra(e,a):
    #apogee
    return a*(1+e)

def v(r,a,mue=mueE):
    return np.sqrt(mue*(2/r-1/a))

def v_circular(z,R=Re,mue=mueE):
    v = np.sqrt(mueE/(z+R))
    return v

def vexit(r,mue=mueE):
    #Vitesse de sortie
    return np.sqrt(2*mue/r)

def vsat(r,mue=mueE):
    #Vitesse de satellisation en orbite circulaire
    return np.sqrt(mue/r)

#Homann Transfert
def HTdv1(r1,r2,mue=mueE):
    #r1 perigee and r2 apogee
    return np.sqrt(2*mue*r2/(r1*(r1+r2)))-np.sqrt(mue/r1)

def HTdv2(r1,r2,mue=mueE):
    #r1 perigee and r2 apogee
    return np.sqrt(mue/r2)-np.sqrt(2*mue*r1/(r2*(r1+r2)))

def v_orbit(v,r,mue=mueE):
    #Test de vitesse pour orbite circulaire
    vs = vsat(r,mue)
    ve = vexit(r,mue)
    if v >= ve : print("The satellite is too fast to be in orbit. \nVexit : ",ve)
    if v < vs : print("The satellite is too slow to be in orbit, it will crash !\nTry a speed higher than :",vs)
    if v<ve and v > (vs+0.01) : print("The satellite can be orbiting in an elliptic trajectory.")
    if v >= (vs -0.01) and v<=(vs+0.01):print("The satellite can be orbiting in a circular trajectory.")
   
    
def theta_init(r,a,e):
    #np.arccos((-1+r*(v**2)*(np.cos(sig)**2)/mue)*(1/e))
    #A verifier c'est chat qui a dit
    return np.arccos((a*(1-e**2)-r)/(e*r))

def theta_t(dt,T):
    return 2*np.pi*dt/T

def T(a,mue=mueE):
    #Période
    return 2*np.pi*np.sqrt(a**3/mue)

def T_format(t):
    # Calculer les composantes de la durée
    heures, reste = divmod(t, 3600)
    minutes, secondes = divmod(reste, 60)
    # Formater la durée sous forme de chaîne de caractères
    duree_formatee = "{:02d}:{:02d}:{:02f}".format(int(heures), int(minutes), secondes)
    
    return duree_formatee

def Fdrag(Cd,area,rho,v):
    return 1/2*Cd*area*rho*v**2

def Adrag(Cd,area,rho,v,m):
    return (1/m)*Fdrag(Cd,area,rho,v)

def radius(e,a,gamma,mue=mueE):

    # Paramètres de l'orbite
    Time =  T(a,mue)# Période orbitale (loi de Kepler)
    # Temps
    t = np.linspace(0, int(Time), 1000)
    
    # Calcul du rayon en fonction du temps
    theta0 = theta_init(a,a,e)
    theta = theta0 + (2 * np.pi * t) / Time
    r = (a * (1 - e**2)) / (1 + e * np.cos(theta - gamma)) / np.sqrt(mue)
    
    # Tracer le rayon en fonction du temps
    plt.figure(figsize=(8, 8))
    plt.plot(t, r)
    plt.title("Rayon au cours du temps")
    plt.xlabel("Temps (s)")
    plt.ylabel("km")
    plt.grid(True)
    plt.show()
    
    return r

def orbit(e,a,gamma,mue=mueE):
    period =  T(a,mue)
    t = np.linspace(0, int(period), 1000)
    theta0 = theta_init(a,a,e)
    theta = theta0 + (2 * np.pi * t) / period
    r = (a * (1 - e**2)) / (1 + e * np.cos(theta - gamma)) / np.sqrt(mue)
    
    x = r*np.cos(theta_init(r,a,e))
    y = r*np.sin(theta_init(r,a,e))
    plt.figure(figsize=(8, 8))
    plt.plot(x, y)
    plt.title("Rayon au cours du temps")
    plt.xlabel("Temps (s)")
    plt.ylabel("km")
    plt.grid(True)
    plt.show()
    
#---------------------------------TD2-----------------------------
#à refaire
def Mi(nb,stage,me,ms,mu,madd=0):
    mass = 0
    for i in range(nb):
        mass += me[i]+ms[i]
    return mass+mu+madd

def Mf(nb,stage,me,ms,mu,madd=0):
    massf= Mi(nb,stage,me,ms,mu,madd)-me[stage-1]
    return massf

def Dvp(ISP,Mi,Mf, go = g):
    return go*ISP*np.log(Mi/Mf)

def mprop(mu,Isp,dv,go=g):
    #mu mass of the payload, dv in m/s
    return mu*(1-np.exp(-dv/(go*Isp)))

#---------------------------------TD3----------------------------

def optimal_staging(k,ISP,dv,b2=3):
    #Lagrange multiplier method
    
    #definition du nombre d'étages
    n_stage = len(k)
    
    #Calcul des omega
    omg = []
    for i in range(n_stage):
        omg.append(k[i]/(k[i]+1))

    #Faire varier b2 jusqu'à avoir le bon dv = dv1+dv2. On part de b2 = 3 par défaut
    b = [b2]
    for i in range(1,n_stage):
        print(b[i-1])
        b.append( (1/omg[i-1])*(1-ISP[i]/ISP[i-1]*(1-omg[i]*b[i-1])) )
    b = b.reverse()

    
    #Calculer les dv
    Dv = np.zeros(n_stage)
    #a = np.zeros(n_stage)
    for i in range(n_stage):

        Dv[i] = 9.81*ISP[i]*np.log(b[i])
        #a[i] = (1+k[i])/b[i]-k[i]
     
    dvf = Dv[0]+Dv[1]
    
    return dvf



def opt_graph(bmax,k,ISP,dv,b=1,pas=0.1):
    dv = []
    bx = np.linspace(b,bmax,int(abs(b-bmax)/pas))
    
    while b<=bmax :
        dv.append(optimal_staging(k,ISP,dv,b))
        b +=0.1
    
    #plot graph
    plt.plot(bx,dv)
    plt.show()
    
    return 0

