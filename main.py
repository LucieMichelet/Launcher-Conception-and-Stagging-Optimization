# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 16:46:14 2023

Members : Lucie Michelet, Audrey Duthil, Kim Lerabo

"""

import webbrowser
import spacemechanics as sp
import numpy as np

#%% Bibliothèque
def mass_check(mu,m,nb_stage):
    #Structural mass per stage (kg)
    min_ms = [500,200,200]
    max_ms = [100000,80000,50000]
    res = True
    
    #Mass in intervals per stage
    for i in range(nb_stage):
        if m[i] > min_ms[i] and m[i]<max_ms[i]:
            res = True
        else :
            res = False       
            print("Stage ",i+1,"mass not in the correct interval : ",m[i],"kg")
            break

    if res == True :            
        #tot mass > tot mass - stage(i)
        for i in range(1,nb_stage):  
            if m[i-1]>m[i]:
                #print("Stage ",i+1,"to",i+1,":",m[0],">",m[1])
                res = True
            else : 
                res = False
                            
        #tOT MASS > 1000 tons
        if m[0] > 1000000:
            res = False
            print("Total mass to heavy : ",m[0],"kg > 1000 000 kg")
            
        if res == True:
            print("Mass validated !")
    else : 
        print("Mass not validated.")

    #Ajouter le check pour un stage = 1 ou 3

def losses(z):
    #m/s
    if z > 200 and z < 1800:
        result = 2.452*(10**-3)*z**2+1.051*z+1387.5
        return result 
    else :
        print("The altitude z is not in the range [200,800]. Please try again.")
        
        
def mass_calcul(k,mu,a):
    if len(k)== 3 :
        Mi = [0,0,mu/a[2]]
        Mi[1] = Mi[2]/a[1]
        Mi[0] = Mi[1]/a[0]

    if len(k) == 2:
        Mi = [0,mu/a[1]]
        Mi[0] = Mi[1]/a[0]
        
    me = np.zeros(n_stage)
    ms = np.zeros(n_stage)

    for j in range(n_stage):
       me[j] = ((1-a[j])/(1+k[j]))*Mi[j] 
       ms[j] = k[j]*me[j]
   
    return Mi,me,ms

#%% Characteristic of stages and propellant 

ISP_solid = [266,295,295] #s        stage 1 only
ISP_petrol = [285,320,320] #s      stage 1,2,3
ISP_liquid = [0,450,450] #s           stage 2 and 3

K = [0.12,0.16,0.25]


#%%------------------------------Mission 1------------------------------------
#                                  CNES
#%% First step : Injection Requirement

#altitude at injection : zp
zp = 340
rp = sp.Re+zp

za = 340
ra = sp.Re+za
i = sp.rad(90)
mu = 230
lat = sp.rad(60.8)

ax = sp.a_z(zp,za)

#Velocity at injection : vp 
vp = sp.v(rp,ax)
print(vp,"km/s")

#launcher azimut is given by cos(i) = sin(az)*cos(lat)
az = np.arcsin(np.cos(i)/np.cos(lat))
print(sp.deg(az))

#%% Second step : Loses, initial velocity and dv requiered
vsat = sp.vsat(sp.Re+zp,sp.mueE)*1000

loss = losses(zp)
vi = sp.omgE*sp.Re*np.cos(i)*1000 #m/s
dv = vsat - vi + loss #m/s

print(loss,"m/s\t",vi,"m/s\t",dv,"m/s")

#%% Third step : Optimal Staging

#Chossing propelant  : solid + petrol
ISP = [ISP_solid[0],ISP_petrol[1]]

k = [K[0],K[1]]

n_stage = len(k)
omg = []

#Lagrange multiplier method
for i in range(n_stage):
    omg.append(k[i]/(k[i]+1))


#Faire varier b2 jusqu'à avoir le bon dv = dv1+dv2
b2 = 5.056
b1 = (1/omg[0])*(1-ISP[1]/ISP[0]*(1-omg[1]*b2))
b = [b1,b2]
Dv = np.zeros(n_stage)
a = np.zeros(n_stage)

for i in range(n_stage):
    Dv[i] = sp.g*ISP[i]*np.log(b[i])
    a[i] = (1+k[i])/b[i]-k[i]
 
dvf = Dv[0]+Dv[1]
print("\ndv required\t",dv,"\ndvf\t\t\t",dvf)


#%% 

Mi,me,ms = mass_calcul(k,mu,a)

mass_check(mu,ms,n_stage)

print("\nMass propellant stage 1 : ",me[0],"kg\nMass propellant stage 2 : ",me[1])
#%% Fourth step : Expected Orbit 

#semi-major axis : ax
#eccentricity : circular
e = 0
#inclination : i
#apogee et perigee : za and z^p
#Period
t = sp.T(ax)
print("semi-major axis : ",ax,"\neccentricity : ",e,"\nApogee alt : ",za,"km\nPerigee alt : ",zp,"km\nInclination : ",sp.deg(i),"°\nPeriod : ",sp.T_format(t))

#%% Fifth step : Check results 

webbrowser.open("http://josselin.desmars.free.fr/work/teaching/launcher/v23/")

# MISSION REUSSIE !!!

#%%------------------------------Mission 2------------------------------------
#                                Roscosmos

zp = 410
rp = sp.Re+zp
za = 410
ra = sp.Re+za
i = sp.rad(51.6)
mu = 2620
lat = sp.rad(46)

#%%
ax = sp.a_z(zp,za)

#Velocity at injection : vp 
vp = sp.v(rp,ax)
print("vp :",vp,"km/s")

#launcher azimut is given by cos(i) = sin(az)*cos(lat)
az = np.arcsin(np.cos(i)/np.cos(lat))
print("azimut :",sp.deg(az),"°")

#%% Second step : Loses, initial velocity and dv requiered
vsat = sp.vsat(sp.Re+zp,sp.mueE)*1000

loss = losses(zp)
vi = sp.omgE*sp.Re*np.cos(i)*1000 #m/s
dv = vsat - vi + loss #m/s

print(loss,"m/s\t",vi,"m/s\t",dv,"m/s")

#%% Third step : Optimal Staging

#Chossing propelant  : solid + petrol = petrol
ISP = [ISP_petrol[0],ISP_petrol[1]]

k = [K[1],K[1]]

n_stage = len(k)
omg = []

#Lagrange multiplier method
for i in range(n_stage):
    omg.append(k[i]/(k[i]+1))

#Faire varier b2 jusqu'à avoir le bon dv = dv1+dv2. On part de b2 = 3.2
b2 = 5.169
#b2 = (1/omg[1])*(1-ISP[2]/ISP[1]*(1-omg[2]*b3))
b1 = (1/omg[0])*(1-ISP[1]/ISP[0]*(1-omg[1]*b2))
b = [b1,b2]

Dv = np.zeros(n_stage)
a = np.zeros(n_stage)

for i in range(n_stage):
    Dv[i] = sp.g*ISP[i]*np.log(b[i])
    a[i] = (1+k[i])/b[i]-k[i]
 
dvf = Dv[0]+Dv[1]
print("dv :",dv,"\ndvf :",dvf)

#%%

Mi,me,ms = mass_calcul(k,mu,a)

mass_check(mu,ms,n_stage)

#%% Fourth step : Expected Orbit 

#semi-major axis : ax
#eccentricity : circular
e = 0.7
#inclination : i
#apogee et perigee : za and z^p
#Period
t = sp.T(ax)
print("semi-major axis : ",ax,"\neccentricity : ",e,"\nApogee alt : ",za,"km\nPerigee alt : ",zp,"km\nInclination : ",sp.deg(i),"°\nPeriod : ",sp.T_format(t))

#%% Fifth step : Check results 

webbrowser.open("http://josselin.desmars.free.fr/work/teaching/launcher/v23/")

#MISSION FAILED : PB PERIGEE INCLINAISON

#%%------------------------------Mission 3------------------------------------
#                                Eutelsat

zp = 300
rp = sp.Re+zp
za = sp.zGEO
ra = sp.rGEO


i = sp.rad(5.2)
mu = 1500
lat = sp.rad(5.2)


#Velocity at injection : vp  jusqu'à 300
vp = sp.v_circular(zp)
print("vp :",vp,"km/s")


#launcher azimut is given by cos(i) = sin(az)*cos(lat)
az = np.arcsin(np.cos(i)/np.cos(lat))
print("azimut :",sp.deg(az),"°")

#%% Second step : Loses, initial velocity and dv requiered
vsat = sp.vsat(sp.Re+zp,sp.mueE)*1000

loss = losses(zp)
vi = sp.omgE*sp.Re*np.cos(i)*1000 #m/s
dv = vsat - vi + loss


print(loss,"m/s\t",vi,"m/s\t",dv,"m/s")

#%% Third step : Optimal Staging

#Chossing propelant  : solid + petrol = petrol
ISP = [ISP_solid[0],ISP_petrol[1],ISP_petrol[2]]

k = [K[0],K[1],K[1]]

n_stage = len(k)
omg = []

#Lagrange multiplier method
for i in range(n_stage):
    omg.append(k[i]/(k[i]+1))

#Faire varier b2 jusqu'à avoir le bon dv = dv1+dv2. On part de b2 = 3.2
b3 = 2.904
b2 = (1/omg[1])*(1-ISP[2]/ISP[1]*(1-omg[2]*b3))
b1 = (1/omg[0])*(1-ISP[1]/ISP[0]*(1-omg[1]*b2))
b = [b1,b2,b3]

a = np.zeros(n_stage)
Dv = np.zeros(n_stage)

for i in range(n_stage):
    Dv[i] = sp.g*ISP[i]*np.log(b[i])
    a[i] = (1+k[i])/b[i]-k[i]
 
dvf = Dv[0]+Dv[1]+Dv[2]
print("dv :",dv,"\ndvf :",dvf)

#%%

Mi,me,ms = mass_calcul(k,mu,a)

mass_check(mu,ms,n_stage)

#%% Fourth step : Expected Orbit 

#semi-major axis : ax
ax = sp.a_z(zp,za)
#eccentricity : circular
e = sp.e_r(rp,ra)
#inclination : i
#apogee et perigee : za and z^p
#Period
t = sp.T(ax)
print("semi-major axis : ",ax,"\neccentricity : ",e,"\nApogee alt : ",za,"km\nPerigee alt : ",zp,"km\nInclination : ",sp.deg(i),"°\nPeriod : ",sp.T_format(t))

#%% Fifth step : Check results 

webbrowser.open("http://josselin.desmars.free.fr/work/teaching/launcher/v23/")

#%%------------------------------Mission 4------------------------------------
#                                  NASA

zp = 1681
rp = sp.Re+zp
za = 1681
ra = sp.Re+za
inc = sp.rad(103)
mu = 170
lat = sp.rad(28.5)

#Velocity at injection : vp = np.sqrt(mueE/(z+R))
vp = sp.v_circular(zp)
print("vp :",vp*1000,"m/s")

#launcher azimut is given by cos(i) = sin(az)*cos(lat)
az = np.arcsin(np.cos(inc)/np.cos(lat))
print("azimut :",sp.deg(az),"°")


#%% Second step : Loses, initial velocity and dv requiered
vsat = sp.vsat(sp.Re+zp,sp.mueE)*1000

loss = losses(zp)
vi = sp.omgE*sp.Re*np.cos(inc)*1000 #m/s
dv = vsat - vi + loss #m/s

print(loss,"m/s\t",vi,"m/s\t",dv,"m/s")

#%% Third step : Stagging optimization

#Solid - liquid - liquid: 3 stages
Isp = [ISP_solid[0],ISP_liquid[1],ISP_liquid[2]]
k = [K[0],K[2],K[2]]
n_stage = len(k)

omg = []

#Lagrange multiplier method
for i in range(n_stage):
    omg.append(k[i]/(k[i]+1))

b3 = 4.0674
b2 = (1/omg[1])*(1-Isp[2]/Isp[1]*(1-omg[2]*b3))
b1 = (1/omg[0])*(1-Isp[1]/Isp[0]*(1-omg[1]*b2))
b = [b1,b2,b3]

Dv = np.zeros(n_stage)
a = np.zeros(n_stage)
dvf = 0

for i in range(n_stage):
    Dv[i] = sp.g*Isp[i]*np.log(b[i])
    a[i] = (1+k[i])/b[i]-k[i]
    dvf += Dv[i]


print("dv :",dv,"\ndvf :",dvf)

#%% Test mass

Mi,me,ms = mass_calcul(k,mu,a)

mass_check(mu,ms,n_stage)

print("\nMass propellant stage 1 : ",me[0],"kg\nMass propellant stage 2 : ",me[1],"kg\nMass propellant stage 3 : ",me[2],"kg")

#%% Fourth step : Expected Orbit 

#semi-major axis : ax
ax = sp.a_z(zp,za)
#eccentricity : circular
e = sp.e_r(rp,ra)
#inclination : i
#apogee et perigee : za and z^p
#Period
t = sp.T(ax)
print("semi-major axis : ",ax,"\neccentricity : ",e,"\nApogee alt : ",za,"km\nPerigee alt : ",zp,"km\nInclination : ",sp.deg(inc),"°\nPeriod : ",sp.T_format(t))

#%% Fifth step : Check results 

webbrowser.open("http://josselin.desmars.free.fr/work/teaching/launcher/v23/")
