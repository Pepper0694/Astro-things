import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



ZAMS5 = pd.read_csv('5.0MZAMS.csv')
n15 = pd.read_csv('n1.5.csv')
n30 = pd.read_csv('n3.0.csv')


#Problem 1
#-------------------
logR = ZAMS5['logR']
logRho = ZAMS5['logRho']
R = (7*10**10)*(10**logR)
Rho = 10**logRho
Rhoc = max(Rho)
alpha15 = max(R)/3.65375
alpha3 = max(R)/6.89685

R15 = alpha15*n15['xi']
Rho15 = Rhoc*(n15['theta']**1.5)
R3 = alpha3*n30['xi']
Rho3 = Rhoc*(n30['theta']**3)


plt.loglog(R,Rho,label='5M Star')
plt.loglog(R15,Rho15,label = 'n = 1.5')
plt.loglog(R3,Rho3,label = 'n = 3.0')

plt.xlabel('Radius (cm)')
plt.ylabel('Density (g/cm^3)')
plt.title('Radius vs. Density')
plt.legend(loc = 'best')
plt.show()


#Problem 2
#------------------

Pgas =10**ZAMS5['logPgas']
P = 10**ZAMS5['logP']
Prad = P - Pgas

plt.loglog(R, Pgas, label = 'Gas Pressure')
plt.loglog(R, Prad, label = 'Radiation Pressure')
plt.xlabel('Radius (cm)')
plt.ylabel('Pressure (dyn/cm^2)')
plt.title('Radius vs. Pressure')
plt.legend(loc='best')
plt.show()


#Problem 3
#----------------------

pp = ZAMS5['pp']
cno = ZAMS5['cno']
total = pp + cno

plt.loglog(R,pp,label='PP Energy Generation Rate')
plt.loglog(R,cno,label='CNO Energy Generation Rate')
plt.loglog(R,total,label='Total Energy Generation Rate')
plt.xlabel('Radius (cm)')
plt.ylabel('Energy Generation Rate (erg/sec/g)')
plt.title('Radius vs. Energy Generation')
plt.legend(loc='best')
plt.show()


#Problem 4
#----------------------

Aladin = ZAMS5['gamma1']

plt.plot(R,Aladin)
plt.xscale('log')
plt.xlabel('Radius (cm)')
plt.ylabel('Adiabatic Index (dlnP_dlnRho at constant S)')
plt.title('Radius vs. Adiabatic Index')
plt.show()



#Problem 5
#-----------------------

rad = ZAMS5['gradr']
conv = ZAMS5['gradT']
adia = ZAMS5['grada']
act = ZAMS5['actual_gradT']

plt.plot(R,act,lw =8,label = 'Actual Temp. Grad.')
plt.plot(R,rad,lw=4, label = 'Radiation Temp. Grad.')
plt.plot(R,adia,lw=3,label = 'Adiabatic Temp. Grad.')
plt.plot(R,conv,'m--',lw=4,label = 'Convection Temp. Grad.')

plt.xscale('log')
plt.xlabel('Radius (cm)')
plt.ylabel('Temperatre Gradient Values')
plt.title('Radius vs. Temperature Gradients')
plt.legend(loc='best')
plt.show()


