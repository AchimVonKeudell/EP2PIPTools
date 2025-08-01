# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 13:56:15 2023

@author: Achim
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import optimize

datafilev = 'C2--75V-2ndclosest-fullspec.txt'
datafilei = 'C1--75V-2ndclosest-fullspec.txt'

datafilev = 'C2.dat'
datafilei = 'C1.dat'

fexcitation = 20000
tstart = 70e-6

datastart = 200
dataend = 750

datavorg = np.loadtxt(datafilev,usecols=(0,1),skiprows=5)
dataiorg = np.loadtxt(datafilei,usecols=(0,1),skiprows=5)
datat = datavorg[datastart:dataend,0]
datav = datavorg[datastart:dataend,1]
datai = dataiorg[datastart:dataend,1]
dt = datat[10]-datat[9]
print(' Time step:',dt)

fig, ax = plt.subplots(2,1)
ax[0].plot(datat/1e-6,datav,label='voltage (kV)')
ax[0].plot([tstart/1e-6,tstart/1e-6],[-8000,8000],linestyle='dashed',color='b')
ax[0].set(xlabel='time (microseconds)',ylabel='Voltage (V)')

ax[0].legend(loc=1, fontsize = 6)
axp = ax[0].twinx()
axp.plot(datat/1e-6,datai/1e-3,label='current (mA)',color='orange')
axp.set(ylabel='current (mA)')
axp.set(ylim=[-50,50])

uderiv = savgol_filter(datav,51,4,deriv=1,delta=dt)

def f(c):
    sum = 0
    for i in range(len(datat)):
        if datat[i]>tstart:
           sum += abs(datai[i]-uderiv[i]*c)
    return sum       

sol = optimize.minimize(f,0, method='Nelder-Mead',tol=1e-8) 
csys = sol.x

isystem = uderiv*csys

axp.plot(datat/1e-6,isystem/1e-3,label='current model (A)',color='m')
axp.legend(loc=4,fontsize = 6)

power = datav*(datai-isystem)
powerint = sum(power*dt)*fexcitation

ax[1].plot(datat/1e-6,abs(power),label='power (W)',color='r')
ax[1].legend(fontsize=6)
ax[1].set(ylabel='power (W)',xlabel='time (microseconds)',ylim=(0,200))

print('C_system: ',"{:.6f}".format(csys[0]/1e-12),' pF')
print(' V_pp: ',"{:.0f}".format(max(datav)-min(datav)),' V')
print(' <Power>: ',"{:.3f}".format(abs(powerint)),' W')

fig.suptitle('Power: ' + "{:.3f}".format(abs(powerint))+' W, '+
             'V_pp: ' + "{:.0f}".format(max(datav)-min(datav))+' V, '+
             'C_system: ' + "{:.2f}".format(csys[0]/1e-12)+' pF, '+
             'f: ' + "{:.1f}".format(fexcitation/1e3)+' kHz')

np.savetxt('powervectors.dat',np.transpose([datat/1e-6,abs(power),datav,datai,isystem]),header='t P V I Isystem',comments='')

