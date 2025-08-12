#
# Cavitation
# # Ruhr University Bochum, 2022
#

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# parameter water 
sigmawater = 0.072  # surface tension in Nm^-2 
etawater = 0.001  # visocsity in Pa s 
pressureeq = 3e5 # pressure infty in Pa 
rhoeq = 1000 # density liquid in kg m^-3 
gamma = 1.95 # adiabatic exponent gas 
epsilon0 = 8.854e-12 # dielectric constant 
Bwater = 3000*1e5 # equation of state water 
nwater = 7 # equation of state water 
cwater0 = 1435 # sound velocity in ms^-1 

# start conditions 
R0 = 25e-6 # initial radius HV tip 
U0 = 0 # initial velocity 

pressure0 =  1e9 # pressure in the plasma after ignition *)
T0 = 1 # initial temperature a.u. 
pplasma = 3/(1.6e-19)*1e25;
T0 = pressure0/(3e28*1.38e-23)

# 10 ns voltage pulse *)
UkV0 = 20000 # (* voltage HV pin in V *)

def Efield(UkV):
    return UkV/R0

def pEfield(UkV,t):
    if t <1e-8:
        return 0.5*epsilon0*Efield(UkV)**2
    else:
        return 0
     

# pressure at the interface
def pgas(R,Rdot):
    return pressure0*(R0/R)**(3*gamma) - 2*sigmawater/R - 4*etawater/R*Rdot

# sound velcoiyt water 
Bwater = 3000*1e5
nwater = 7
def cwater(R,Rdot): 
  return cwater0*((pgas(R, Rdot) + Bwater)/(Bwater + pressureeq))**((nwater - 1)/(2*nwater))

# enthalpy 
def enthalpy(R,Rdot,t): 
  return (1/rhoeq*nwater/(nwater - 1)*(Bwater + 
     pressureeq)*(((pEfield(UkV0,t) + pgas(R,Rdot) + 
          Bwater)/(Bwater + pressureeq))**((nwater - 1)/nwater) - 1))

# solving the RP equation 
def func(y,t):
    R = y[0]
    Rp = y[1]
    Rpp = 1/R*1/(1-Rp/cwater(R,Rp))*(-1.5*(1-Rp/(3*cwater(R,Rp)))*Rp**2+(1+Rp/cwater(R,Rp))*enthalpy(R,Rp,t))
    return [Rp,Rpp]

y0 = [R0,U0]
t = np.linspace(1e-10,100e-6,1000)
y = odeint(func,y0,t)

gammaad = 1.95
R0 = 25e-6
T0 = 4000
p0 = 1e9
Temp = T0*(R0/y[:,0])**(3*(gammaad - 1))
p = p0*(R0/y[:,0])**(3*(gammaad))

fig = plt.figure()
ax = []
ax.append(fig.add_subplot())
ax[0].plot(t/1e-6,y[:,0]/1e-6)
ax[0].set(xlabel='t (s)',ylabel='R (micrometer)',xlim=(0,1e-4/1e-6),ylim=(0,500))

