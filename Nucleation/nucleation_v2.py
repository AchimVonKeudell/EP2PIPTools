# #############################
#
# Cluster Coagulation Simulation
# after eq. 22 in 
# Marian von Smoluchowski
# Zeitschrift fuer physikalische Chemie 92, S. 129-168 (1917) 
#  > 7200 Zitate :)
# 
# ############################

from numba import njit
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation


# read temperature dependence from cavitation model
datatemp = np.loadtxt('cavitation.dat',usecols=(0,3),skiprows=1).T
print(datatemp[1])

r0 = 0.3    # lattice constant tungsten   
maxparticler = 5 # max particlediameter in nm
D = 1*1e-4 # diffusion constant metal in vapor m^2/s
DR = D*r0*1e-9 # Transportfaktor 
print('DR:',DR)

# hot spot estimate
# number of atoms in a volume 12 micrometer diameter x 5 nm height
# footprint of the plasma at the surface which is removed per pulse
nbspot = 6e-6**2*np.pi*5e-9/(r0*1e-9)**3
# esteimate of the starting seed particles
# this number exisits initially in a sphere with radius 25 micrometer at the tip
n0 = nbspot/(4*np.pi/3*(25e-6)**3)
print('Initial density after pulse : ',n0,' m^-3')             

# supersaturation value estimate
# tungsten partial pressure after pulse
pspot = n0/2.4e25 # atm
print(' pressure W: ',pspot, ' (bar)')

# tungsten vapor pressure tungsten @ 4000 K
# J Res Natl Bur Stand A Phys Chem. 1965 Sep-Oct; 69A(5): 417–421.
# doi: 10.6028/jres.069A.044
# Vapor Pressure and Heat of Sublimation of Tungsten
# R. Szwarc,2 E. R. Plante, and J. J. Diamond
pvapor = 10**(7.93-45087/4000) # atm
print(' vapor pressure W: ',pvapor, ' (bar)')

# definition of vectors
maxnb = int(4*np.pi/3*(maxparticler*0.5)**3/r0**3)
print('length of particle i vector: ',maxnb)
n = np.zeros(maxnb) # densities cluster
r = np.ones(maxnb) # radii cluster
for i in range(maxnb):
  r[i] = i**(1/3)*r0  
dn = np.zeros(maxnb)
n[1] = n0


# coagulation algorithm
# @njit(nopython=True)
def dn_dt(n,t):

    temp = np.interp(t,datatemp[0],datatemp[1])  
    # print('Temp: ',temp,' (K)')

    # total density of particles    
    particles = 0
    for j in range(1,maxnb):
          particles += n[j]
 
    # cogaulation sequence
    # smalles bin i=1
    dn[1] = -n[1]*particles*DR
    # bins i=2 -> maxnb
    for i in range(2,maxnb):
      prod = 0      
      j = 1
      # prodcution of new particle via coagulation
      while j < i:
          # two smaller particles combine and build one new particle at size i
          # with probaility 1 per time step
          prod += 0.5*n[j]*n[i-j]
          j += 1
      # particles i are lost via reactions with all other particles     
      dn[i] = (prod - n[i]*particles)*DR
    return dn

# Define plot
fig, ax = plt.subplots(1,1,figsize=(6,6))  
line1 = ax.loglog(r[1:],n[1:],lw=1,linestyle='solid',color='magenta',label='N(r) (norm.)')[0] 
rm=[0.3,0.3]
nrm=[1e-4,2]
line2 = ax.loglog(rm,nrm,lw=1,linestyle='dashed',color='b',label='r @max')[0] 
ax.set(ylim=(1e-4,2),xlim=(1e-1,10))
ax.set(ylabel='N(r) (norm.)',xlabel='r radius cluster (nm)')
ax.legend()


# Define timing
dt = 1e-11
t = 0
tend = 1e-7
microsteps = 20
framesNB=int(tend/(dt*microsteps))

def animate(k):
     global dn,n,t
  
     for f in range(microsteps):
       t += dt
       dn = dn_dt(n,0)
       n += dn*dt
 
     natoms = 0
     for i in range(1,maxnb):
        natoms += i*n[i] 
  
     print('Time: ',int(t/1e-12),' (ps), density total: ',"{:.2e}".format(natoms),' (m^-3)')
     maxn = max(n[1:])
     rmax = r[np.argmax(n[1:])+1]
     # print('rmax ',rmax,' (nm)')
     line1.set_ydata(n[1:]/maxn)
     line2.set_xdata([rmax,rmax])
     
     fig.suptitle('time: ' + "{:.2f}".format(t/1e-9)+' (ns), n_max: '+ "{:.0e}".format(maxn)+' m^-3, r_@max: '+ "{:.1f}".format(rmax)+' nm')

        
anim = animation.FuncAnimation(fig,animate,interval=0,frames=framesNB)
anim.save('nucleation.gif',fps=25,dpi=300)

#np.savetxt('nucleation.dat',np.transpose([r[1:],n[1:]/max(n[1:])]),fmt='%.12f',header='r (n(r))',comments='')
            