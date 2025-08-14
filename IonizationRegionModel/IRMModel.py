#   Ionization Region Modell
#
#   from Gudmundsson et al. PSST 25, 065004 (2016)
#
#   IRMGlobale Variablen

import numpy as np

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

import matplotlib.pyplot as plt

#global beta,Tg,UIR,VIR,SRT,SBP,L,ret,tend,tpulse
#global mar,mti,ID0,UD0,ID,UD,pgas,kb,ngas,el

model = 'Ti'
#model = 'Al'
#model = 'Cr'

if model == 'Ti':
    mamu = 48
    
    gammamet = 0.0 # sekundärelektronen Ti+ on Ti
    gammaarp = 0.06 # sekundärelektronen Ar+ on Ti from Hagstrum model
    
    metname = 'Ti'
    metpname = 'Ti+'
    
    # Fitting Parameter IRM
    ret = 0.5
    beta = 0.9
    UIR = 25
    #ID0 = 55 # applied maximum current
    #UD0 = 500 # applied voltage

    tpulse = 100*1e-6 # pulse length in microseconds 
    tbegin = 15*1e-6 # pulse length in microseconds 
    
    ionizationdegreestart = 5e-2
    
    outputfile = 'IRMModelTi.dat'
    datane = 'DataTi/probe-Ti.txt'
    datapower = 'DataTi/IV-Ti.txt'
    
    # Dimension IRM
    z1 = 0.001
    z2 = 0.025 # hier größer 
    L = z2-z1;
    #rc1 = 0.005
    #rc2 = 0.020
    widthrc = 0.010
    rcpos = 0.013
    #rc1 = 0.009
    #rc2 = 0.015
    rc1 = rcpos-widthrc/2
    rc2 = rcpos+widthrc/2
    
    ymax = 2e20
    
if model == 'Al':
    mamu = 27
    
    gammamet = 0.0 # sekundärelektronen Al+ on Al
    gammaarp = 0.049 # sekundärelektronen Ar+ on Al from Hagstrum
    
    metname = 'Al' 
    metpname = 'Al+'
    
    # Fitting Parameter IRM
    ret = 0.85
    beta = 0.9
    UIR = 14
    #ID0 = 38 # applied maximum current
    #UD0 = 600 # applied voltage

    tpulse = 100*1e-6 # pulse length in microseconds     
    tbegin = 15*1e-6 # pulse length in microseconds 
    
    ionizationdegreestart = 8e-2
    
    outputfile = 'IRMModelAl.dat'
    datane = 'DataAl/probe-Al.txt'
    datapower = 'DataAl/IV-Al.txt'
    
     # Dimension IRM
    z1 = 0.001
    z2 = 0.015
    L = z2-z1;
    rc1 = 0.007
    rc2 = 0.015
    widthrc = 0.007
    rcpos = 0.013
    #rc1 = 0.009
    #rc2 = 0.015
    rc1 = rcpos-widthrc/2
    rc2 = rcpos+widthrc/2
    
    
    ymax = 1.5e20
    
if model == 'Cr':
    mamu = 52
    
    gammamet = 0.0 # sekundärelektronen Cr+ on Cr
    gammaarp = 0.1 # sekundärelektronen Ar+ on Cr ??
    
    metname = 'Cr'
    metpname = 'Cr+'
    
    # Fitting Parameter IRM
    ret = 0.5
    beta = 0.90
    UIR = 15

    #ID0 = 30 # applied maximum current
    #UD0 = 750 # applied voltage
    
    tpulse = 150*1e-6 # pulse length in microseconds 
    tbegin = 15*1e-6 # pulse length in microseconds 
    
    ionizationdegreestart = 3e-2
    
    outputfile = 'IRMModelCr.dat'
    datane = 'DataCr/probe-Cr.txt'
    datapower = 'DataCr/IV-Cr.txt'

    # Dimension IRM
    z1 = 0.001
    z2 = 0.015
    L = z2-z1;
    widthrc = 0.006
    rcpos = 0.013
    #rc1 = 0.009
    #rc2 = 0.015
    rc1 = rcpos-widthrc/2
    rc2 = rcpos+widthrc/2
    
    ymax = 5e20


# ionization arte from Lennon et al.
# J. Phys Chem. Ref. Data 17, 1285, (1988)
def ionizationrate(I,a0,a1,a2,a3,a4,a5,Te):
    return (np.exp(-I/Te)*(Te/I)**0.5*
                        (a0*(np.log10(Te/I))**0+
                         a1*(np.log10(Te/I))**1+
                         a2*(np.log10(Te/I))**2+
                         a3*(np.log10(Te/I))**3+
                         a4*(np.log10(Te/I))**4+
                         a5*(np.log10(Te/I))**5))

def kizmetal(te):
    if model == 'Ti':
        return 2.8278e-13*te**-0.0579*np.exp(-8.7163/te)
    if model == 'Al':
        return ionizationrate(6,3.0764e-7,-2.8504e-7,1.0518e-7,-2.1274e-8,-2.5008e-8,2.0889e-8,te)*1e-6
        #return 2.8278e-13*te**-0.0579*np.exp(-8.7163/te)
    if model == 'Cr':
        return ionizationrate(6,1.7476e-7,-9.4892e-8,-3.9137e-8,2.2035e-8,1.0542e-8,-4.5512e-9,te)*1e-6
        #return 2.8278e-13*te**-0.0579*np.exp(-8.7163/te)
    
def eizmetal():
    if model == 'Ti':
        return 6.628
    if model == 'Al':
        return 5.98
    if model == 'Cr':
        return 6.758
    
# yields from SRIM    
def yieldmetal():
    if model == 'Ti':
       return 0.51 # Trim Yield Ti on Ti @ 500 eV
    if model == 'Al':
       return 1.12 # Trim Yield Al on Al @ 600 eV 
    if model == 'Cr':
       return 1.65 # Trim Yield Cr on Cr @ 750 eV

def yieldargon():
    if model == 'Ti':
       return 0.65 # Trim Yield Ar on Ti @ 500 eV
    if model == 'Al':
       return 0.84 # Trim Yield Ar on Al @ 600 eV
    if model == 'Cr':
       return 1.79 # Trim Yield Ar on Cr @ 750 eV


Tg = 300
kb = 1.3806e-23
el = 1.6021e-19
mar = 40*1.6e-27
mmet = mamu*1.6e-27


# thermal velocities
var = np.sqrt(2*kb*Tg/(np.pi*mar))
vmet = np.sqrt(2*kb*Tg/(np.pi*mmet))

# Dimensions IRM, volumes and surface areas
VIR = L*rc2**2*np.pi-L*rc1**2*np.pi
SRT = rc2**2*np.pi-rc1**2*np.pi
SBP = L*2*rc2*np.pi + SRT + L*2*rc1*np.pi

# experimental parameters
tend = 200*1e-6 # time window simulation in microseconds


# Anfangsbedingungen
pgas = 0.5 # in Pa
ngas = pgas/(kb*Tg)
c0 = np.zeros(5)

# initial condition
c0[0] = ngas # Ar
c0[1] = ngas*ionizationdegreestart # Arp
c0[2] = ngas*1e-12 # Metal
c0[3] = ngas*1e-12 # Metal Ions
c0[4] = 1 # Te in (eV)


tm, Vm, Im, _ = np.loadtxt(datapower).T

Im = savgol_filter(Im, 51, 3) # smooth current signal
Im[Im<0] = 0 # remove negative numbers (due to noise)
Im[tm<7] = 0 # remove switiching noise around t = 0
If = interp1d(tm*1e-6, Im, bounds_error=False, fill_value=0)

Vm = -savgol_filter(Vm, 11, 3)
Vm[Vm<0] = 0 # remove negative numbers (due to noise)
Uf = interp1d(tm*1e-6, Vm, bounds_error=False, fill_value=0)


# -----------------------------------------------------------
# Balance equations for the IRM
# -----------------------------------------------------------
def IRMFunc(t,c):
    nar = c[0]
    narp = c[1]
    nmet = c[2] 
    nmetp = c[3]
    te = c[4]

    # assuming a linear increase in current during the pulse as input,
    # here data could be directly chosen   
    # at the end the calculated current is also an output quantity
    # in terms of self consistency I_calc shoud coincide with ID0*t/tpulse
    
    # if t<tpulse and t>tbegin: 
    #   ID = ID0*t/tpulse
    #   UD = UD0
    #else:
    #   ID = 0
    #   UD = UD0
       
    ID = If(t)
    UD = Uf(t)
   
    # Fixing quasineutrality 
    ne = narp + nmetp

    # ionization rates argon
    kizar = 2.34e-14*te**0.59*np.exp(-17.44/te);
    eizar = 15.8
    
    # ionization rates metal defined via functions above
    kizmet = kizmetal(te)
    eizmet = eizmetal()
    
    # mean free path
    sigma = 10e-19
    lambdag = 1/(nar*sigma)
    Fcoll = 1- np.exp(-L/lambdag)

    #flux densities ar ions towards racetrack and bulk plasma
    fluxarrt = beta * narp * np.sqrt(el*UIR/mar)
    fluxarbp = SRT/SBP*(1/beta-1)*fluxarrt

    # dar
    dar  = (0.5*var*(ngas-nar)/VIR*SBP - ne*nar*kizar 
            + fluxarrt/VIR*SRT- 
            0.5*vmet/L*Fcoll*mmet/mar*(nmet+nmetp)/(nar+narp)*nar) 

    # darp
    darp  = ne*nar*kizar - 1/VIR*(fluxarrt*SRT+fluxarbp*SBP)

    # flux densities metal ions towards racetrack and bulk plasma
    fluxmetrt = beta * nmetp * np.sqrt(el*UIR/mmet)
    fluxmetbp = SRT/SBP*(1/beta-1)*fluxmetrt
 
    # dmetal ions
    dmetp  = ne*nmet*kizmet - 1/VIR*(fluxmetrt*SRT+fluxmetbp*SBP) 

    if t<tpulse and t>tbegin:
       yieldar = yieldargon() # Trim Yield Ar on Ti
       yieldmet = yieldmetal()
    else:
       yieldar = 0 # Trim Yield Ar on Ti
       yieldmet = 0 # Trim Yield Ti on Ti
 
    # dmetal neutrals
    dmet = - ne*nmet*kizmet - 0.5*vmet*nmet/VIR*SBP + SRT/VIR*(fluxarrt*yieldar+fluxmetrt*yieldmet) 

    # balance equation for electron temperature
    dte = 2/3*1/(el*ne)*(0.5*UIR/UD*ID*UD/VIR-el*ne*nar*kizar*eizar-el*ne*nmet*kizmet*eizmet)

    return np.array([dar,darp,dmet,dmetp,dte])


# Solving equation system
tspan = (tbegin,tend)
sol = solve_ivp(IRMFunc,tspan,c0,'RK45') #,rtol=1e-7,atol=1e-7
y = sol.y
t = sol.t

# plot densities
fig, ax = plt.subplots(2,2,figsize=(9,9))
fig.tight_layout(pad=5.0)
#ax[0].set_yscale('log')
ax[0,0].plot(t/1e-6,y[0,:],label='Ar')
ax[0,0].plot(t/1e-6,y[1,:],label='Ar+')
ax[0,0].plot(t/1e-6,y[2,:],label=metname)
ax[0,0].plot(t/1e-6,y[3,:],label=metpname)
ax[0,0].plot(t/1e-6,y[1,:]+y[3,:],label='ne',linestyle='dashed')

data = np.loadtxt(datane,skiprows=1).T     
ax[0,0].plot(data[0,:],data[1,:],label='n_e (data)',lw=0,marker='o',markersize=1)
ax[0,0].set(xlabel='time (microseconds)',ylabel='densities (m^-3)',xlim=(-20,tend/1e-6+100),ylim=(0,ymax))
ax[0,0].legend(fontsize = 6,loc = 1)

# plot temperature
data = np.loadtxt(datane,skiprows=1).T     
ax[0,1].plot(data[0,:],data[2,:],label='T_e (data)',lw=0,marker='o',markersize=1)
ax[0,1].plot(t/1e-6,y[4,:],label='Te (eV)')
ax[0,1].set(xlabel='time (microseconds)',ylabel='temperature (eV)',ylim=(0,12),xlim=(-20,tend/1e-6+100))
ax[0,1].legend(fontsize = 6, loc = 1)

# info box
infobox = ''
infobox += 'model: ' + model +'\n'
infobox += 'beta: ' + "{:.2f}".format(beta) + '\n'    
infobox += 'ret: ' + "{:.2f}".format(ret) + '\n'    
infobox += 'UIR: ' + "{:.1f}".format(UIR) + ' (V) \n'    
infobox += 'EfieldIR: ' + "{:.3f}".format(UIR/L) + ' (V/m) \n'
infobox += 'S/V: ' + "{:.3f}".format((SBP+SRT)/VIR) + ' (1/m)'
#infobox += 'Id0: ' + "{:.1f}".format(ID0) + ' (A)'    
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax[0,1].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0,1].transAxes)

# plot modelled current and power
#UD = UD0
#ID = ID0
fluxarrt = beta * y[1,:] * np.sqrt(el*UIR/mar);
fluxmetrt = beta * y[3,:] * np.sqrt(el*UIR/mmet);
Icalc = el*fluxarrt*SRT*(1 + (1 - ret)*gammaarp)+el*fluxmetrt*SRT*(1 + (1 - ret)*gammamet)
Pcalc = (fluxarrt*SRT*(1-ret)*gammaarp*(Uf(t)-UIR)*el + 
         fluxmetrt*SRT*(1-ret)*gammamet*(Uf(t)-UIR)*el + 
         If(t)*UIR + fluxarrt*SRT*(Uf(t)-UIR)*el + fluxmetrt*SRT*(Uf(t)-UIR)*el)

data = np.loadtxt(datapower,skiprows=1).T     
ax[1,0].plot(data[0,:],data[2,:],label='I (data)',lw=0,color='m',marker='o',markersize=1)
ax[1,0].plot(data[0,:],data[3,:],label='P (data)',lw=0,color='orange',marker='o',markersize=1)
ax[1,0].plot(t/1e-6,Icalc,label='I_calc (A)',color='m')
ax[1,0].plot(t/1e-6,Pcalc/1e3,label='P_calc (kW)',color='orange')
ax[1,0].set(xlabel='time (microseconds)',ylabel='I_calc (A), P_calc (kW)',xlim=(-20,tend/1e-6+100),ylim=(0,60))
ax[1,0].legend(fontsize = 6, loc = 4)

ax2 = ax[1,0].twinx()
ax2.set(ylabel='voltage')
ax2.plot(data[0,:],data[1,:],label='V (data)',lw=0,color='g',marker='o',markersize=1)


IRMrightx = [rc1,rc2]
IRMleftx = [-rc1,-rc2]
IRMtopy = [z2,z2]
IRMbottomy = [z1,z1]

symmetryx = [0,0]
symmetryy = [-5*1e-3,50*1e-3]
ax[1,1].plot(symmetryx,symmetryy,color='b',lw=1,linestyle='dashed')

racetrackleftx = [-rcpos,-rcpos]
racetrackrightx = [rcpos,rcpos]
racetracky = [0,z2]
ax[1,1].plot(racetrackleftx,racetracky,color='r',lw=1,linestyle='dashed',label='pos. racetrack')
ax[1,1].plot(racetrackrightx,racetracky,color='r',lw=1,linestyle='dashed')

targetx = [-25.4*1e-3,25.4*1e-3]
targety = [0,0]
ax[1,1].plot(targetx,targety,color='black',lw=3,label='2" target')

ax[1,1].fill_between(IRMrightx,IRMtopy,IRMbottomy,color='none',hatch='X',edgecolor='b',label='IR')
ax[1,1].fill_between(IRMleftx,IRMtopy,IRMbottomy,color='none',hatch='X',edgecolor='b')
ax[1,1].set_yticks([])
ax[1,1].set(xlabel='r (m)',ylabel='z (m)',xlim=(-30*1e-3,30*1e-3),ylim=(-5*1e-3,+50*1e-3))
ax[1,1].text(rc1-0.003,z2+0.001,'rc_1')
ax[1,1].text(rc2-0.003,z2+0.001,'rc_2')
ax[1,1].text(rc2+0.001,z2-0.002,'z_2')
ax[1,1].text(rc2+0.001,z1+0.001,'z_1')
ax[1,1].legend(fontsize=8)

plt.show()

np.savetxt(outputfile,np.transpose([t/1e-6,y[0,:],y[1,:],y[2,:],y[3,:],y[1,:]+y[3,:],y[4,:],Icalc,Pcalc/1e3]),fmt='%.8e',header='t Ar Ar+ Met Met+ ne Te Icalc Pcalc',comments='')


