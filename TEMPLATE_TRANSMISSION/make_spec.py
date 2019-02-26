import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.pyplot import *
import pickle
from matplotlib.ticker import FormatStrFormatter
from fm import *
rc('font',family='serif')

#load crosssections between wnomin and wnomax (in cm-1)
xsects=xsects_HST(6000,9400)

#the parameters
#planet/star system params--typically not free parameters in retrieval
Rp= 1.10#0.930#*x[4]# Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)
Rstar=1.20#0.598   #Stellar Radius in Solar Radii
M =1.0#1.78    #Mass in Jupiter Masses
#TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
Tirr=1400#1500#x[0]#544.54 #terminator **isothermal** temperature--if full redistribution this is equilibrium temp
logKir=-0.5  #TP profile IR opacity controlls the "vertical" location of the gradient
logg1=-1     #single channel Vis/IR opacity. Controls the delta T between deep T and TOA T
#Composition parameters---assumes "chemically consistnat model" described in Kreidberg et al. 2015
logMet=0.0#x[1]#1.5742E-2 #.   #Metallicity relative to solar log--solar is 0, 10x=1, 0.1x = -1 used -1.01*log10(M)+0.6
logCtoO=-0.26#x[2]#-1.97  #log C-to-O ratio: log solar is -0.26
logPQCarbon=-5.5  #CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to constant value at quench pressure value
logPQNitrogen=-5.5  #N2, NH3 Quench pressure--forces N2 and NH3 to ""  --ad hoc for chemical kinetics--reasonable assumption
#A&M Cloud parameters
logKzz=9 #log Rayleigh Haze Amplitude (relative to H2)
fsed=1.0 #haze slope--4 is Rayeigh, 0 is "gray" or flat.  
logPbase=1.5  #gray "large particle" cloud opacity (-35 - -25)
logCldVMR=-15.0 #cloud fraction
xRp=1.0
#simple 'grey+rayleigh' parameters
logKcld = -40
logRayAmp = -30
RaySlope = 0



#setting up "fake data" wavelength grid
wlgrid=np.arange(1.1, 1.6, 0.035)


#seting up input state vector. Must be in this order as indicies are hard wired in fx inside fm
x=np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M, logKzz, fsed,logPbase,logCldVMR, logKcld, logRayAmp, RaySlope])


#calling forward model
#thermochemical gas profile scaling factors
# 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18  19
#H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He  mmw
gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) #can be made free params if desired (won't affect mmw)
#gas_scale[13]=0.
#gas_scale[14]=0.
y_binned,y_mod,wnocrop,atm=fx(x,wlgrid,gas_scale, xsects)  #returns model spectrum, wavenumber grid, and vertical abundance profiles from chemistry


#pickle.dump([wlgrid, y_binned, wnocrop, y_mod, atm],open('Model_full_cold_02_10x_10gp.pic','wb'))


#creating fake data
y_meas=y_binned
err=np.zeros(len(y_meas))+50E-6
#pickle.dump([wlgrid, y_meas, err],open('data.pic','wb'))

ymin=np.min(y_binned)*1E2*0.9
ymax=np.max(y_binned)*1E2*1.1
fig1, ax=subplots()
xlabel('$\lambda$ ($\mu$m)',fontsize=18)
ylabel('(R$_{p}$/R$_{*}$)$^{2} \%$',fontsize=18)
minorticks_on()
errorbar(wlgrid, y_meas*100, yerr=err*100, xerr=None, fmt='Dk')
plot(wlgrid, y_binned*1E2,'ob')
plot(1E4/wnocrop, y_mod*1E2)
ax.set_xscale('log')
ax.set_xticks([0.3, 0.5,0.6,0.8,1, 2, 3, 4, 5, 7, 10])
ax.axis([0.3,5,ymin,ymax])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(length=10,width=1,labelsize='large',which='major')
savefig('spectrum.pdf',fmt='pdf')
show()
close()

#PLOTTING ATMOPHERE.............
#unpacking variables
#P is in bars
#T is in K
#H2O, CH4,CO,CO2,NH3,Na,K,TiO,VO,C2H2,HCN,H2S,FeH,H2,He are gas mixing ratio profiles
#qc is the condensate abundance profile given an "f_sed" value and cloud base pressure
#r_eff is the effective cloud droplet radius given (see A&M 2001 or Charnay et al. 2017)
#f_r is the mixing ratio array for each of the cloud droplet sizes.
P,T, H2O, CH4,CO,CO2,NH3,Na,K,TiO,VO,C2H2,HCN,H2S,FeH,H2,He,qc,r_eff,f_r=atm


fig2, ax1=subplots()
#feel free to plot whatever you want here....
ax1.semilogx(H2O,P,'b',ls='--',lw=2,label='H2O')
ax1.semilogx(CH4,P,'black',ls='--',lw=2,label='CH4')
ax1.semilogx(CO,P,'g',ls='--',lw=2,label='CO')
ax1.semilogx(CO2,P,'orange',ls='--',lw=2,label='CO2')
ax1.semilogx(NH3,P,'darkblue',ls='--',lw=2,label='NH3')
ax1.semilogx(Na,P,'b',lw=2,label='Na')
ax1.semilogx(K,P,'g',lw=2,label='K')
ax1.semilogx(TiO,P,'k',lw=2,label='TiO')
ax1.semilogx(VO,P,'orange',lw=2,label='VO')
ax1.set_xlabel('Mixing Ratio',fontsize=20)
ax1.set_ylabel('Pressure [bar]',fontsize=20)
ax1.semilogy()
ax1.legend(loc=4,frameon=False)
ax1.axis([1E-9,1,100,1E-7])

#plotting TP profile on other x-axis
ax2=ax1.twiny()
ax2.semilogy(T,P,'r-',lw='4',label='TP')
ax2.set_xlabel('Temperature [K]',color='r',fontsize=20)
ax2.axis([0.8*T.min(),1.2*T.max(),100,1E-6])
for tl in ax2.get_xticklabels(): tl.set_color('r')
savefig('atmosphere.pdf',fmt='pdf')
show()
close()




pdb.set_trace()






