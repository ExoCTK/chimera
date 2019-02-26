from matplotlib import pyplot as plt
from dynesty import plotting as dyplot
import pickle


#PLOTTING UP CORNER PLOT----------------------------------------------------------------
truth=[1300,0,-0.26,9, 1, 1.5, -15, 1]  #values to overplot on histograms
labels=['Tirr', 'logMet', 'logCtoO', 'logKzz', 'fsed' ,'logPbase','logCldVMR','xRp']

#import past run
#samples=pickle.load(open('dyn_output_100LP.pic','rb')) #an example 100 live point run
samples=pickle.load(open('dyn_output_500LP.pic','rb'))
#samples=pickle.load(open('dyn_output_1000LP.pic','rb'))  #an example 1000 live point run

#printing evidence:
print('ln(Z)= ', samples.logz[-1])

# corner plot
fig, axes = dyplot.cornerplot(samples,smooth=0.05, color='red',show_titles=True, labels=labels, truths=truth,title_kwargs={'y': 1.04}, fig=plt.subplots(8, 8, figsize=(12, 12)))
plt.savefig('stair_pairs.pdf',fmt='pdf')
plt.show()
plt.close()

#GENERATING RANDOMLY SAMPLES SPECTRA-----------------------------------------------------
import numpy as np
xsecs=xsects_HST(2000, 30000)

Nspectra=200

#loading in data again just to be safe
wlgrid, y_meas, err=pickle.load(open('data.pic','rb'))

#setting up default parameter values--SET THESE TO SAME VALUES AS IN LOG-LIKE FUNCTION
#planet/star system params--xRp is the "Rp" free parameter, M right now is fixed, but could be free param
Rp= 1.10   # Planet radius in Jupiter Radii--this will be forced to be 1 bar radius--arbitrary (scaling to this is free par)
Rstar=1.20  #Stellar Radius in Solar Radii
M =1.0   #Mass in Jupiter Masses
xRp=1.0  #scaling factor to Radius

#TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
Tirr=1400     #Irradiation temperature as defined in Guillot 2010
logKir=-0.5  #TP profile IR opacity controlls the "vertical" location of the gradient
logg1=-1     #single channel Vis/IR opacity. Controls the delta T between deep T and TOA T

#Composition parameters---assumes "chemically consistent model" described in Kreidberg et al. 2015
logMet=0.0  #.   #Metallicity relative to solar log--solar is 0, 10x=1, 0.1x = -1: valid range is -1.5 - 3.0
logCtoO=-0.26   #log C-to-O ratio: log solar is -0.26: valid range is -1.0 - 0.3
logPQCarbon=-5.5  #CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to constant value at quench pressure value: valid range -6.0 - 1.5
logPQNitrogen=-5.5  #N2, NH3 Quench pressure--forces N2 and NH3 to ""

#Ackerman & Marley 2001 Cloud parameters--physically motivated with Mie particles
logKzz=9 #log Kzz (cm2/s)--valid range: 2 - 11 -- higher values make larger particles
fsed=1.0 #sediminetation efficiency--valid range: 0.5 - 5--lower values make "puffier" more extended cloud
logPbase=1.5  #cloud base pressure--valid range: -6.0 - 1.5
logCldVMR=-15.0 #cloud condensate base mixing ratio (e.g, see Fortney 2005)--valid range: -15 - -2.0

#simple 'grey+rayleigh' parameters just in case you don't want to use a physically motivated cloud
#(most are just made up anyway since we don't really understand all of the micro-physics.....)
logKcld = -40  #uniform in altitude and in wavelength "grey" opacity (it's a cross-section)--valid range: -50 - -10
logRayAmp = -30  #power-law haze amplitude (log) as defined in des Etangs 2008 "0" would be like H2/He scat--valid range: -30 - 3
RaySlope = 0  #power law index 4 for Rayleigh, 0 for "gray".  Valid range: 0 - 6


#weighting the posterior samples for appropriate random drawing
from dynesty import utils as dyfunc
samp, wts = samples.samples, np.exp(samples.logwt - samples.logz[-1])
samples2 = dyfunc.resample_equal(samp, wts)

#choosing random indicies to draw from properly weighted posterior samples
draws=np.random.randint(0, samples2.shape[0], Nspectra)
Nwno_bins=xsecs[2].shape[0]
y_mod_array=np.zeros((Nwno_bins, Nspectra))
y_binned_array=np.zeros((len(wlgrid), Nspectra))

for i in range(Nspectra):
    print(i)
    #make sure this is the same as in log-Like
    Tirr, logMet, logCtoO, logKzz, fsed ,logPbase,logCldVMR,xRp=samples2[draws[i],:]
    x=np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M, logKzz, fsed,logPbase,logCldVMR, logKcld, logRayAmp, RaySlope])
    print(samples.samples[draws[i],:])
    gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1., 1.])
    y_binned,y_mod,wno,atm=fx(x,wlgrid,gas_scale,xsecs)
    y_mod_array[:,i]=y_mod
    y_binned_array[:,i]=y_binned

#saving these arrays since it takes a few minutes to generate
pickle.dump([wlgrid, y_meas, err, y_binned_array, wno, y_mod_array],open('spectral_samples.pic','wb'))



#PLOTTING SPECTRAL SPREAD-----------------------------------------------------
wlgrid, y_meas, err, y_binned_array, wno, y_mod_array=pickle.load(open('spectral_samples.pic','rb'))
y_median=np.zeros(wno.shape[0])
y_high_1sig=np.zeros(wno.shape[0])
y_high_2sig=np.zeros(wno.shape[0])
y_low_1sig=np.zeros(wno.shape[0])
y_low_2sig=np.zeros(wno.shape[0])

for i in range(wno.shape[0]):
    percentiles=np.percentile(y_mod_array[i,:],[4.55, 15.9, 50, 84.1, 95.45])
    y_low_2sig[i]=percentiles[0]
    y_low_1sig[i]=percentiles[1]
    y_median[i]=percentiles[2]
    y_high_1sig[i]=percentiles[3]
    y_high_2sig[i]=percentiles[4]


from matplotlib.pyplot import *
from matplotlib.ticker import FormatStrFormatter

ymin=np.min(y_meas)*1E2*0.98
ymax=np.max(y_meas)*1E2*1.02
fig1, ax=subplots()
xlabel('$\lambda$ ($\mu$m)',fontsize=14)
ylabel('(R$_{p}$/R$_{*}$)$^{2} \%$',fontsize=14)
minorticks_on()


#for i in range(20): plot(wlgrid, y_binned_array[:,i]*100.,alpha=0.5,color='red')
#for i in range(20): plot(1E4/wno, y_mod_array[:,i]*100.,alpha=0.5,color='red')

fill_between(1E4/wno[::-1],y_low_2sig[::-1]*100,y_high_2sig[::-1]*100,facecolor='r',alpha=0.5,edgecolor='None')
fill_between(1E4/wno[::-1],y_low_1sig[::-1]*100,y_high_1sig[::-1]*100,facecolor='r',alpha=1.,edgecolor='None')


errorbar(wlgrid, y_meas*100, yerr=err*100, xerr=None, fmt='Dk')
plot(1E4/wno, y_median*1E2)
ax.set_xscale('log')
ax.set_xticks([0.3, 0.5,0.8,1,1.4, 2, 3, 4, 5])
ax.axis([0.3,5.0,ymin,ymax])

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(length=5,width=1,labelsize='small',which='major')
savefig('spectrum_fits.pdf',fmt='pdf')
show()
close()


