import numpy as np
from matplotlib import pyplot as plt
import pdb
import dynesty
from multiprocessing import Pool
from fm import *
import pickle
#load crosssections between wnomin and wnomax
xsects=xsects(6000,9400)


# log-likelihood
def loglike(theta):

    #setting default parameters---will be fixed to these values unless replaced with 'theta'
    #planet/star system params--typically not free parameters in retrieval
    Rp= 1.10#0.930#*x[4]# Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)
    Rstar=1.20#0.598   #Stellar Radius in Solar Radii
    M =1.0#1.78    #Mass in Jupiter Masses
    xRp=1.0
    #TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
    Tirr=1300#1500#x[0]#544.54 #terminator **isothermal** temperature--if full redistribution this is equilibrium temp
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


    #unpacking parameters to retrieve
    Tirr, logMet, logCtoO, logKzz, fsed ,logPbase,logCldVMR,xRp=theta
   

    ##all values required by forward model go here--even if they are fixed
    x=np.array([Tirr, logKir,logg1, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M, logKzz, fsed,logPbase,logCldVMR, logKcld, logRayAmp, RaySlope])
    gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1., 1.])
    y_binned,y_mod,wnocrop,atm=fx(x,wlgrid,gas_scale,xsects)


    #y_mod=y_meas+x[1]
    loglikelihood=-0.5*np.sum((y_meas-y_binned)**2/err**2)
    return loglikelihood


# prior transform
def prior_transform(utheta):
    Tirr, logMet, logCtoO, logKzz, fsed ,logPbase,logCldVMR,xRp=utheta
    #prior ranges...
    Tirr = 2600 * Tirr + 400
    logMet= 4.5*logMet-1.5
    logCtoO=1.3*logCtoO-1
    logKzz=6*logKzz+5
    fsed=3.5*fsed+0.5
    logPbase=7.5*logPbase-6.0
    logCldVMR=8*logCldVMR-10
    xRp=1*xRp+0.5

    return Tirr, logMet, logCtoO, logKzz, fsed ,logPbase,logCldVMR,xRp



#####loading in data##########
wlgrid, y_meas, err=pickle.load(open('data.pic','rb'))
outname='dyn_output_500LP.pic'  #dynesty output file name (saved as a pickle)
Nparam=8  #number of parameters--make sure it is the same as what is in prior and loglike
Nproc=4  #number of processors for multi processing--best if you can run on a 12 core+ node or something
Nlive=500 #number of nested sampling live points

pool = Pool(processes=Nproc)
dsampler = dynesty.NestedSampler(loglike, prior_transform, ndim=Nparam,
                                 bound='multi', sample='auto', nlive=Nlive,
                                 update_interval=3., pool=pool, queue_size=Nproc)

#extracting results from sampler object
dres = dsampler.results
#dumping as a pickle
pickle.dump(dres,open(outname,'wb'))  

pdb.set_trace()



