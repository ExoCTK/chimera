import math
import numpy as np
import scipy as sp
from array import *
from scipy import interpolate
from scipy import signal
from scipy import special
from scipy import interp
from scipy import ndimage
import pdb
from pickle import *
from numba import jit
import time
from scipy.interpolate import RegularGridInterpolator

#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
@jit(nopython=True)
def kcoeff_interp(logPgrid, logTgrid, logPatm, logTatm, wnogrid, kcoeff):
    Ng, NP, NT, Nwno, Nord=kcoeff.shape
    Natm=len(logTatm)
    kcoeff_int=np.zeros((Natm,Nwno,Ng,Nord))

    for i in range(Natm):  #looping through atmospheric layers

        y=logPatm[i]
        x=logTatm[i]

        p_ind_hi=np.where(logPgrid>=y)[0][0]
        p_ind_low=np.where(logPgrid<y)[0][-1]
        T_ind_hi=np.where(logTgrid>=x)[0][0]
        T_ind_low=np.where(logTgrid<x)[0][-1]

        y2=logPgrid[p_ind_hi]
        y1=logPgrid[p_ind_low]
        x2=logTgrid[T_ind_hi]
        x1=logTgrid[T_ind_low]
    
        for j in range(Ng): #looping through gases
            for k in range(Nwno): #looping through wavenumber
                for l in range(Nord): #looping through g-ord
                    arr=kcoeff[j,:,:,k,l]
                    Q11=arr[p_ind_low,T_ind_low]
                    Q12=arr[p_ind_hi,T_ind_low]
                    Q22=arr[p_ind_hi,T_ind_hi]
                    Q21=arr[p_ind_low,T_ind_hi]
                    fxy1=(x2-x)/(x2-x1)*Q11+(x-x1)/(x2-x1)*Q21
                    fxy2=(x2-x)/(x2-x1)*Q12+(x-x1)/(x2-x1)*Q22
                    fxy=(y2-y)/(y2-y1)*fxy1 + (y-y1)/(y2-y1)*fxy2
                    kcoeff_int[i,k,j,l]=fxy

    return kcoeff_int

#########


#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
@jit(nopython=True)
def CalcTauXsecCK(kcoeffs,Z,Pavg,Tavg, Fractions, r0,gord, wts, Fractions_Continuum, xsecContinuum):
    ngas=Fractions.shape[0]
    nlevels=len(Z)
    nwno=kcoeffs.shape[1]
    trans=np.zeros((nwno, nlevels))+1.
    dlarr=np.zeros((nlevels,nlevels))
    ncont=xsecContinuum.shape[-1]
    uarr=np.zeros((nlevels,nlevels))
    kb=1.38E-23
    kbTavg=kb*Tavg
    Pavg_pascal=1E5*Pavg
    for i in range(nlevels-2):
        for j in range(i):
            index=i-j-1
            r1=r0+Z[i]
            r2=r0+Z[i-j]
            r3=r0+Z[index]
            dlarr[i,j]=(r3**2-r1**2)**0.5-(r2**2-r1**2)**0.5
            uarr[i,j]=dlarr[i,j]*Pavg_pascal[index]/kbTavg[index]
    
    for v in range(nwno):
        for i in range(nlevels-2):
            transfull=1.
            #for CK gases--try to do ALL gases as CK b/c of common interpolation
            for k in range(ngas):
                transtmp=0.
                for l in range(len(wts)):
                    tautmp=0.
                    for j in range(i):
                        index=i-j-1
                        tautmp+=2.*Fractions[k,index]*kcoeffs[index,v,k,l]*uarr[i,j]
                    transtmp+=np.exp(-tautmp)*wts[l]/2.
                transfull*=transtmp
            #for continuum aborbers (gas rayligh, condensate scattering etc.--nlayers x nwno x ncont
            #'''
            for k in range(ncont):
                tautmp=0.
                for j in range(i):

                    index=i-j-1
                    tautmp+=2.*Fractions_Continuum[k,index]*xsecContinuum[index,v,k]*uarr[i,j]

                transfull*=np.exp(-tautmp)
            #'''
            trans[v,i]=transfull
    return trans


#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
def tran(xsects, T, P, mmw,Ps,CldOpac,alphaH2O,alphaCH4,alphaCO,alphaCO2,alphaNH3,alphaNa,alphaK,alphaTiO,alphaVO, alphaC2H2, alphaHCN, alphaH2S,alphaFeH,fH2,fHe,amp,power,f_r,M,Rstar,Rp):
    t1=time.time()

    #renaming variables, bbecause why not
    fH2=fH2
    fHe=fHe
    fH2O=alphaH2O
    fCH4=alphaCH4
    fCO=alphaCO
    fCO2=alphaCO2
    fNH3=alphaNH3
    fNa=alphaNa
    fK=alphaK
    fTiO=alphaTiO
    fVO=alphaVO
    fC2H2=alphaC2H2
    fHCN=alphaHCN
    fH2S=alphaH2S
    fFeH=alphaFeH
    mmw=mmw
    #pdb.set_trace()
    #Na and K are fixed in this model but can be made free parameter if desired
    #If T < 800 K set these equal to 0!!!--they condense out below this temperature (roughly)
   
    
    Fractions = np.array([fH2,fHe,fH2O, fCH4, fCO, fCO2, fNH3,fNa,fK,fTiO,fVO,fC2H2,fHCN,fH2S,fFeH])  #gas mole fraction profiles
                        #H2Ray, HeRay  Ray General,
    Frac_Cont = np.array([fH2,fHe,fH2*0.+1.,fH2*0.+1])  #continuum mole fraction profiles
    Frac_Cont=np.concatenate((Frac_Cont, f_r),axis=0)

    #Load measured cross-sectional values and their corresponding
    #T,P,and wno grids on which they were measured
    Pgrid, Tgrid, wno, gord, wts, xsecarr, radius, Mies=xsects[0:8]
    '''
    Pgrid = restore.xsects[0]
    Tgrid = restore.xsects[1]
    wno = restore.xsects[2]
    gord=restore.xsects[3]
    wts=restore.xsects[4]
    xsecarr = restore.xsects[5]
    radius=restore.xsects[6]
    Mies=restore.xsects[7]
    '''
    
    #Calculate Temperature, Pressure and Height grids on which
    #transmissivity will be computed
    n = len(P)
    nv = len(wno)
    
    
    Z=np.zeros(n)  #level altitudes
    dZ=np.zeros(n)  #layer thickness array
    r0=Rp*69911.*1.E3  #converting planet radius to meters
    mmw=mmw*1.660539E-27  #converting mmw to Kg
    kb=1.38E-23
    G=6.67384E-11
    M=M*1.898E27

    
    #Compute avg Temperature at each grid
    Tavg = np.array([0.0]*(n-1))
    Pavg = np.array([0.0]*(n-1))
    for z in range(n-1):
        Pavg[z] = np.sqrt(P[z]*P[z+1])
        Tavg[z] = interp(np.log10(Pavg[z]),sp.log10(P),T)
    #create hydrostatic altitutde grid from P and T
    Phigh=P.compress((P>Ps).flat)  #deeper than reference pressure
    Plow=P.compress((P<=Ps).flat)   #shallower than reference pressure
    for i in range(Phigh.shape[0]):  #looping over levels above ref pressure
        i=i+Plow.shape[0]-1
        g=G*M/(r0+Z[i])**2#g0*(Rp/(Rp+Z[i]/(69911.*1E3)))**2
        H=kb*Tavg[i]/(mmw[i]*g)  #scale height
        dZ[i]=H*np.log(P[i+1]/P[i]) #layer thickness, dZ is negative
        Z[i+1]=Z[i]-dZ[i]   #level altitude
        #print(P[i], H/1000, Z[i]/1000, g)
    for i in range(Plow.shape[0]-1):  #looping over levels below ref pressure
        i=Plow.shape[0]-i-1
        g=G*M/(r0+Z[i])**2#g0*(Rp/(Rp+Z[i]/(69911.*1E3)))**2
        H=kb*Tavg[i]/(mmw[i]*g)
        dZ[i]=H*np.log(P[i+1]/P[i])
        Z[i-1]=Z[i]+dZ[i]
        #print(P[i], H/1000., Z[i]/1000, g)

    #pdb.set_trace()
    #Interpolate values of measured cross-sections at their respective
    #temperatures pressures to the temperature and pressure of the
    #levels on which the optical depth will be computed
    t2=time.time()
    #print('Setup', t2-t1)
    #make sure   200 <T <4000 otherwise off cross section grid
    TT=np.zeros(len(Tavg))
    TT[:]=Tavg
    TT[Tavg < 300] = 300.
    TT[Tavg > 3000] = 3000.
    PP=np.zeros(len(Pavg))
    PP[:]=Pavg
    PP[Pavg < 3E-6]=3E-6
    PP[Pavg >=300 ]=300


    kcoeffs_interp=10**kcoeff_interp(np.log10(Pgrid), np.log10(Tgrid), np.log10(PP), np.log10(TT), wno, xsecarr)
    t3=time.time()
    #print('Kcoeff Interp', t3-t2)
    #continuum opacities (nlayers x nwnobins x ncont)***********
    xsec_cont=kcoeffs_interp[:,:,0,0]
    wave = (1/wno)*1E8
    sigmaH2 = xsec_cont*0.+1*((8.14E-13)*(wave**(-4.))*(1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4  #H2 gas Ray
    sigmaHe = xsec_cont*0.+1*((5.484E-14)*(wave**(-4.))*(1+(2.44E5)*(wave**(-2.))))*1E-4   #He gas Ray
    #Rayleigh Haze from des Etangs 2008
    wno0=1E4/0.43
    sigmaRay=xsec_cont*0.+2.E-27*amp*(wno/wno0)**power*1E-4
    #grey cloud opacity
    sigmaCld=xsec_cont*0.+CldOpac

    #mie scattering 
    xsecMie=Mies[0,0,:,:].T
    sigmaMie=np.repeat(xsecMie[np.newaxis,:,:],len(Pavg),axis=0)

    xsecContinuum=np.array([sigmaH2.T,sigmaHe.T,sigmaRay.T,sigmaCld.T]).T #building continuum xsec array (same order as cont_fracs)
    xsecContinuum=np.concatenate((xsecContinuum, sigmaMie),axis=2)
    #(add more continuum opacities here and in fractions)
    t4=time.time()
    #print("Continuum Xsec Setup ", t4-t3)
    #********************************************
    #Calculate transmissivity as a function of
    #wavenumber and height in the atmosphere
    t=CalcTauXsecCK(kcoeffs_interp,Z,Pavg,Tavg, Fractions, r0,gord,wts,Frac_Cont,xsecContinuum)
    t5=time.time()
    #print('Transmittance', t5-t4)    

    #Compute Integral to get (Rp/Rstar)^2 (equation in brown 2001, or tinetti 2012)
    F=((r0+np.min(Z[:-1]))/(Rstar*6.955E8))**2+2./(Rstar*6.955E8)**2.*np.dot((1.-t),(r0+Z)*dZ)
    t6=time.time()
    #print('Total in Trans', t6-t1)
    #pdb.set_trace()
    return wno, F, Z#, TauOne
#**************************************************************

#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
def instrument_tran_non_uniform(wlgrid,wno, Fp):
    szmod=wlgrid.shape[0]
    delta=np.zeros(szmod)
    Fratio=np.zeros(szmod)
    for i in range(szmod-1):
        delta[i]=wlgrid[i+1]-wlgrid[i]  
    delta[szmod-1]=delta[szmod-2] 

    for i in range(szmod-1):
        i=i+1
        loc=np.where((1E4/wno > wlgrid[i]-0.5*delta[i-1]) & (1E4/wno < wlgrid[i]+0.5*delta[i]))
        Fratio[i]=np.mean(Fp[loc])

    loc=np.where((1E4/wno > wlgrid[0]-0.5*delta[0]) & (1E4/wno < wlgrid[0]+0.5*delta[0]))
    Fratio[0]=np.mean(Fp[loc])
    
    Fratio_int=Fratio
    return Fratio_int, Fp



#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
def xsects_HST(wnomin, wnomax):
    ### Read in CK arrays
    # H2H2
    file='../ABSCOEFF_CK/H2H2_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2H2=10**(kcoeff-4.)
    # H2He
    file='../ABSCOEFF_CK/H2He_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2He=10**(kcoeff-4.)
    # H2O
    file='../ABSCOEFF_CK/H2O_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2O=10**(kcoeff-4.)
    # CH4
    file='../ABSCOEFF_CK/CH4_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrCH4=10**(kcoeff-4.)
    # CO
    file='../ABSCOEFF_CK/CO_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrCO=10**(kcoeff-4.)
    # CO2
    file='../ABSCOEFF_CK/CO2_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrCO2=10**(kcoeff-4.)
    # NH3
    file='../ABSCOEFF_CK/NH3_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrNH3=10.**(kcoeff-4.)
    # Na
    file='../ABSCOEFF_CK/Na_allard_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrNa=10.**(kcoeff-4.)
    # K
    file='../ABSCOEFF_CK/K_allard_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrK=10.**(kcoeff-4.)
    # TiO
    file='../ABSCOEFF_CK/TiO_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrTiO=10.**(kcoeff-4.)
    # VO
    file='../ABSCOEFF_CK/VO_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrVO=10.**(kcoeff-4.)
    # C2H2
    file='../ABSCOEFF_CK/C2H2_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrC2H2=10.**(kcoeff-4.)
    # HCN
    file='../ABSCOEFF_CK/HCN_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrHCN=10.**(kcoeff-4.)
    # H2S
    file='../ABSCOEFF_CK/H2S_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2S=10.**(kcoeff-4.)
    # FeH
    file='../ABSCOEFF_CK/FeH_CK_STIS_WFC3_10gp_2000_30000wno.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrFeH=10.**(kcoeff-4.)
    
    xsecarr = np.log10(np.array([xsecarrH2H2,xsecarrH2He, xsecarrH2O, xsecarrCH4, xsecarrCO,xsecarrCO2,xsecarrNH3,xsecarrNa, xsecarrK, xsecarrTiO, xsecarrVO, xsecarrC2H2,xsecarrHCN,xsecarrH2S,xsecarrFeH]))
    
    #loading mie coefficients-----------------------------
    cond_name='MgSiO3'
    file='../MIE_COEFFS/'+cond_name+'_r_0.01_300um_wl_0.3_200um_interp_STIS_WFC3_2000_30000wno.pic'
    wno,radius,Mies = load(open(file,'rb'), encoding='latin1')  #radius is in mixrons!
    SSA=Mies[1,:,:]/Mies[0,:,:]#single scatter albedo
    Mies[1,:,:]=SSA  #single scatter albedo
    Mg2SiO4=Mies  #Mies = Qext, Qs, asym
    xxsec=Mg2SiO4[0,:,:].T*np.pi*radius**2*1E-12 #scattering cross-section
    Mg2SiO4[0,:,:]=xxsec.T

    mies_arr=np.array([Mg2SiO4])

    #cropping wavenumber range
    #wnomin=6000#2000#5000#2000
    #wnomax=9400#30000#10000#30000
    loc=np.where((wno <= wnomax) & (wno >= wnomin))[0]
    wno=wno[loc]
    xsecarr=xsecarr[:,:,:,loc,:]
    mies_arr=mies_arr[:,:,:,loc]    


    #loading chemistry grid
    logCtoO, logMet, Tarr, logParr, gases=load(open('chem_full.pic','rb'), encoding='latin1')
    print('Cross-sections Loaded')
    return P,T,wno,g,wts,xsecarr,radius*1E-6,mies_arr,logCtoO, logMet, Tarr, logParr, np.log10(gases)



#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
def xsects_JWST(wnomin, wnomax):
    ### Read in CK arrays
    # H2H2
    file='../ABSCOEFF_CK/CK_H2H2_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2H2=10**(kcoeff-4.)
    # H2He
    file='../ABSCOEFF_CK/CK_H2He_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2He=10**(kcoeff-4.)
    # H2O
    file='../ABSCOEFF_CK/CK_H2O_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2O=10**(kcoeff-4.)
    # CH4
    file='../ABSCOEFF_CK/CK_CH4_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrCH4=10**(kcoeff-4.)
    # CO
    file='../ABSCOEFF_CK/CK_CO_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrCO=10**(kcoeff-4.)
    # CO2
    file='../ABSCOEFF_CK/CK_CO2_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrCO2=10**(kcoeff-4.)
    # NH3
    file='../ABSCOEFF_CK/CK_NH3_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrNH3=10.**(kcoeff-4.)
    # Na
    file='../ABSCOEFF_CK/CK_Na_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrNa=10.**(kcoeff-4.)
    # K
    file='../ABSCOEFF_CK/CK_K_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrK=10.**(kcoeff-4.)
    # TiO
    file='../ABSCOEFF_CK/CK_TiO_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrTiO=10.**(kcoeff-4.)
    # VO
    file='../ABSCOEFF_CK/CK_VO_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrVO=10.**(kcoeff-4.)
    # C2H2
    file='../ABSCOEFF_CK/CK_C2H2_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrC2H2=10.**(kcoeff-4.)
    # HCN
    file='../ABSCOEFF_CK/CK_HCN_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrHCN=10.**(kcoeff-4.)
    # H2S
    file='../ABSCOEFF_CK/CK_H2S_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrH2S=10.**(kcoeff-4.)
    # FeH
    file='../ABSCOEFF_CK/CK_FeH_JWST_R100.pic'
    P,T,wno,kcoeff,g,wts = load(open(file,'rb'), encoding='latin1')
    xsecarrFeH=10.**(kcoeff-4.)
    
    xsecarr = np.log10(np.array([xsecarrH2H2,xsecarrH2He, xsecarrH2O, xsecarrCH4, xsecarrCO,xsecarrCO2,xsecarrNH3,xsecarrNa, xsecarrK, xsecarrTiO, xsecarrVO, xsecarrC2H2,xsecarrHCN,xsecarrH2S,xsecarrFeH]))
    
    #loading mie coefficients-----------------------------
    cond_name='MgSiO3'
    file='../MIE_COEFFS/'+cond_name+'_r_0.01_300um_wl_0.3_200um_interp_JWST_R100.pic'
    wno,radius,Mies = load(open(file,'rb'), encoding='latin1')  #radius is in mixrons!
    SSA=Mies[1,:,:]/Mies[0,:,:]#single scatter albedo
    Mies[1,:,:]=SSA  #single scatter albedo
    Mg2SiO4=Mies  #Mies = Qext, Qs, asym
    xxsec=Mg2SiO4[0,:,:].T*np.pi*radius**2*1E-12 #scattering cross-section
    Mg2SiO4[0,:,:]=xxsec.T

    mies_arr=np.array([Mg2SiO4])

    #cropping wavenumber range
    #wnomin=6000#2000#5000#2000
    #wnomax=9400#30000#10000#30000
    loc=np.where((wno <= wnomax) & (wno >= wnomin))[0]
    wno=wno[loc]
    xsecarr=xsecarr[:,:,:,loc,:]
    mies_arr=mies_arr[:,:,:,loc]    


    #loading chemistry grid
    logCtoO, logMet, Tarr, logParr, gases=load(open('chem_full.pic','rb'), encoding='latin1')
    print('Cross-sections Loaded')
    return P,T,wno,g,wts,xsecarr,radius*1E-6,mies_arr,logCtoO, logMet, Tarr, logParr, np.log10(gases)


#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
@jit(nopython=True)
def cloud_profile(fsed,cloud_VMR, Pavg, Pbase):
    cond=cloud_VMR
    loc0=np.where(Pbase >= Pavg)[0][-1]
    cond_mix=np.zeros(len(Pavg))+1E-50
    cond_mix[0:loc0+1]=cond*(Pavg[0:loc0+1]/Pavg[loc0])**fsed  #A&M2001 eq 7., but in P-coordinates (using hydrostatic) and definition of f_sed
    return cond_mix

#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
#computes verticle distribution of particles assuming a log normal distribution
#and the balance between uplift and sedimentation
@jit(nopython=True)
def particle_radius(fsed,Kzz,mmw,Tavg, Pavg,g, rho_c,mmw_c, qc,rr):
    dlnr=np.abs(np.log(rr[1])-np.log(rr[0]))
    kb=1.38E-23  #boltzman constant
    mu0=1.66E-27  #a.m.u.
    d=2.827E-10  #bath gas molecule diameter (m)
    alpha=1.4  #alpha factor from A&M 2001 (don't need to change this)
    sig_eff=2  #log-normal particle size distribution width
    
    #atmosphere properties
    H=kb*Tavg/(mmw*mu0*g)  #scale height
    rho_a=Pavg*mmw*mu0*1E5/(kb*Tavg)  #atmospheric mass density
    
    wmix=Kzz/H  #vertical mixing velocity
    mfp=kb*Tavg/(2**0.5*np.pi*d**2*Pavg*1E5)   #mean free path
    eta=5./16.*np.sqrt(np.pi*2.3*mu0*kb*Tavg)*(Tavg/59.7)**.16/(1.22*np.pi*d**2) #dynamic viscosity of bath gas
    
    #computing varius radius profiles
    r_sed=2./3.*mfp*((1.+10.125*eta*wmix*fsed/(g*(rho_c-rho_a)*mfp**2))**.5-1.)  #sedimentation radius
    #r_eff=r_sed*np.exp(-0.5*(alpha+1)*np.log(1.+sig_eff))  #charnay formula--"effective" radius of particle distribution
    r_eff=r_sed*fsed**(1./alpha)*np.exp(-0.5*(alpha+1)*np.log(sig_eff)**2)  #A&M2011 equation 17 effective radius
    r_g=r_sed*fsed**(1./alpha)*np.exp(-0.5*(alpha+6.)*np.log(sig_eff)**2) #A&M formula (13)--lognormal mean (USE THIS FOR RAD)
    
    #droplet VMR
    f_drop=3.*mmw_c*mu0*qc/(4.*np.pi*rho_c*r_g**3)*np.exp(-4.5*np.log(sig_eff)**2)  #
    prob_lnr=np.zeros((len(rr),len(r_g)))
    for i in range(len(prob_lnr)): prob_lnr[i,:]=1./((2.*np.pi)**0.5*np.log(sig_eff))*np.exp(-0.5*np.log(rr[i]/r_g)**2/np.log(sig_eff)**2)*dlnr
    f_r=prob_lnr*f_drop
    return r_sed, r_eff, r_g, f_r


'''
#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
class restore():
    xsects = xsects()

'''
#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
def TP(Teq, Teeff, g00, kv1, kv2, kth, alpha):
    
    
    Teff = Teeff
    f = 1.0  # solar re-radiation factor
    A = 0.0  # planetary albedo
    g0 = g00
    
    # Compute equilibrium temperature and set up gamma's
    T0 = Teq
    gamma1 = kv1/kth
    gamma2 = kv2/kth
    
    # Initialize arrays
    logtau =np.arange(-10,20,.1)
    tau =10**logtau
    
    
    #computing temperature
    T4ir = 0.75*(Teff**(4.))*(tau+(2.0/3.0))
    f1 = 2.0/3.0 + 2.0/(3.0*gamma1)*(1.+(gamma1*tau/2.0-1.0)*sp.exp(-gamma1*tau))+2.0*gamma1/3.0*(1.0-tau**2.0/2.0)*special.expn(2.0,gamma1*tau)
    f2 = 2.0/3.0 + 2.0/(3.0*gamma2)*(1.+(gamma2*tau/2.0-1.0)*sp.exp(-gamma2*tau))+2.0*gamma2/3.0*(1.0-tau**2.0/2.0)*special.expn(2.0,gamma2*tau)
    T4v1=f*0.75*T0**4.0*(1.0-alpha)*f1
    T4v2=f*0.75*T0**4.0*alpha*f2
    T=(T4ir+T4v1+T4v2)**(0.25)
    P=tau*g0/(kth*0.1)/1.E5
    
    
    # Return TP profile
    return T, P


#*******************************************************************
# FILE: xsects.py
#
# DESCRIPTION: 
#
# USAGE: 
#
# RETURNS: 
#*******************************************************************
def fx(x,wlgrid,gas_scale, xsects):
    #print(x)
    p1=time.time()
   
    #Unpacking Guillot 2010 TP profile params (3 params)
    Tiso=x[0]
    logKir=x[1]
    logg1=x[2]
    #Unpacking Chemistry Parms
    Met=10.**x[3]  #metallicity
    CtoO=10.**x[4] #C/O
    logPQC=x[5]  #carbon quench pressure
    logPQN=x[6]  #nitrogen quench pressure
    #unpacking planet params
    Rp=x[7]  #planet radius (in jupiter)
    Rstar=x[8]   #stellar radius (in solar)
    M=x[9]   #planet mass (in jupiter)
    #unpacking and converting A&M cloud params
    Kzz=10**x[10]*1E-4  #Kzz for A&M cloud
    fsed=x[11]  #sedimentation factor for A&M cloud
    Pbase=10.**x[12]  #cloud top pressure
    Cld_VMR=10**x[13]  #Cloud Base Condensate Mixing ratio
    #unpacking and converting simple cloud params
    CldOpac=10**x[14]
    RayAmp=10**x[15]
    RaySlp=x[16]

    #Setting up atmosphere grid****************************************
    logP = np.arange(-6.8,1.5,0.1)+0.1
    P = 10.0**logP
    g0=6.67384E-11*M*1.898E27/(Rp*69911.*1.E3)**2
    kv=10.**(logg1+logKir)
    kth=10.**logKir
    tp=TP(Tiso, 100,g0 , kv, kv, kth, 0.5)
    T = interp(logP,np.log10(tp[1]),tp[0])
    t1=time.time()

    #interpolation chem
    logCtoO, logMet, Tarr, logParr, loggas=xsects[8:]
    Ngas=loggas.shape[-2]
    gas=np.zeros((Ngas,len(P)))
    #capping T at bounds
    TT=np.zeros(len(T))
    TT[:]=T[:]
    TT[TT>2800]=2800
    TT[TT<400]=400
    for j in range(Ngas):
        gas_to_interp=loggas[:,:,:,j,:]
        IF=RegularGridInterpolator((logCtoO, logMet, np.log10(Tarr),logParr),gas_to_interp,bounds_error=False)
        for i in range(len(P)):
            gas[j,i]=10**IF(np.array([np.log10(CtoO), np.log10(Met), np.log10(TT[i]), np.log10(P[i])]))*gas_scale[j]

    #H2O, CH4, CO, CO2, NH3, N2, HCN, H2S,PH3, C2H2, C2H6, Na, K, TiO, VO, FeH, H,H2, He, mmw=gas
    H2Oarr, CH4arr, COarr, CO2arr, NH3arr, N2arr, HCNarr, H2Sarr,PH3arr, C2H2arr, C2H6arr, Naarr, Karr, TiOarr, VOarr, FeHarr, Harr,H2arr, Hearr, mmw=gas
    #Super simplified non-self consistent quenching based on quench pressure
    #Carbon
    PQC=10.**logPQC
    loc=np.where(P <= PQC)
    CH4arr[loc]=CH4arr[loc][-1]
    COarr[loc]=COarr[loc][-1]
    H2Oarr[loc]=H2Oarr[loc][-1]
    CO2arr[loc]=CO2arr[loc][-1]

    #Nitrogen
    PQN=10.**logPQN
    loc=np.where(P <= PQN)
    NH3arr[loc]=NH3arr[loc][-1]
    N2arr[loc]=N2arr[loc][-1]
    t2=time.time()

    #hacked rainout (but all rainout is...)....if a mixing ratio profile hits '0' (1E-12) set it to 1E-20 at all layers above that layer
    rain_val=1E-11
    loc=np.where(TiOarr <= rain_val)[0]
    if len(loc>1): TiOarr[0:loc[-1]-1]=1E-20
    #loc=np.where(VOarr <= rain_val)[0]
    if len(loc>1):VOarr[0:loc[-1]-1]=1E-20 #VO and TiO rainout togather
    loc=np.where(Naarr <= rain_val)[0]
    if len(loc>1): Naarr[0:loc[-1]-1]=1E-20
    loc=np.where(Karr <= rain_val)[0]
    if len(loc>1):Karr[0:loc[-1]-1]=1E-20
    loc=np.where(FeHarr <= rain_val)[0]
    if len(loc>1):FeHarr[0:loc[-1]-1]=1E-20

    #print('Chemistry', t2-t1)

    #ackerman & Marley cloud model here
    mmw_cond=140.961#74.5#140.691
    rho_cond=3270#1984#3270.
    rr=10**(np.arange(-2,2.6,0.1))  #Droplet radii to compute on: MUST BE SAME AS MIE COEFF ARRAYS!!!!!!!!! iF YOU CHANGE THIS IT WILL BREAK
    qc=cloud_profile(fsed,Cld_VMR, P,Pbase)
    r_sed, r_eff, r_g, f_r=particle_radius(fsed,Kzz,mmw,T, P,g0, rho_cond,mmw_cond,qc, rr*1E-6)
   

    Pref=1.1#10.1  #reference pressure bar-keep fixed
    #computing transmission spectrum-----------
    #RayAmp=1E-20
    #RaySlp=0.
    #CldOpac=1E-50

    spec = tran(xsects, T,P,mmw, Pref,CldOpac, H2Oarr, CH4arr,COarr,CO2arr,NH3arr,Naarr,Karr,TiOarr,VOarr,C2H2arr,HCNarr,H2Sarr,FeHarr,H2arr,Hearr,RayAmp,RaySlp,f_r, M, Rstar, Rp)
    wno = spec[0]
    F = spec[1]
    
    y_binned,junk=instrument_tran_non_uniform(wlgrid,wno, F)
    p3=time.time()
    #print("Total ", p3-p1)    
    chemarr=np.array([P,T, H2Oarr, CH4arr,COarr,CO2arr,NH3arr,Naarr,Karr,TiOarr,VOarr,C2H2arr,HCNarr,H2Sarr,FeHarr,H2arr,Hearr,qc,r_eff,f_r])

    return y_binned,F,wno,chemarr


