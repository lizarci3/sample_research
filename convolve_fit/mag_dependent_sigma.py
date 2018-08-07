# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 10:49:37 2017

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import timeit
import math
#from scipy.stats import chisquare
import scipy
from scipy import optimize, interpolate
import astropy
from scipy import integrate
from astropy.convolution import Gaussian1DKernel
import os
import datetime
from extra_functions import hist, mf

start=timeit.time.time()
now = datetime.datetime.now()
    

def find_perturbed_masses(bks, bmsl1, funcsigks):
    
    array = np.c_[bks, bmsl1]
    
    mslpt = []
    
    for i in range(len(array)):
        
        tsig = funcsigks(array[i][0]) #for each ks mag find associated 1sgima 
        
        uerrb = (10**(array[i][1]+tsig))-(10**(array[i][1]))
        lerrb = (10**(array[i][1]))-(10**(array[i][1]-tsig))
        sigg = ((uerrb)+(lerrb))/2
        
        
        mcmass1 = (10**(array[i][1]))+(np.random.normal(loc=0, scale=(sigg), size=1)) #boo random.normal std can only be a float
        
        if(mcmass1!='nan'):        
            mslp = np.log10(mcmass1)
            mslpt.append(mslp)
        
        
        i = i + 1
    
    return(mslpt)
    
    
def convolve_each_point(smfbs, sch, funcsigks):
    
    '''
    I tried, for each point in the SMF, to define a given gaussian kernel for that point (magnitude dependent kernel)
    and then to convolve that point with that kernel
    '''
    schc = []
    
    for i in xrange(len(sch)): ## Convolve Schechter function in at each point in the function with a different gaussian kernel
        gauss1 = gausssian_for_point(smfbs[i][0], sch[i], funcsigks)
        schcp = astropy.convolution.convolve(sch, gauss1, boundary='extend', normalize_kernel=True)
        schc.append(schcp[i])
        
    return(schc)
                    
def gausssian_for_point(massbin, schbin, sigksfunc):
    '''
    Create a gaussian for each specific binned mass, that has a gaussian
    that is dependent on magnitude
    massbin = in dex, mass of this particular bin
    schbin = Schchter value in this point (linearspace)
    sigksfunc = function that obtains sigma for any Ks mag
    '''
    
    ## Define gaussian to convolve
    gx = np.arange(-rangeg, rangeg, (max(smf_or[:,0])-min(smf_or[:,0]))/len(smf_or[:,0]))
    
    ## Find the equivalent magnitude at that point
    ksp = (massbin-bo)/mo 
    
    ## Find the associated 1 sigma at that magnitude
    sigks = sigksfunc(ksp)
    
    #gaussa1 = np.exp((-(gx/sigks)**2)/2)
    #a1 = integrate.trapz(gaussa1, gx) # Normalize, area under curve
    #gauss1 = gaussa1/a1
    
    gauss2 = np.exp(-np.power(gx - 0, 2.) / (2 * np.power(sigks, 2.))) #Essentially the same as above
    
    return(gauss2)
    
################################################################################


def interpolate_sigkms(sigksm):
    
    
    frac = interpolate.interp1d(sigksm[:,0], sigksm[:,1], kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=(5.268905329969642004e-02), assume_sorted=False)  #If below interpolation limit give smallest sigma that we can safely calculate, i.e., Ks~19
    '''
    plt.figure()
    plt.title('Interpolation of Ks vs 1sig')
    plt.scatter(np.arange(17,22,0.1), frac(np.arange(17,22,0.1)), color='b')
    '''
    return(frac)
    
    
################################################################################

flag = 0
print('With or without double cores?')
flag = raw_input()
maglim = 23.5
binsize = 0.2

if(flag=='without'):
    data = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_dw_cd_nostr_maglim_v_mu_2s01z05k_clean.dat') #artifacts have been removed as well as double cores and it starts at 18.5 (to remove artifacts but to be extra sure I added the 18.5 line anyway)
    outfilec=open('test_convolved_magdep2_'+str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'.dat', 'w')

else:
    data = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/smf/250716/eddington/bs+convolved/final/one_sig_magdependent/test_cat_wdoublecores.dat') #Double cores have been left in this version, FOR TESTING PURPOSES ONLY
    outfilec=open('test_convolved_magdep2_'+str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'_wdcores.dat', 'w')

data = data[data[:,6]<maglim]
data = data[data[:,6]>18.5] #already in place in catalog, this is just paranoia really
 

### Set up your parameter space
astart = 0.3
aend = 0.8
astep = 0.01

mstart = 10.54
mend = 10.65
mstep = 0.01

pstart = 0.000215
pend = 0.000248
pstep = 0.000001


rangeg = 21



## Define Ks-Mass conversion

mo=-0.348270925509 #calculated in ks_mass_220716.py
bo=18.2842292843

ks = data[:,6]
kse = data[:,7]
v = data[:,9] #survey volume, Mpc^3 p(z,m)
w = data[:,8]

## Define limits for binning
mlow=(-0.348270925509)*maglim+18.2842292843 #completeness of the sample
mhigh=13
#print('Mlow', mlow, 'Mhigh', mhigh)


## Find masses
msl = (mo)*(ks)+bo

#print('Min Max, Len mass', min(msl), max(msl), len(msl))
  

sigksm = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/ks-mass/220716/ks_1sigma.dat')
funcsigks = interpolate_sigkms(sigksm)

bs = 3 # how many bootstrap re-sampling to do
i = 0
x = np.arange(9,12,0.1)
   
  
smf_or=mf(msl, mlow, mhigh, np.ones(len(msl))*(w/v), width=binsize)
smf_or = smf_or[smf_or[:,1]>0]
 


### Test the output of your SMF
'''
plt.figure(figsize=(10,7), facecolor='white', edgecolor='b')
plt.ylim(-8,-3)
plt.xlim(10, 12.5)
   
if(flag=='without'):
    plt.title('Original Stellar Mass Function PEGs-NOdobulecores', fontsize=25)
else:
    plt.title('Original Stellar Mass Function PEGs-wdobulecores', fontsize=25)
    
plt.xlabel(r'$\mathrm{log(M_*/M_{\odot})}$',fontsize=22,labelpad=10)
plt.ylabel(r'$\Phi\,\mathrm{(dex^{-1}\,Mpc^{-3})}$',fontsize=22)
y_low=np.log10(smf_or[:,1])-np.log10(smf_or[:,1]-smf_or[:,2])
y_high=np.log10(smf_or[:,1]+smf_or[:,2])-np.log10(smf_or[:,1])
plt.errorbar(smf_or[:,0], np.log10(smf_or[:,1]), yerr=(y_low, y_high),color='#e0301e', linestyle='None', marker='o', alpha=0.3, markersize=10,  elinewidth=3, markeredgewidth=1, ecolor='black', label='Deep+Wide')
'''

   

while(i<bs):
    print('Re-sampling #', i)
    
    
    ## Randomly use a sample from the original data
    n = np.random.randint(0,len(ks),(len(ks))) # Take a randomize sample of galaxies from the original catalog        
    
    
    bks=ks[n]
    bv=v[n]
    bw=w[n]
    
    ## Find masses
    #mslp = (mo)*(bks)+bo #to start use un-perturbed masses, to test the Schechter-convolved fit
    bmsl1 = (mo)*(bks)+bo
    
    ### Pertub masses around sigma
    
    
    mslp = find_perturbed_masses(bks, bmsl1, funcsigks)
   
    #### After doing some testing in the jupyter notebook from log to linear space
    #### I realized I do not need to change w or v. The reasoning behind w is obvious
    #### slight mass estimation but object still comes from same field. And for v, the
    #### change in v is motivated for a change in Ks but what we are saying is that we are
    #### perturbing M inside the uncertainty in M not affecting Ks original measurements
    
    ## Bin you data
    smfbs=mf(mslp, mlow, mhigh, np.ones(len(mslp))*(bw/bv), width=binsize)
    smfbs = smfbs[smfbs[:,1]>0]
  
       
    
    ### Test the output of your SMF
    '''
    plt.figure(figsize=(10,7), facecolor='white', edgecolor='b')
    plt.ylim(-8,-3)
    plt.xlim(10, 12.5)
    
    y_lowbs=np.log10(smfbs[:,1])-np.log10(smfbs[:,1]-smfbs[:,2])
    y_highbs=np.log10(smfbs[:,1]+smfbs[:,2])-np.log10(smfbs[:,1])
    
    if(flag=='without'):
        plt.title('BSMC Stellar Mass Function PEGs-NOdoublecores', fontsize=25)
    else:
        plt.title('BSMC Stellar Mass Function PEGs-Wdoublecores', fontsize=25)
        
        
    plt.xlabel(r'$\mathrm{log(M_*/M_{\odot})}$',fontsize=22,labelpad=10)
    plt.ylabel(r'$\Phi\,\mathrm{(dex^{-1}\,Mpc^{-3})}$',fontsize=22)
    plt.errorbar(smfbs[:,0], np.log10(smfbs[:,1]), yerr=(y_lowbs, y_highbs),color='black', linestyle='None', marker='p', alpha=0.5, markersize=10,  elinewidth=3, markeredgewidth=1, ecolor='black', label='BS')
    plt.errorbar(smf_or[:,0], np.log10(smf_or[:,1]), yerr=(y_low, y_high),color='#e0301e', linestyle='None', marker='o', alpha=0.4, markersize=10,  elinewidth=3, markeredgewidth=1, ecolor='black', label='Original')
    '''
    
    ### 'Hand-made' Chi2 Fit, it agrees with Scipy curve fit
    
    # Now we build nested loops to change each parameter
    # Start with the minimum values and add each step until reaching their max value
    # Find min X2 only
    
    alpha = astart
    mstar = mstart
    phistar = pstart
    
    
    
    chitc = 99999999
    altc = 0
    mstc = 0 
    pstc = 0
    
    
   
    
    while (alpha<=aend):
        mstar = mstart
        while (mstar<=mend):
            phistar=pstart
            while (phistar<=pend):
                #print('model', alpha, mstar, phistar)
                schc = []
                ## Define the 'expected' function given the previous point in parameter space
                sch = phistar * np.log(10.) * 10**((smfbs[:,0]-mstar)*(alpha+1)) * np.exp(-(10**((smfbs[:,0]-mstar))))
                
                ## First option- convolve each point with corresponding kernel
                schc = convolve_each_point(smfbs, sch, funcsigks)
                
                ## Second option -  convolve all the Schechter function with an array of gaussians: The number of dimensions of gaussian
                ## should match those for the array, and the dimensions should be odd in all directions.
                
                ## Chi2 computed with this convolved Schechter point in model space
                chic2 = np.sum(((smfbs[:,1]-schc)/(smfbs[:,2]))**2)
                
                if(chic2<chitc):
                    chitc = chic2
                    altc = alpha
                    mstc = mstar
                    pstc = phistar
                    
   
                
                phistar=phistar+pstep
                
            mstar=mstar+mstep
            
        alpha=alpha+astep
    
    
    print (mstc, altc, pstc, chitc, file=outfilec)
    #print (mstc, altc, pstc, chitc)
    
    '''
    x = np.arange(9,12,0.1)        
    schcp = pstc * np.log(10.) * 10**((x-mstc)*(altc+1)) * np.exp(-(10**((x-mstc))))
    plt.plot(x, np.log10(schcp), linestyle=':', color='black', alpha=0.6)
    plt.text(10.2, -6, 'a: '+str(altc)+' ms: '+str(mstc)+' p*: '+str(pstc), color='black')
    '''
    
        
    
    i = i + 1
    
outfilec.close()

 
print('t_total=',timeit.time.time()-start, 's') 
