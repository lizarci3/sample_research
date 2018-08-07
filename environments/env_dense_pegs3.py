# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 15:46:46 2016

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import math
from scipy import interpolate
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from extra_functions_heat import hist2


'''
Script created to evaluate the environment of bright passive galaxies-luminosity functions instead of mass functions
'''
#######################################################################
 # Auxiliary functions
 
def make_region_vecindario(array, data_full, rad):
    
    
    # ID that bright galaxy's environment
    for item in array: #Make a region file of found neighboours
        if(item[2]>0): #only neighbours, not th central
            print ('point('+str(data_full[item[3]][4])+','+str(data_full[item[3]][5])+') #point= circle color= yellow', file=reg)
            print ('text('+str(data_full[item[3]][4])+','+str(data_full[item[3]][5]-0.002)+') #text={nb '+str(item[1])+'} color=yellow', file=reg)
            
        elif(item[2]==0): #If distance to central =0, then this is the central
        
            print ('circle('+str(data_full[item[3]][4])+','+str(data_full[item[3]][5])+', '+str(rad)+') #point= circle color= magenta', file=reg)
        
#######################################################################
 # Auxiliary functions
 
def make_region_fullcat(fullcatalog, radius):
    
    total = open('test_full.reg', 'w')
    
    print('# Region file format: DS9 version 4.1',file=total)
    print('# Environment around PEGs galaxies (Environment is comprised of SF and PE gzK galaxies', file=total)
    print('global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',file=total)
    print('fk5',file=total)
    
    fullcatalog = fullcatalog[fullcatalog[:,6]>artifact]
    fullcatalog = fullcatalog[(fullcatalog[:,12]<-21.63)]
    #fullcat = fullcat[fullcat[:,6]>artifact] #Remove bright artifacts    
    #fullcat = fullcat[(fullcat[:,12]<-21.63)] # Cut the Rest-frame Z mags to the completeness of the Passive Sample (Due to the diff redshifts btw PEGs and SFGs the Ks<20.5 was not enough, when transformed to Z rest frame, the SF population has a brighter mag completeness at 20.5 than the PEGs->See chapter 3 for details on the redz distribution)    
    
    
    #fullcatalog = fullcatalog[fullcatalog[:,0]==13]
    
    for item in fullcatalog:
    
        if(item[8]==0):
            coll = 'pink'
        elif(item[8]==1):
            coll = 'cyan'
            
        print ('point('+str(item[4])+','+str(item[5])+') #point= circle color= '+coll, file=total)
        print ('text('+str(item[4])+','+str(item[5]+0.002)+') #text={all '+str(item[1])+'} color='+coll, file=total) 
        print ('circle('+str(item[4])+','+str(item[5])+', '+str(radius/10)+') #point= circle color= '+coll, file=total)
        
    total.close()
        
#################################################################################################### 
            
def make_region_allpeg(pewd):
    
    regall = open('test_all.reg', 'w')


    print('# Region file format: DS9 version 4.1',file=regall)
    print('# Environment around PEGs galaxies (Environment is comprised of SF and PE gzK galaxies', file=regall)
    print('global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',file=regall)
    print('fk5',file=regall)
    
    
    for item in pewd: #Make a region file of found neighboours
        kind = item[8]        
        
        if(kind==0):
            print ('point('+str(item[4])+','+str(item[5])+') #point= circle color= pink', file=regall)
            print ('text('+str(item[4]-0.004)+','+str(item[5])+') #text={pe '+str(item[1])+'} color=pink', file=regall)
        elif(kind==1):
            print ('point('+str(item[4])+','+str(item[5])+') #point= circle color= blue', file=regall)
            print ('text('+str(item[4]-0.004)+','+str(item[5])+') #text={sf '+str(item[1])+'} color=blue', file=regall)
        
    regall.close()
       
####################################################################################################

def make_region_centrals(centrals, maglim, rad):
    
    regc = open('test_centrals.reg','w')
    
    print('# Region file format: DS9 version 4.1',file=regc)
    print('# Central PEGs Brighter than Ks<'+str(maglim), file=regc)
    print('global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',file=regc)
    print('fk5',file=regc)
    
    for cent in central:
    
        print('point('+str(cent[4])+','+str(cent[5])+') #point= circle color= cyan', file=regc) #ID our bright galaxy
        print ('text('+str(cent[4])+','+str(cent[5]+0.002)+') #text={c_'+str(cent[1])+'} color=cyan', file=regc)
            
        print('circle('+str(cent[4])+','+str(cent[5])+','+str(rad)+') #color= green', file=regc) #Radius of search Established with mrad converted to degrees
    
    regc.close()
    
####################################################################################################
    
def  king_of_the_hill(array, fullcat):
    
    '''
    This subroutine should check if there is a more massive companion within the radius of search
    '''
    '''
    At this point array contains the central, and any neighbour inside a radius or mrad
    Each row in array represents an ith neighbour
    Col 0: Field of central
    Col 1: ID of central
    Col 2: Distance to central
    Col 3: Index of neighbour in full cat
    '''
    king = array[array[:,2]==0] #ID the one that was defined as central, i.e., the one that is at a distance 0 from the original central
    #print('Central is', king[:,0], king[:,1], king[:,2], king[:,3])

    ## Find the Ks magnitude of that central
    kking = fullcat[float(king[:][:,3])][6]
    #print('K mag of central,', kking)
    kingf = True # Start assuming that for this particular central, he is the king of his hill
    
    
    # Search within all found neighbours if any of those neighbours has a Ks magnitude that is less than 
    # that of the central. If one is found, the flag will return False, and this central will not
    # be considered since it is not the KING of its HILL
    
    for thing in array:
        
       ks =  fullcat[float(thing[3])][6]
       kind = fullcat[float(thing[3])][8]
       
       if(ks<kking): ## No matter if it is SF or PE if there is a more massive galaxy you are not king of your hill
           #print('KS nb less thank king', kking, ks, kind)
           kingf = False
        
    #print(kingf)
    return(kingf)
        
    
####################################################################################################
    
def test_kdtree(maglim, mrad, fullcat, centrals):
    
    test = open('test_kdtree.dat', 'w')
    
    fullcat_test = fullcat[fullcat[:,0]==13]
    centrals_test = centrals[centrals[:,0]==13]
    
    print('Testing Tree')
    print('I am testing', fullcat_test[:,0], 'with centrals in', centrals_test[:,0])
    
    rad = mrad*1000/8.471/60/60 # from Mpc to deg at z~1.6
    
    
    tree_test = cKDTree(zip(fullcat_test[:,4],fullcat_test[:,5]))
    
    for each in centrals_test:
        print('Test Central', each[0])
        something = tree_test.query_ball_point((each[4], each[5]), r=rad)
        
        for el in something:
            print(fullcat_test[el][0], fullcat_test[el][1], fullcat_test[el][2], fullcat_test[el][3], fullcat_test[el][4], fullcat_test[el][5], fullcat_test[el][6], fullcat_test[el][7], fullcat_test[el][8], fullcat_test[el][9], fullcat_test[el][10], fullcat_test[el][11], fullcat_test[el][12], file=test)
        
        print(something)
        
    test.close()
        
####################################################################################################

def create_cat_nb(maglim, mrad, fullcat, centrals):
    
    
    '''
    In this subroutine we will use a KDTree to search in a given radius
    all neighbours around our centrals. These neighbours are both SF + PEGs
    to the limiting mag of the sample
    '''
    
    ###Transform the radius of search from MPC to deg
    #0.01639 for 0.5 Mpc (example)
    rad = mrad*1000/8.471/60/60 # from Mpc to deg at z~1.6
    print('A radius of ', mrad, ' Mpc gives ', rad, ' in degrees')
    
    
    print('Number of galaxies brighter than ', maglim, ' i.e., # of centrals: ', len(centrals))
    
    ## If you want to see as region files the recovered neighbours uncomment the next lines.
    #make_region_allpeg(fullcat) #For testing purposes
    #make_region_centrals(centrals, maglim, rad) #For testing purposes
    make_region_fullcat(fullcat, rad)
    
    vecindad8 = []
    vearray8 = []
    vearray8u = []
    
    
    ## Make a KDTree of the full catalog, used to search the neighbours
    tree = cKDTree(zip(fullcat[:,4],fullcat[:,5]))## Make a tree of all the available data, not just the  bright objects
    i = 0 #each group will have an ID number to link galaxies from the same group
    
    
    ## For each central, look for galaxies around it, we will look for up to 50 nbs but on avereage there are 3 on the wide, 10 on deep
    
    for item in centrals:
        
        print('Working with central ', item[0], item[1])
        
        
        array = []
        arrayi = []
        arrayc = []
        
        dist, ind = tree.query((item[4], item[5]), k=50, distance_upper_bound=rad)
        
        ### Make an array that contains all the neighbours, each line in the array represents a neighbour. Rmember that in this way
        ### you are also recovering the central and it will have a distance of zero to the central (duh)
        ### 0 Field of central
        ### 1 ID of central
        ### 2 Distance to central
        ### 3 Index in full catalog of this neighbour
        
        arrayi = np.c_[item[0]*np.ones(len(dist)),item[1]*np.ones(len(dist)), dist, ind] # Field, ID, of central, with distance and index of ith neighbour
        '''
        iwant = 35568
        pepitof = arrayi[arrayi[:,0]==13]
        #print(pepitof)
        
        for pep in pepitof:
            ##print(pep)
            #print(fullcat[pep[3]][1])
            pidd = fullcat[pep[3]][1]
            #print(pidd)
            
            if(pidd==iwant):
                print('I am here with dist', pep[2])
        '''
        arrayc = arrayi[arrayi[:,2]==0] # Central galaxy, this is used to check the presence of centrals in the final cat
        array = arrayi [arrayi[:,2]<np.inf] # remove objects outside of your radius: For some reason tree.query flags them as having a distance of inf
        
        # This next line should help you decide how many are 'loners'
        print('For central ', item[0], item[1], ' I found # nbs ', len(array[:,2])) #How many objects were left in the array once you removed objects outside the radius of search    
        
        #Check that you are in fact recovering the central as well. There is an if statement next as well to do it a bit more efficiently, this was used for initial testing
        #print('Central recovered ', arrayc[:,0], arrayc[:,1], arrayc[:,2])
        
        #At this point array contains the central, and any neighbour inside a radius or mrad
        #Each row in array represents an ith neighbour
        #Col 0: Field of central
        #Col 1: ID of central
        #Col 2: Distance to central
        #Col 3: Index of neighbour in full cat
        
        
        ### Make sure you are actually recovering all centrals, so we are comparing the central we are working with with the one recovered in the KDTree with a distance of zero, what will happen if array c is empty?
        if(item[1]!=arrayc[:,1]): ## With this I confirm that all centrals were also recovered
            print('Check me', item[0], item[1])
               
        ## We  are keeping loner centrals 06/02/2017        
        if(len(array)>0): #If you found any neighbours around the galaxy set to 1, if you also want to include centrals without companions set to 0
        
            kingdom = king_of_the_hill(array, fullcat) #make sure there are no other more massive galaxies part of this environment, hence this galaxy is not trully a central
            
            
            if(kingdom==True): #If in fact there are no more massive galaxies in its surroundings
                
                make_region_vecindario(array, fullcat, rad)
                
                for el in array: #save into a vector every index that was classified as a neighbour or central itself (distance=0)
                
                    vec =[]
                    central = -99 #if something goes wrong th -99 flag will show up
                    vecindad8 = np.append(vecindad8, el[3]) #Append only the indices of everyone in array: centrals + neighbours
                    
                    if(el[2]==0): # If the distance to central is zero, this is in fact the central galaxy
                        central = 1 #True
                    elif(el[2]!=0): # if not zero, is a neighbouring galaxy
                        central = 0 #False
                        
                    vec = np.hstack((el, central, i))# for each element found around each central, save it into a vector that will include de info of whether or not this is a central or a nb and which group number this is (stored in the i variable)
                    #print('vec', vec)
                    
                    if(len(vearray8)==0):
                        vearray8 = vec
                    else:
                        vearray8 = np.vstack((vearray8, vec))
                i = i+1
    
          
    vecinos8 = np.unique(vecindad8) #to avoid double-counting neighbouring galaxies, we check that each index is repeated once (An index represents its position in the full catalog, eg. index = 0 represents the first object in the catalog)
    
    #vearray8u = vearray8[vearray8[:,3]==np.array(vecinos8)] #Tried doing an easy way like this but didn't work
    
    
    ## Check that in the full array you are not double-counting nbs by using vecinos8
    
    for ea8 in vecinos8:
        
        for eaa in vearray8:
            
            if(ea8==eaa[3]):
            
                if(len(vearray8u)==0):
                    vearray8u = eaa
                else:
                    vearray8u = np.vstack((vearray8u,eaa))
                
                break #as soon as you find a nb break the cycle so that you don't save more than one
                    
            
    
    ####### At this point vearray8u should contain centrals and its nbs without double counting, whereas vearray8 contains everyone, even those nbs in multiple environments
    ####### See comment below
    ############################ Make full cats out of array
    finala = []
    
    '''
    IMPORTANT:
    For the following loop you can either use vearray8 which includes all found neighbours, or
    you can use vearray8u which makes sure you are not doubly counting
    During last meeting we decided that vearray8 is the way to go
    '''
    
    for eacha in vearray8: # Make full catalog of all neighbours around the central
    
        veca = []
        veca = np.r_[fullcat[eacha[3]][0], fullcat[eacha[3]][1], fullcat[eacha[3]][2], fullcat[eacha[3]][3], fullcat[eacha[3]][4], fullcat[eacha[3]][5], fullcat[eacha[3]][6], fullcat[eacha[3]][7], fullcat[eacha[3]][8], fullcat[eacha[3]][9], fullcat[eacha[3]][11], fullcat[eacha[3]][12], eacha[4], eacha[5]]    
        
        if(len(finala)==0):
            finala = veca
        else:
            finala = np.vstack((finala, veca))
    
    areasa = determine_areas(fullcat, centrals, maglim, rad, mrad, finala)
    finalta = np.c_[finala, areasa*np.ones(len(finala[:,0]))]
    
    
    #### Columns in final array finalta so far
    # (0) FIELD
    # (1) ID
    # (2) X
    # (3) Y
    # (4) RA
    # (5) DEC
    # (6) Ks Mag
    # (7) Ks Mag err
    # (8) Kind 0 = PE, 1 = SF
    # (9) Completeness
    # (10) Redshift (From Muzzin's rainbow plot)
    # (11) Rest-frame Z mag based on given redshift
    # (12) Flag that will ID if the object is a central or not 1 = central 0 = neighbour
    # (13) Group #
    # (14) Areas Mpc^2
    
    cleanc = finalta[finalta[:,12]==1]
    cleannb = finalta[finalta[:,12]==0]
    print('Number of *truly* centrals found,', len(cleanc))
    print('Number of nb around centrals', len(cleannb))
    
    #### Save into a test file the found array. This test file is not in the right folder to avoid over-writing files, so once you are sure this is right move it to the right folder with the right name
    
    #np.savetxt('test2_'+work_with+'_nb_sfpeg_ks'+str(maglim)+'_r'+str(mrad)+'.dat', finalta, fmt='%s', header='Combined catalog of SF+PE galaxies in a neighbouring '+str(mrad)+' Mpc radius of PE galaxies brighter than '+str(maglim)+' in '+work_with+' fields, calculated with env_dense_pegs3.py \n (0) Field \n (1) ID  \n (2) X [pix] \n (3) Y [pix] \n (4) RA [deg] \n (5) DEC [deg] \n (6) Ktot [AB mag] \n (7) Ktot err [AB mag] \n (8) Kind \n (9) Completeness \n (10) Redz (Muzzins rainbow plot) \n (11) Rest-frame Z mag \n (12) Central flag (1=central, 0= neighbour of central) \n (13) Group # (Used to identify objects of the same group) \n (14) Areas Mpc^2 of the circular environments from which they were drawn.')
    
    return(finalta, len(cleanc))
###############################################################################

def determine_areas(fullcat, centrals, maglim, rad, mrad, final):

    ## Number of circles (or environments in the deep fields)   the number of circles*area_each_circle = surface areas of all neighbourhoods. Up until really bright mags there is little overlap 
    n = len(centrals[:,0])#number of circles analyzed in the deep fields
    
    # kind, mag, wd deep wide of all
    
    #Area of each circle in Mpc^2
    areac = np.pi*(mrad)**2
    
    #Number of environments analyzed*area of each = area of all circles
    #Note that this assumes no overlap btw circles which is true for brighter centrals but check for fainter mags
    #This overlap will become important if you are trying to compare with the SMF of background
    
    areat = n*areac
   
    return(areat)
       
####################################################################################################  
    
    
def completeness(cat):
    
    '''
    Find the weight associted to each galaxy, due to completeness correction. 

    '''
    
    wt=[]
    
    pd1 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/bias/deep_completeness/pecomp3-d1.dat')
    pd2 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/bias/deep_completeness/pecomp3-d2.dat')
    pd3 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/bias/deep_completeness/pecomp3-d3.dat')
    pd4 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/bias/deep_completeness/pecomp3-d4.dat')
   
    sd1 = np.genfromtxt('/Users/osejo/Desktop/SEDfit/2011/complsim/d1/sfcomp2-d1.dat')
    sd2 = np.genfromtxt('/Users/osejo/Desktop/SEDfit/2011/complsim/d2/sfcomp2-d2.dat')
    sd3 = np.genfromtxt('/Users/osejo/Desktop/SEDfit/2011/complsim/d3/sfcomp2-d3.dat')
    sd4 = np.genfromtxt('/Users/osejo/Desktop/SEDfit/2011/complsim/d4/sfcomp2-d4.dat')
   
    
    fp1 = interpolate.interp1d(pd1[:,0], pd1[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)    
    fp2 = interpolate.interp1d(pd2[:,0], pd2[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)
    fp3 = interpolate.interp1d(pd3[:,0], pd3[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)
    fp4 = interpolate.interp1d(pd4[:,0], pd4[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)
    
    fs1 = interpolate.interp1d(sd1[:,0], sd1[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)    
    fs2 = interpolate.interp1d(sd2[:,0], sd2[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)
    fs3 = interpolate.interp1d(sd3[:,0], sd3[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)
    fs4 = interpolate.interp1d(sd4[:,0], sd4[:,3], kind='cubic', axis=-1, copy=True, bounds_error=False, assume_sorted=False)
    
    for item in cat:
        
        wi=0
        field = item[0]
        kind = item[8]
        ks = item[6]
        
        if (field==11): #deep fields
            if(kind==0):
                wi = fp1(ks)
            elif(kind==1):
                wi = fs1(ks)
                
        elif (field==12): #deep fields
            if(kind==0):
                wi = fp2(ks)
            elif(kind==1):
                wi = fs2(ks)
                
        elif (field==13): #deep fields
            if(kind==0):
                wi = fp3(ks)
            elif(kind==1):
                wi = fs3(ks)
                
        elif (field==14): #deep fields
            if(kind==0):
                wi = fp4(ks)
            elif(kind==1):
                wi = fs4(ks)
                
        elif (field>14): # wide fields
            wi = 1
        
        if (wi<1):
            wi=1
            
        
        #print(wi, field, ks, kind)
        wt.append(wi)
        
    return (wt)
    
################################################################################ 

def find_areasmpc_wcat(cat):
    
    aw = 15.53+9.56 #areas of w1+w4 in deg^2
    am = aw*(60**2)*(60**2)*(8.471**2)/(1000**2) #areas of w1+w4 in Mpc^2
    print('Areas in deg of field', aw, 'Same area in Mpc**2', am) 
    at = np.ones(len(cat[:,0]))*am
    
    print(len(at))    
    return(at)
    
#################################################################################################### 
    
def find_areasmpc_dcat(cat):
    
    ad = 2.51
    ad2 = 0.911
    at = []
    
    for el in cat:
        
        ai = 0
        ks = el[6]
        kind = el[8]
        
        if(kind==0):
            
            if(ks<23):
                ai = ad
            elif(ks>=23):
                ai = ad-ad2
            
        elif(kind==1):
            
            ai = ad
        
        am = ai*(60**2)*(60**2)*(8.471**2)/(1000**2)
        at.append(am)
   
    print(len(at))    
    return(at)
    
    
#################################################################################################### 

def find_areasmpc_allcat(cat):
    
    aw = 15.53+9.56
    ad = 2.51
    ad2 = 0.911
    at = []
    print(len(cat[:,0]))
    
    for el in cat:
        ai=0
        kind = el[8]
        
        if (kind==0):
            if(el[6]<20.5):
                ai = aw+ad
            elif(el[6]>=20.5)and(el[6]<23):
                ai = ad
            elif(el[6]>=23):
                ai = ad - ad2
                
        elif(kind==1):
            if(el[6]<21.5):
                ai = aw+ad
            elif(el[6]>=21.5)and(el[6]<23.5):
                ai = ad
            
        am = ai*(60**2)*(60**2)*(8.471**2)/(1000**2)
        at.append(am)
    
    print(len(at))    
    return(at)

####################################################################################################


def get_redshifts(cat, kind):
    
   print('there is nothing here')
        
    
###############################################################################
    
def get_redshifs_field(wd):
    
    sfz = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/peakz_ks_sf.dat')
    pez = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/peakz_ks_pe.dat')
    
    sfzi = interpolate.interp1d(sfz[:,0], sfz[:,1], kind='linear', axis=-1, copy=True, bounds_error=False, assume_sorted=False, fill_value="extrapolate")
    pezi = interpolate.interp1d(pez[:,0], pez[:,1], kind='linear', axis=-1, copy=True, bounds_error=False, assume_sorted=False, fill_value="extrapolate")
    
    mags = np.arange(17,24,0.1)
    redzt = []
    
    # Look a the interpolation
    plt.figure()
    plt.xlim(17,24)
    plt.ylim(1.5,1.8)
    plt.plot(mags, pezi(mags), color='red', marker='o', linestyle='--')
    plt.plot(mags, sfzi(mags), color='blue', marker='o', linestyle='--')
    
    for gal in wd:
        
        
        redz = 0
        kind = gal[8]
        ks = gal[6]
        
        if(kind==0):
            
            redz = pezi(ks)
            
        elif(kind==1):
            
            redz = sfzi(ks)
           
            
        redzt.append(redz)
    
    #print(np.unique(redzt))
    return(redzt)
#################################################################################################### 
    
def get_abs_full(wd):
        
    ssp = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/sed_fit/make_bbsed/ssp_lr62_bc2003_z001.bbsed')
    agepe = 9.052783
    ebvpe = 0
    
    ssppe = ssp[ssp[:,2]==ebvpe]#Models that more closely reprsent  our PE population
    ssppe = ssppe[ssppe[:,4]==agepe]
    restfzp = ssppe[ssppe[:,0]==0][:,12]
   
    cons = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/sed_fit/make_bbsed/cons_lr62_bc2003_z001.bbsed')
    agesf = 8.324639
    ebvsf = 0.3
    
    conssf = cons[cons[:,2]==ebvsf]#Models that more closely reprsent  our PE population
    conssf = conssf[conssf[:,4]==agesf]
    restfzs = conssf[conssf[:,0]==0][:,12] #Rest-frame Z magnitude for this model
    
    dmagt = []
    abzt = []
    
    cat = 0
    zmagrest = -99
    
    for gal in wd:
        
        kind = gal[8]
        ks = gal[6]
        redz = gal[11]
        
        if(kind==0):
            cat = ssppe
            zmagrest = restfzp
        elif(kind==1):
            cat = conssf
            zmagrest = restfzs
    

        for item in cat:
            
            if(item[0]==round(redz, 2)):
                
                kmag = item[15] #K magnitude for this model at this peak redshift
                dmag = kmag - zmagrest #Delta Mag = Kmag model - Restframez
                abz = ks - dmag # Absolute z mag = ksmag - deltamag
                #print(redz, kind, item[0], ks, zmagrest, kmag, dmag, abz) ##To check output
                dmagt.append(dmag)
                abzt.append(float(abz))
    
    #array = np.c_[binredz, dmagt, abzt]
    #print(array)
    return(abzt)
    
    
    
###############################################################################


    
def make_cats(work_with, maglim, mrad):
    
    '''
    Only needs to be run once... it will create a full catalog of SF+PE galaxies 
    in the field so in the end you will have 3 catalogs, one for only wide, another
    for only deep and a third one that will combine both
    They will contain an incompleteness correction, areas, redshift and rest-frame Z mags
    '''
    
    wd = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/sfpe_dw_cd_nostr_maglim.dat') #This catalog is already down to the completeness magnitudes of each field
    wd = wd[wd[:,6]<20.5] #For the wide and all catalogs 20.5 for deep 23
    print(wd.shape)
    
    
    # Decide if you want to work with the full cat or only wide or deep, if both comment out next two lines
    if(work_with=='wide'):
        wd = wd[wd[:,0]>14] # work only with wide fields
    elif(work_with=='deep'):
        wd = wd[wd[:,0]<=14] # work only with the deep fields
    elif(work_with=='all'):
        wd = wd
    
    
    # Field LF Full catalog
    
    # Find the completeness All
    wt = completeness(wd) #Completeness for the entire catalog, based on field and magnitude
    wd = np.c_[wd, wt] #Column 9 will represent now the completeness of the sample
    
    
    # Areas full cat
    
    if(work_with=='wide'):
        at = find_areasmpc_wcat(wd)
    elif(work_with=='deep'):
        at = find_areasmpc_dcat(wd)
    elif(work_with=='all'):
        at = find_areasmpc_allcat(wd)
    
    wd = np.c_[wd, at] #Column 10 will be the surface area of the sample in physical Mpc^2
    print(wd.shape, np.unique(at))
    
    
    zfield = get_redshifs_field(wd)
    wd = np.c_[wd, zfield] #Col 11 will be an approx number for redz
    abszm = get_abs_full(wd)
    wd = np.c_[wd, abszm] # Col 12 will be absolute Z magnitudes
    #print(min(wd[:,12]), max(wd[:,12]))
    
    #np.savetxt('test_fullcat.dat', wd, fmt='%s')
    
    
    
###############################################################################

def select_nb(fullcat, maglim, mrad, work_with):
    
    
    '''
    In this script, first you will select all centrals (only passive galaxies) from the fullcatalog
    Once you have ID centrals based on quiescence and a limiting mag, the create_cat_nb
    will create a catalog of the neighbours of these centrals
    '''
    
    #print ('Unique fields ', np.unique(fullcat[:,0]))
    print ('How many SF+PEGs we have in the wide fields: ', len(fullcat[:,0]),'To a max Ks of',  max(fullcat[:,6]))
    
    ### Select Centrals
    pec = fullcat[fullcat[:,8]==0] #only passive galaxies will be considered as centrals
    central = pec[(pec[:,6]<maglim) & (pec[:,6]>artifact)] #massive centrals given a limiting magnitude established in global variable maglim and also making extra sure there are no  artifacts (even though they were removed from full cat anyway)
    
    print('Central min mag', min(central[:,6]), 'Central max mag ', max(central[:,6]))
    print('This is how many centrals I am working with ', len(central))
    
    
    ## Find neighbours around all centrals    
    neighbours, ncentrals = create_cat_nb(maglim, mrad, fullcat, central) #Returns a catalog of neighbours and the number of trully centrals
      
    
    print(neighbours.shape)
    return(neighbours, ncentrals)
    
###############################################################################

def make_fraction_graph(neib):
    
    '''
    Make a graph with the fractional absolute magnitudes
    of neighbouring galaxies with respect to their centrals
    '''
    i = 0 # Go through all groups
    ngroups = np.unique(neib[:,13])
    print('how many groups', len(ngroups))
    
    array = []
    
    while(i<=len(ngroups)): #Go through all the groups
                
        groupi = neib[neib[:,13]==i] #ID all the elements in group ith
        cent = groupi[groupi[:,12]==1] #Identify the central in group ith
        zcent = cent[:,11] # ID the rest-frame Z mag of the central in group ith
        
        #print('Group i', i, groupi) #for testing purposes
        #print('Central', i, cent, kscent)
        
        for item in groupi: #for every item in this group
            vec = []
            f = -99 #f represents fraction not flux. So we want a fractional representation 
            
            # In my opinion the following if is redundant
            
            if(item[12]==1): #If you were ID as the central
                f = 1
            elif(item[12]==0):
                #f = 10**(-0.4*item[11])/10**(-0.4*zcent) # I wrote the same thing in a different way to easily remember where it comes from
                f = 10**(-0.4*(item[11]-zcent))
                
            vec = np.hstack((item, f))
            
            if(len(array)==0):
                array = vec
            else:
                array = np.vstack((array, vec))
                
        i = i+1
    
        
    #np.savetxt('test3_'+work_with+'_nb_sfpeg_ks'+str(maglim)+'_r'+str(mrad)+'.dat', array, fmt='%s', header='Combined catalog of SF+PE galaxies in a neighbouring '+str(mrad)+' Mpc radius of PE galaxies brighter than '+str(maglim)+' in '+work_with+' fields, calculated with env_dense_pegs3.py \n (0) Field \n (1) ID  \n (2) X [pix] \n (3) Y [pix] \n (4) RA [deg] \n (5) DEC [deg] \n (6) Ktot [AB mag] \n (7) Ktot err [AB mag] \n (8) Kind \n (9) Completeness \n (10) Redz (Muzzins rainbow plot) \n (11) Absolute Z mag \n (12) Central flag (1=central, 0= neighbour of central) \n (13) Group # (Used to identify objects of the same group) \n (14) Areas Mpc^2 of the circular environments from which they were drawn \n (15) Fraction (Z/Zcentral)')
    
    return(array)

###############################################################################

def background_correction(fullcat, arrayf, mrad, iterations): 
    
    '''
    NO LONGER USED
    Take random centrals (regardless of magnitude) and find neighbouring
    galaxies around this random central in the established radius
    Found galaxies will help to determine foreground/background contamination
    
    Soooo the problem with this approach is that passive galaxies are possibly already
    somewhat clustered, so the sample is not a baseline. Also the central was not necessarily
    more massive than the rest so if I tried to correct the fractional graph it won't work.
    You need to correct for background/foreground contamination, it needs to be done BEFORE
    you find the fractional 
    '''
    i = 1
    n = iterations #number of iterations
    j = 1 #ID of each group, same number of groups as centrals
    ncentrals = len(arrayf[arrayf[:,12]==1]) # how many centrals we got from the REAL SAMPLE 
    rad = mrad*1000/8.471/60/60 #from mpc to deg radius of search
    
    #print('There are ', ncentrals, ' centrals in the sample') #for testing purposes
    treer = cKDTree(zip(fullcat[:,4],fullcat[:,5])) #Make KDTree to search for nbs
    arrayf = []    
    
    while(i<=n):
        #Select random centrals, the same number of centrals as the real ones
        nrand = np.random.randint(0, ncentrals+1, ncentrals) #ncentrals
        rcentrals = fullcat[nrand][:]#random galaxies that will serve as centrals
        #print('Centrals', rcentrals[:,1])
        
        for each in rcentrals:
            
            zcr = each[11] # Z abs magnitude of central
            distr, indr = treer.query((each[4], each[5]), k=50, distance_upper_bound=rad)
            result = np.c_[distr, indr]
            #print('Nbs around this central ')
            #print(result)
            
            ## Make a full catalog of all neighbouring galaxies found
            
            for el in result:
                vec = []
                area = 51.05088062
                #print(el[1])
                centralf = -99 #Central flag 1= yes, 0 =no
                f = -99 #fractional Zrest-frame mag compared with the central
                
                if(el[0]==0):
                    centralf = 1
                    f = 1
                    #print('Himself', el[0], fullcat[el[1]][1])
                    vec = np.hstack((fullcat[el[1]][0], fullcat[el[1]][1], fullcat[el[1]][2], fullcat[el[1]][3], fullcat[el[1]][4], fullcat[el[1]][5], fullcat[el[1]][6], fullcat[el[1]][7], fullcat[el[1]][8], fullcat[el[1]][9], fullcat[el[1]][11], fullcat[el[1]][12], centralf, j, area, f))
                    
                    #if(len(arrayf)==0):
                        #arrayf = vec
                    #else:
                        #arrayf = np.vstack((arrayf, vec))
                        
                elif(el[0]>0)and(el[0]!=np.inf): #Do not select central AND galaxies outside of rad (tree.query flags as inf galaxies outside the radius of search)
                    centralf = 0
                    f = (10**(-0.4*fullcat[el[1]][11]))/(10**(-0.4*zcr))
                    #print('Companions', el[0], fullcat[el[1]][1])
                    vec = np.hstack((fullcat[el[1]][0], fullcat[el[1]][1], fullcat[el[1]][2], fullcat[el[1]][3], fullcat[el[1]][4], fullcat[el[1]][5], fullcat[el[1]][6], fullcat[el[1]][7], fullcat[el[1]][8], fullcat[el[1]][9], fullcat[el[1]][11], fullcat[el[1]][12], centralf, j, area, f))
                
                    if(len(arrayf)==0):
                        arrayf = vec
                    else:
                        arrayf = np.vstack((arrayf, vec))
                    

            j = j+1
            
            
        i = i+1
    print('Shape random NBs', arrayf.shape, 'Number of iteratios', n)   
    #np.savetxt('test_random_'+work_with+'_nb_sfpeg_ks_r'+str(mrad)+'.dat', arrayf, fmt='%s', header='Combined catalog of SF+PE galaxies in a neighbouring '+str(mrad)+' Mpc radius of random galaxies in '+work_with+' fields, calculated with env_dense_pegs3.py \n (0) Field \n (1) ID  \n (2) X [pix] \n (3) Y [pix] \n (4) RA [deg] \n (5) DEC [deg] \n (6) Ktot [AB mag] \n (7) Ktot err [AB mag] \n (8) Kind \n (9) Completeness \n (10) Redz (Muzzins rainbow plot) \n (11) Absolute Z mag \n (12) Central flag (1=central, 0= neighbour of central) \n (13) Group # (Used to identify objects of the same group) \n (14) Areas Mpc^2 of the circular environments from which they were drawn. \n (15) Z/Zcentral')
    return(arrayf)
    
###############################################################################
def test_fraction(neib, bcg, mrad, binsf): #This test was done with background correction fo all, not divided
    numgroup = 22 #test group it has 3 elements, including its central

    group = neib[neib[:,13]==numgroup]
    nbgroup = group[group[:,12]<1] #Neighbours are ONLY those that are not centrals
    central = group[group[:,12]==1]
    magcentral = central[11]
    areag = (np.pi)*(mrad**2) #area of each group in approx numbers, this is a test, but for the real one it will have the ACTUAL Good pix inside of radius of search
    
    print('Len bins', len(binsf))
    print('Shape of group '+str(numgroup)+' is '+str(group.shape))
    
    print('Group', group[:,6], group[:,12], group[:,11])
    print('NBs', nbgroup[:,6], nbgroup[:,12], nbgroup[:,11])
    print('Central', central[:,6], central[:,12], 'mag ', magcentral)

    ## As a hist    
    centerg = (binsf[:-1] + binsf[1:]) / 2 #Use for plotting

    hg, hgedges = np.histogram(nbgroup[:,11], binsf, weights=nbgroup[:,9]/areag) #Delta Histogram of counts in group
    hge, hgedgese = np.histogram(nbgroup[:,11], binsf, weights=nbgroup[:,9]) #Delta Histogram of counts in group not per area, since these will be used for the SQRT(N) errors
    
    harray = np.c_[hg, centerg]
    harray = harray[harray[:,0]>0]
    bghg = harray[:,0] - bcg(harray[:,1]) 
    print(harray[:,1], bghg)
    
    #np.savetxt('test.dat', np.c_[hg, centerg]) #if you want to check the delta hist of neighbours
    
    ## As a #
    ones  = np.ones(len(nbgroup[:,6]))/areag #No need for completeness correction when working with the wide field, as is the case in this test. But apply if working with Deep
    nbarray = np.c_[nbgroup, ones]    
    #bgn = nbarray[:,15]-bcg(round(nbarray[:,11], abs(np.log10(binsize))))
    bgn = nbarray[:,15]-bcg(nbarray[:,11])
    print(nbarray[:,11], bgn)
    
    
#########################################################################################

    
#########################################################################################
    
    
def apply_corr_redz(f, fs, fp):
    print('NOthing here')
    
    '''
    ncf = 1 #don't have rainbow plot for these yet
    ncfs = 1
    
    sigma = 2.0 # width of 2d gaussian kernel for smoothing #2.0 for PEGs, 3.0 for SFGs
    zBinWidth = 0.1
    kBinWidth = 0.5


    # define bin edges
    zbins=np.arange(0.0,3.5,zBinWidth)#1-3 original values 0.5-3.5 for PEGs 0-3.5 for SFGs
    #kbins=np.arange(min(d2k),max(d2k)+0.5,kBinWidth) #gives it a bit more freedom
    kbins=np.arange(17, 24, kBinWidth) #24 0or 21

    d2=hist2(d2z,d2k,zbins,kbins)
    zmids=d2[1]
    kmids=d2[2]
    
    pzmpe = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/fraction/pzm_correction.dat')

    pcentrals = pzmpe[:,5] #Column that represents 19.5, the most likely magnitude of our centrals
    pfcent = interp1d(zmids,centrals)
    pfcent2 = interp1d(zmids, centrals**2)
    
    
    

    plt.figure()
    plt.title('Redshift distribution galaxies at Ks= 19.5')
    plt.plot(zmids, pfcent(zmids), color='red', linestyle='--')

    print(zmids[0]*1.0001,zmids[-1]*0.9999)
    a = quad(fcent,zmids[0]*1.0001,zmids[-1]*0.9999)
    a2 = quad(fcent2,zmids[0]*1.0001,zmids[-1]*0.9999)
    print(a[0]/a2[0])
 
 
 
    return(ncf, ncfs, ncfp)
    '''

###############################################################################

def fractional(neib, bcgo, bcgs, bcgp, mrad, binsf, binsp, binss, ncentrals):
    
    '''
    25/01/17
    Make a catalog that will allow us to compare the relative sizes
    of the companions compared to centrals
    '''
    
    #print('Len bins', len(binsf))
    areage = (np.pi)*(mrad**2) #area of each group in approx numbers, this is a test, but for the real one it will have the ACTUAL Good pix inside of radius of search
    areag = areage*ncentrals 
    
    # Rough estimate less areas--for testing only
    areag1 = areage*(ncentrals-5)
    areag2 = areage/3*(5)
    #areag = (areag1+areag2)
        
    
    print('Area of each group', areage, 'Mpc^2') 
    print('Area all groups', areag, 'Mpc^2') 
    
    groups = np.unique(neib[:,13])
    i = 0
    tarray = []
    harraysft = []
    harraysft2 = [] 
    
    harraypft = []
    harraypft2 = []
    
    harraysf = []
    harraypf = []
    
    while(i<len(groups)):
        #print('Group', i)
        
        ##Initialize each array for each group
        vec = 0
        vecn = 0
        harrays = []
        harrayp = []
        harraysc = []
        harraypc = []
        group = []
        nbgroup = []
        central = []
        sfnb = []
        penb = []
        
        group = neib[neib[:,13]==i]
        nbgroup = group[group[:,12]<1] #Neighbours are ONLY those that are not centrals
        central = group[group[:,12]==1] #ID the central 
        vec = np.hstack((central[:,11], 1/areag, 0, 1, 1))
        
        magcentral = central[:,11]
        #print('Mag central', magcentral)
        
        if(len(tarray)==0):
            tarray = vec
        else:
            tarray = np.vstack((tarray,vec))
        
        sfnb = nbgroup[nbgroup[:,8]==1] #SF Neighbours 
        penb = nbgroup[nbgroup[:,8]==0] #PE Neighbours
        
        
        ##For testing
        
        #print('Shape of group '+str(i)+' is '+str(group.shape))
    
        #print('Group', group[:,6], group[:,12], group[:,11])
        #print('NBs', nbgroup[:,6], nbgroup[:,12], nbgroup[:,11])
        #print('Central', central[:,6], central[:,12])
        #print('SF NBS', sfnb[:,6], sfnb[:,8], sfnb[:,11])
        #print('PE NBS', penb[:,6], penb[:,8], penb[:,11])
        
        
        ## Make a narrow histogram of SF and PE neighbours
        centernbs = (binss[:-1] + binss[1:]) / 2 #Use for plotting
        centernbp = (binsp[:-1] + binsp[1:]) / 2 #Use for plotting
        
        import correct_frac_setup as correction
    
        ncfp2, ncfsa2, ncfa2 = correction.get_values()
        #ncfp2 = 1

        hgs, hgedgess = np.histogram(sfnb[:,11], binss, weights=ncfp2*sfnb[:,9]/areag) #Delta Histogram of counts in group
        hges, hgedgeses = np.histogram(sfnb[:,11], binss, weights=ncfp2*sfnb[:,9]) #Delta Histogram of counts in group not per area, since these will be used for the SQRT(N) errors
    
        hgp, hgedgesp = np.histogram(penb[:,11], binsp, weights=ncfp2*penb[:,9]/areag) #Delta Histogram of counts in group
        hgep, hgedgesep = np.histogram(penb[:,11], binsp, weights=ncfp2*penb[:,9]) #Delta Histogram of counts in group not per area, since these will be used for the SQRT(N) errors
    
        
        
        ## Remove bins that are empty 
        harrays = np.c_[centernbs, hgs]
        harrays = harrays[harrays[:,1]>0]
        
        harrayp = np.c_[centernbp, hgp]
        harrayp = harrayp[harrayp[:,1]>0]
        
        ###For testing
        
        #print('HIST SF', harrays[:,0], harrays[:,1])
        #print('HIST PE', harrayp[:,0], harrayp[:,1])
        
        
        
        ## Do a background correction
        bgsf = func4(harrays[:,0], bcgs[0][0], bcgs[0][1], bcgs[0][2], bcgs[0][3], bcgs[0][4])
        harraysc = harrays[:,1] - bgsf 
    
        
        bgpe = func4(harrayp[:,0], bcgp[0][0], bcgp[0][1], bcgp[0][2], bcgp[0][3], bcgp[0][4])
        harraypc = harrayp[:,1] - bgpe 
        
        
        harraysf = np.c_[harrays[:,0], harraysc] 
        harraypf = np.c_[harrayp[:,0], harraypc]
        
        ### For testing
        
        #print('SF Corr hist', harrays[:,0], harrays[:,1], bgsf, harraysc)
        #print('PE Corr hist', harrayp[:,0], harrayp[:,1], bgpe, harraypc)
        
        
        
        fs = 10**(-0.4*(harraysf[:,0]-magcentral)) #for SF
        fp = 10**(-0.4*(harraypf[:,0]-magcentral)) #for PE
        
        ## Add to final array the fractional comparison of each NB
        harraysft = np.c_[harraysf, np.ones(len(harraysf)), fs, np.zeros(len(harraysf))] 
        harraypft = np.c_[harraypf, np.zeros(len(harraypf)), fp, np.zeros(len(harraypf))]
        
        tarray = np.vstack((tarray, harraysft, harraypft))
        
        
        i=i+1
        
    
        #tarray cols
        #(0) Rest-frame Z mag
        #(1) #/area/bin # In Mpc^2
        #(2) SF or PE flag, SF = 1, PE = 0
        #(3) Fractional flux of each nb relative to their central (In rest-frame Z)
        #(4) Central = 1 Neighbour = 0
    
    
    #np.savetxt('test_tarray.dat', tarray, fmt='%s')
    
    ### Scatter plot, for testing purposes
    #plt.figure()
    #plt.yscale('log')
    #plt.xlim(0.02,1.4)
    #plt.ylim(1,1.5)
    #plt.scatter(tarray[:,3], tarray[:,1])
    
    #### Bin SF+PE
    
    bin_centers = np.arange(0.0,1.1,0.1)
    bin1=0
    bin2=0
    bin3=0
    bin4=0
    bin5=0
    bin6=0
    bin7=0
    bin8=0
    bin9=0
    bin10=0
    bin11=0
    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    i6 = 0
    i7 = 0
    i8 = 0
    i9 = 0
    i10 = 0
    i11 = 0
     
    
    for el in tarray:
        
        # To start, a simple multiple ifs
        if(el[3]>0.95): #basically centrals
            bin1 = bin1+el[1]
            i1 = i1 + 1
        
        elif(el[3]>0.85)and(el[3]<=0.95):
            bin2 = bin2+el[1]
            i2 = i2 + 1
            
        
        elif(el[3]>0.75)and(el[3]<=0.85):
            bin3 = bin3+el[1]
            i3 = i3 + 1
            
        
        elif(el[3]>0.65)and(el[3]<=0.75):
            bin4 = bin4+el[1]
            i4 = i4 + 1
            
        
        elif(el[3]>0.55)and(el[3]<=0.65):
            bin5 = bin5+el[1]
            i5 = i5 + 1
            
        
        elif(el[3]>0.45)and(el[3]<=0.55):
            bin6 = bin6+el[1]
            i6 = i6 + 1
            
        
        elif(el[3]>0.35)and(el[3]<=0.45):
            bin7 = bin7+el[1]
            i7 = i7 + 1
            
        
        elif(el[3]>0.25)and(el[3]<=0.35):
            bin8 = bin8+el[1]
            i8 = i8 + 1
            
        
        elif(el[3]>0.15)and(el[3]<=0.25):
            bin9 = bin9+el[1]
            i9 = i9 + 1
            
            
        elif(el[3]>0.05)and(el[3]<=0.15):
           bin10 = bin10+el[1]
           i10 = i10 + 1
           
        elif(el[3]>-0.05)and(el[3]<=0.05):
           bin11 = bin11+el[1]
           i11 = i11 + 1
            
    
    #### Bin SF
    
    
    sbin1=0
    sbin2=0
    sbin3=0
    sbin4=0
    sbin5=0
    sbin6=0
    sbin7=0
    sbin8=0
    sbin9=0
    sbin10=0
    sbin11=0
    
    si1 = 0
    si2 = 0
    si3 = 0
    si4 = 0
    si5 = 0
    si6 = 0
    si7 = 0
    si8 = 0
    si9 = 0
    si10 = 0
    si11 = 0
    
   
    
     
    harraysft2 = tarray[tarray[:,2]==1]
    #print('SF', harraysft2)
    
    for els in harraysft2:
        
        # To start, a simple multiple ifs
        if(els[3]>0.95): #basically centrals
            sbin1 = sbin1+els[1]
            si1 = si1 + 1       
        
        elif(els[3]>0.85)and(els[3]<=0.95):
            sbin2 = sbin2+els[1]
            si2 = si2 + 1
        
        elif(els[3]>0.75)and(els[3]<=0.85):
            sbin3 = sbin3+els[1]
            si3 = si3 + 1
        
        elif(els[3]>0.65)and(els[3]<=0.75):
            sbin4 = sbin4+els[1]
            si4 = si4 + 1
        
        elif(els[3]>0.55)and(els[3]<=0.65):
            sbin5 = sbin5+els[1]
            si5 = si5 + 1
        
        elif(els[3]>0.45)and(els[3]<=0.55):
            sbin6 = sbin6+els[1]
            si6 = si6 + 1
        
        elif(els[3]>0.35)and(els[3]<=0.45):
            sbin7 = sbin7+els[1]
            si7 = si7 + 1
        
        elif(els[3]>0.25)and(els[3]<=0.35):
            sbin8 = sbin8+els[1]
            si8 = si8 + 1
        
        elif(els[3]>0.15)and(els[3]<=0.25):
            sbin9 = sbin9+els[1]
            si9 = si9 + 1
            
        elif(els[3]>0.05)and(els[3]<=0.15):
            sbin10 = sbin10+els[1]
            si10 = si10 + 1
           
        elif(els[3]>-0.05)and(els[3]<=0.05):
            sbin11 = sbin11+els[1]
            si11 = si11 + 1
    
    
    pbin1=0
    pbin2=0
    pbin3=0
    pbin4=0
    pbin5=0
    pbin6=0
    pbin7=0
    pbin8=0
    pbin9=0
    pbin10=0
    pbin11=0
    
    pi1 = 0
    pi2 = 0
    pi3 = 0
    pi4 = 0
    pi5 = 0
    pi6 = 0
    pi7 = 0
    pi8 = 0
    pi9 = 0
    pi10 = 0
    pi11 = 0
    
    
    harraypft2 = tarray[tarray[:,2]==0]
    #print('PE', harraypft2)
    
    for elp in harraypft2:
        
        # To start, a simple multiple ifs
        if(elp[3]>0.95): #basically centrals
            pbin1 = pbin1+elp[1]
            pi1 = pi1 + 1
        
        elif(elp[3]>0.85)and(elp[3]<=0.95):
            pbin2 = pbin2+elp[1]
            pi2 = pi2 + 1
        
        elif(elp[3]>0.75)and(elp[3]<=0.85):
            pbin3 = pbin3+elp[1]
            pi3 = pi3 + 1
        
        elif(elp[3]>0.65)and(elp[3]<=0.75):
            pbin4 = pbin4+elp[1]
            pi4 = pi4 + 1
        
        elif(elp[3]>0.55)and(elp[3]<=0.65):
            pbin5 = pbin5+elp[1]
            pi5 = pi5 + 1
        
        elif(elp[3]>0.45)and(elp[3]<=0.55):
            pbin6 = pbin6+elp[1]
            pi6 = pi6 + 1
        
        elif(elp[3]>0.35)and(elp[3]<=0.45):
            pbin7 = pbin7+elp[1]
            pi7 = pi7 + 1
        
        elif(elp[3]>0.25)and(elp[3]<=0.35):
            pbin8 = pbin8+elp[1]
            pi8 = pi8 + 1
        
        elif(elp[3]>0.15)and(elp[3]<=0.25):
            pbin9 = pbin9+elp[1]
            pi9 = pi9 + 1
            
        elif(elp[3]>0.05)and(elp[3]<=0.15):
           pbin10 = pbin10+elp[1]
           pi10 = pi10 + 1
           
        elif(elp[3]>-0.05)and(elp[3]<=0.05):
           pbin11 = pbin11+elp[1]
           pi11 = pi11 + 1
    

    
    #import correct_frac_setup as correction #Doesn't make sense to do it here
    
    #ncfp, ncfsa, ncfa = correction.get_values()
    
    ncfp = 1
    

    f = np.hstack((bin11*ncfp, bin10*ncfp, bin9*ncfp, bin8*ncfp, bin7*ncfp, bin6*ncfp, bin5*ncfp, bin4*ncfp, bin3*ncfp, bin2*ncfp, bin1))
    fs = np.hstack((sbin11*ncfp, sbin10*ncfp, sbin9*ncfp, sbin8*ncfp, sbin7*ncfp, sbin6*ncfp, sbin5*ncfp, sbin4*ncfp, sbin3*ncfp, sbin2*ncfp, sbin1))
    fp = np.hstack((pbin11*ncfp, pbin10*ncfp, pbin9*ncfp, pbin8*ncfp, pbin7*ncfp, pbin6*ncfp, pbin5*ncfp, pbin4*ncfp, pbin3*ncfp, pbin2*ncfp, pbin1))
    
    
    fe = np.hstack((np.sqrt(i11*ncfp)/areag, np.sqrt(i10*ncfp)/areag, np.sqrt(i9*ncfp)/areag, np.sqrt(i8*ncfp)/areag, np.sqrt(i7*ncfp)/areag, np.sqrt(i6*ncfp)/areag, np.sqrt(i5*ncfp)/areag, np.sqrt(i4*ncfp)/areag, np.sqrt(i3*ncfp)/areag, np.sqrt(i2*ncfp)/areag, np.sqrt(i1)/areag))
    fes = np.hstack((np.sqrt(si11*ncfp)/areag, np.sqrt(si10*ncfp)/areag, np.sqrt(si9*ncfp)/areag, np.sqrt(si8*ncfp)/areag, np.sqrt(si7*ncfp)/areag, np.sqrt(si6*ncfp)/areag, np.sqrt(si5*ncfp)/areag, np.sqrt(si4*ncfp)/areag, np.sqrt(si3*ncfp)/areag, np.sqrt(si2*ncfp)/areag, np.sqrt(si1)/areag))
    fep = np.hstack((np.sqrt(pi11*ncfp)/areag, np.sqrt(pi10*ncfp)/areag, np.sqrt(pi9*ncfp)/areag, np.sqrt(pi8*ncfp)/areag, np.sqrt(pi7*ncfp)/areag, np.sqrt(pi6*ncfp)/areag, np.sqrt(pi5*ncfp)/areag, np.sqrt(pi4*ncfp)/areag, np.sqrt(pi3*ncfp)/areag, np.sqrt(pi2*ncfp)/areag, np.sqrt(pi1)/areag))
    
    
    
    plt.figure(figsize=(10,7))
    plt.yscale('log')
    plt.xlim(-0.1,1.05)
    
    if(work_with=='wide'):
        plt.ylim(10**-2,1.5*10**0)
    elif(work_with=='deep'):
        #plt.ylim(10**-1,2*10**1)
        
        limits = np.c_[f, fs, fp]
        #print(np.c_[bin_centers, limits])
        limits[limits == 0] = max(limits[:,1]) ### JUST FOR THE LIMITS!! TO make sure the min value is not -inf!
        plt.ylim((min(limits[:,1]))*0.5, (max(limits[:,1]))*2)
    
    
        
    plt.yscale('log')
    plt.title('Neighbouring galaxies around PEGs 18.5<Ks<'+str(maglim), fontsize=21)
    plt.xlabel('(L/L$_C$)', fontsize=19)
    plt.ylabel(r'$\Sigma$ (#/$Mpc^2$/(0.1 (L/L$_C$) bin)', fontsize=19)
    
    plt.xticks(np.arange(0,1.1,0.1))
    plt.scatter(bin_centers, f, marker='D', s=73, color='green', alpha=0.5) #with errorbars sqrt(N) after you remove area dependence
    plt.scatter(bin_centers, fs, marker='s', s=73, color='blue', alpha=0.5)
    plt.scatter(bin_centers, fp, marker='p', s=83, color='red', alpha=0.5)
    
    plt.errorbar(bin_centers, f, yerr =  fe, color='green', linestyle ='None', alpha=0.5, marker='D', ms=10)
    plt.errorbar(bin_centers, fs, yerr =  fes, color='blue', linestyle ='None', alpha=0.5, marker='s', ms=10)
    plt.errorbar(bin_centers, fp, yerr =  fep, color='red', linestyle ='None', alpha=0.5, marker='p', ms=12)
    
    #print(bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9)
    #print(sbin1, sbin2, sbin3, sbin4, sbin5, sbin6, sbin7, sbin8, sbin9)
    
    
#################################################################################
    
###############################################################################

def fractional2(neib, bcgo, bcgs, bcgp, mrad, binsf, binsp, binss, ncentrals):
    
    '''
    02/04/17
    Make a catalog that will allow us to compare the relative sizes
    of the companions compared to centrals.
    Unlike fractional this subroutine has a different background correction
    '''
    
    #print('Len bins', len(binsf))
    areage = (np.pi)*(mrad**2) #area of each group in approx numbers, this is a test, but for the real one it will have the ACTUAL Good pix inside of radius of search
    areag = areage*ncentrals 
    
    # Rough estimate less areas--for testing only
    areag1 = areage*(ncentrals-5)
    areag2 = areage/3*(5)
    #areag = (areag1+areag2)
        
    
    print('Area of each group', areage, 'Mpc^2') 
    print('Area all groups', areag, 'Mpc^2') 
    
    groups = np.unique(neib[:,13])
    i = 0
    tarray = []
    harraysft = []
    harraysft2 = [] 
    
    harraypft = []
    harraypft2 = []
    
    harraysf = []
    harraypf = []
    
    while(i<len(groups)):
        #print('Group', i)
        
        ##Initialize each array for each group
        vec = 0
        vecn = 0
        harrays = []
        harrayp = []
        harraysc = []
        harraypc = []
        group = []
        nbgroup = []
        central = []
        sfnb = []
        penb = []
        
        group = neib[neib[:,13]==i]
        nbgroup = group[group[:,12]<1] #Neighbours are ONLY those that are not centrals
        central = group[group[:,12]==1] #ID the central 
        vec = np.hstack((central[:,11], 1/areag, 0, 1, 1))
    
###############################################################################
  
def kind_of_neighborhood(neib, mrad):
    
    
    '''
    Created to solve questions like how many centrals have 0,1,2 or 3 companions
    And are these companions SF or PE?
    '''

    groups = np.unique(neib[:,13])
    print('There are ', len(groups), 'groups')
    i = 0
    alone = 0
    one = 0
    dos = 0
    tres = 0
    
    onesf = 0
    onepe = 0
    
    dossf = 0
    dospe = 0
    
    tressf = 0
    trespe = 0
    
    ### For each group of galaxies: A group represents the central+all found neighbouring galaxies in the defined radius: mrad
    while(i<len(groups)):
        #print('Group', i)
        kind = -99
        sfgroup = []
        pegroup = []
        
        sfgroup3 = []
        pegroup3 = []
        
        group = neib[neib[:,13]==i]
        nbgroup = group[group[:,12]<1] #Neighbours are ONLY those that are not centrals
        central = group[group[:,12]==1] #ID the central
        
        ## For testing purposes
        
        #print('Number of NBs', len(nbgroup))
        #print('Group', group[:,6], group[:,8], group[:,12])
        #print('Central', central[:,6], central[:,8], central[:,12])
        #print('NBs', nbgroup[:,6], nbgroup[:,8], nbgroup[:,12])
        
        
        if(len(nbgroup)==0):
            #print('I am alone')
            alone = alone + 1
            
            
        elif(len(nbgroup)==1):
            one = one + 1
            kind = np.unique(nbgroup[:,8])
            
            if(kind==1):
                onesf = onesf + 1
            elif(kind==0):
                onepe = onepe + 1
            else:
                print('Problem!')
            
        elif(len(nbgroup)==2):
            dos = dos + 1
            #print(nbgroup[:,6], nbgroup[:,8], nbgroup[:,12])
            sfgroup = nbgroup[nbgroup[:,8]==1]
            pegroup = nbgroup[nbgroup[:,8]==0]
            
            dossf = dossf + len(sfgroup)
            dospe = dospe + len(pegroup)
                
            
        elif(len(nbgroup)==3):
            tres = tres + 1
            
            sfgroup3 = nbgroup[nbgroup[:,8]==1]
            pegroup3 = nbgroup[nbgroup[:,8]==0]
            
            tressf = tressf + len(sfgroup3)
            trespe = trespe + len(pegroup3)
        
        i = i + 1
        
        
    #print('Total', alone+one+dos+tres, 'Alone', alone, 'One', one, 'Two', dos, 'Three', tres)
    #print('One', one, 'SF', onesf, 'PE', onepe)
    #print('Two', dos, 'SF', dossf, 'PE', dospe)
    #print('Three', tres, 'SF', tressf, 'PE', trespe)
    
    
    vecnb = np.arange(0,4,1)
    vecnb2 = np.r_[alone, one, dos, tres]
    vecnbs = np.r_[-99, onesf, dossf, tressf]
    vecnbp = np.r_[-99, onepe, dospe, trespe]
    
    
    nbarray = np.c_[vecnb, vecnb2]
    snbarray =  np.c_[vecnb, vecnbs]
    pnbarray =  np.c_[vecnb, vecnbp]
    
    ### Express SF and PE in %%
    total = nbarray[:,0]*nbarray[:,1] #these numbers will set the 100% of the stacked bar
    
    psf = snbarray[:,1]/total
    npsf = psf*nbarray[:,1]
    print(psf[1])
    ppe = pnbarray[:,1]/total
    nppe = ppe*nbarray[:,1]
    
    #print('Initial total array', nbarray)
    #print('Total # of gals', total)
    #print('Percentage SF', psf)
    #print('What that percentage will be in #', npsf)
    
    
    
    ### Numbers figure
    #plt.figure()
    #plt.xlim(-0.5,3.5)
    #plt.ylim(-5,40) #number of centrals
    #plt.xticks(np.arange(0,4,1))
    #plt.ylabel('# of centrals')
    #plt.xlabel('# of neighbours')
    #plt.title('Kinds of neighbours')
    #plt.scatter(nbarray[:,0], nbarray[:,1], marker=u'o', color='green', s=70)
   
    ### Bar graph
    plt.figure(figsize=(10,7))
    plt.xlim(-0.5,3.5)
    plt.xticks(np.arange(0,4,1))
    plt.ylabel('# of centrals', size=22)
    plt.xlabel('# of neighbours', size=22)
    plt.title('Given number of neighbours for centrals, Wide Fields', size=23)
    plt.ylim(0,40) #number of centrals
    plt.bar(nbarray[:,0]-0.2, nbarray[:,1], width=0.4, color='gray', alpha=0.5, bottom=0)
    
    plt.bar(snbarray[:,0]-0.2, npsf, width=0.4, color='blue', alpha=0.5, bottom=0)
    plt.bar(pnbarray[:,0]-0.2, nppe, width=0.4, color='red', alpha=0.5, bottom=npsf)
    
    
    plt.text(0.3,0.5*(nbarray[1][1]),'SF '+str(int(psf[1]*100))+' %', color='blue')
    plt.text(0.3,0.6*(nbarray[1][1]),'PE '+str(int(ppe[1]*100))+' %', color='red')
    
    plt.text(1.3,0.5*(nbarray[2][1]),'SF '+str(int(psf[2]*100))+' %', color='blue')
    plt.text(1.3,1.2*(nbarray[2][1]),'PE '+str(int(ppe[2]*100))+' %', color='red')
    
    
    plt.text(2.3,0.5*(nbarray[3][1]),'SF '+str(int(psf[3]*100))+' %', color='blue')
    plt.text(2.3,2*(nbarray[3][1]),'PE '+str(int(ppe[3]*100))+' %', color='red')
    
    '''
    ### 3D graph
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylim(-5,40)
    ax.set_xlim(-0.5,3.5)
    ax.set_zlim(-0.5,3.5)
    ax.set_xlabel('# of nbs')
    ax.set_ylabel('# of centrals')
    ax.set_zlabel('# of SF')
    ax.scatter(vecnb, vecnb2, vecnbs, c='blue', marker='o')
    '''
  
###############################################################################

def kind_of_neighborhood_deep(neib, mrad):
    
    
    '''
    Created to solve questions like how many centrals have 0,1,2 or 3 companions
    And are these companions SF or PE?
    '''

    groups = np.unique(neib[:,13])
    print('There are ', len(groups), 'groups')
    i = 0
    array = 0
    
    buckets = np.zeros((41,3))
    num = np.arange(0,41,1)
    array = np.c_[num, buckets]
    
    
    #### Final array
    # (0) Number of neighbours
    # (1) ALL= number of galaxies that have col 0 number of neighbours
    # (2) SF = Of those galaxies that have col 0 number of neighbours how many are SF
    # (3) PE = Of those galaxies that have col 0 number of neighbours how many are PE
    
    ### For each group of galaxies: A group represents the central+all found neighbouring galaxies in the defined radius: mrad
    while(i<len(groups)):
        print('Group', i)
        
        group = neib[neib[:,13]==i]
        nbgroup = group[group[:,12]<1] #Neighbours are ONLY those that are not centrals
        
        ## For testing purposes
        
        print('Number of NBs', len(nbgroup))
       
        
        ## ALL
        array[len(nbgroup)][1] = array[len(nbgroup)][1]+1
        
        ## SF and PE
        sfgroup = nbgroup[nbgroup[:,8]==1]
        pegroup = nbgroup[nbgroup[:,8]==0]
        
        array[len(nbgroup)][2] = array[len(nbgroup)][2] + len(sfgroup) ## SF
        array[len(nbgroup)][3] = array[len(nbgroup)][3] + len(pegroup) ## PE
        
        
        i = i + 1
        
        
    print(array)
 
   
    ### Express SF and PE in %%
    total = array[:,0]*array[:,1] #these numbers will set the 100% of the stacked bar
    
    psf = array[:,2]/total
    npsf = psf*array[:,1]
    
    
    ppe = array[:,3]/total
    nppe = ppe*array[:,1]
    
    
    #print('Total # of gals', total)
    #print('Percentage SF', psf)
    #print('What that percentage will be in #', npsf)
    
    
    ### Bar graph
    
    wbin = 0.7
    plt.figure(figsize=(10,7))
    #plt.xlim(-1,20)
    plt.xlim(7,23)
    #plt.ylim(0,2.5) #number of centrals
    plt.ylim(0,max(array[:,1])+1)    
    
    plt.xticks(np.arange(8,23,1))
    plt.ylabel('# of centrals', size=22)
    plt.xlabel('# of neighbours', size=22)
    plt.title('Given number of neighbours for centrals, Deep Fields', size=23)
    
    plt.bar(array[:,0]-wbin/2, array[:,1], width=wbin, color='green', alpha=0.5, bottom=0)
    
    plt.bar(array[:,0]-wbin/2, npsf, width=wbin, color='blue', alpha=0.5, bottom=0)
    plt.bar(array[:,0]-wbin/2, nppe, width=wbin, color='red', alpha=0.5, bottom=npsf)
    
    
    plt.text(array[8][0]-0.4,1.2*(array[8][1]),'SF '+str(int(psf[8]*100))+' %', color='blue')
    plt.text(array[8][0]-0.4,1.4*(array[8][1]),'PE '+str(int(ppe[8]*100))+' %', color='red')
    
    plt.text(array[12][0]-1,1.2*(array[12][1]),'SF '+str(int(psf[12]*100))+' %', color='blue')
    plt.text(array[12][0]-1,1.4*(array[12][1]),'PE '+str(int(ppe[12]*100))+' %', color='red')
    
    plt.text(array[14][0]-1,1.2*(array[14][1]),'SF '+str(int(psf[14]*100))+' %', color='blue')
    plt.text(array[14][0]-1,1.4*(array[14][1]),'PE '+str(int(ppe[14]*100))+' %', color='red')
    
    
    plt.text(array[15][0]-0.4,1.1*(array[15][1]),'SF '+str(int(psf[15]*100))+' %', color='blue')
    plt.text(array[15][0]-0.4,1.2*(array[15][1]),'PE '+str(int(ppe[15]*100))+' %', color='red')
    
    plt.text(array[18][0]-0.4,1.1*(array[18][1]),'SF '+str(int(psf[18]*100))+' %', color='blue')
    plt.text(array[18][0]-0.4,1.2*(array[18][1]),'PE '+str(int(ppe[18]*100))+' %', color='red')
    
    
    plt.text(array[22][0]-0.7,1.2*(array[22][1]),'SF '+str(int(psf[22]*100))+' %', color='blue')
    plt.text(array[22][0]-0.7,1.4*(array[22][1]),'PE '+str(int(ppe[22]*100))+' %', color='red')
    
    
    
    
###############################################################################
  
def func1(x, a, b):
    return a*x + b
    
    
###############################################################################
  
def func2(x, a, b, c):
  return a*(x**2) + b*x + c
  

###############################################################################
  
def func3(x, a, b, c):
  return a*(x**2) + b*x + c
    

############################################################################### 
  
def func4(x, a, b, c, d, e):
  return a*(x**4) + b*(x**3) + c*(x**2) + d*x + e
  
############################################################################### 

def scatter1(neib, ncentrals):
    

    centrals1 = neib[neib[:,12]==1] #ID the centrals of all groups
    #print(centrals1[:,4], centrals1[:,5])
    nbs1 = neib[neib[:,12]==0] # ID the neighbouring galaxies
    xx = np.arange(-26,-22,1)
    xx50 = (np.log10(0.5)+0.4*xx)/0.4
    #print(len(centrals1), len(nbs1))
    
    plt.figure(figsize=(10,7))
    plt.title('Correlation between centrals and their neighbours', size=23)
    plt.xlabel(r'Z$_{rest-frame}$ of centrals', size=21)
    plt.ylabel(r'Z$_{rest-frame}$ of neighbours', size=21)
    
    if(work_with=='wide'):
        plt.xlim(-26,-24.5)
        plt.ylim(-26,-23.5)
        
    elif(work_with=='deep'):
        #plt.xlim(-25.5,-24.8)
        plt.ylim(-25.5,-21.5)
        plt.xlim(min(centrals1[:,11])-0.1, max(centrals1[:,11])+0.1)
        
        
    plt.plot(xx, xx, color='black', linestyle=':', label='1:1', linewidth=3)
    plt.plot(xx, xx50, color='gray', linestyle='--', label='50%')
    plt.plot(xx50, xx, color='gray', linestyle='--')
    
    plt.fill_between(xx, xx50, xx, alpha=0.1, color='gray')
    plt.fill_between(xx50, xx50, xx, alpha=0.1, color='gray')
    plt.axvspan(min(centrals1[:,11]), max(centrals1[:,11]), color='y', alpha=0.2, lw=0)
    
    i = 0
    flag = 'True'
    
    while(i<=ncentrals):
        
        #print('Working with group', i)
        group = neib[neib[:,13]==i] # isolate one group at a time
        #print('Group', i, group)
        
        cent = group[group[:,12]==1] # ID the central of this group
        nb = group[group[:,12]==0] # ID the nbs
        #print('My central is', cent[:,6], cent[:,8], cent[:,12], cent[:,11])
        #print('And I have ', len(nb), 'nbs')
        centz = cent[:,11]
        centkse = cent[:,7]
        
        
        #plt.errorbar(centz, centz, yerr=centkse, xerr=centkse, fmt='o', ecolor='magenta', mec='magenta', mfc='None', mew=0.5, ms=8, alpha=1)
        plt.plot(centz, centz, linestyle='None', marker='o',  mec='magenta', mfc='None', mew=0.5, ms=8, alpha=1)
        
        if(len(nb)>0):
            
            for el in nb:
                col = 'black'
                #print('I have nb', el[6], el[8], el[12])
                
                if(el[8]==0):
                    col = 'red'
                elif(el[8]==1):
                    col = 'blue'
                    
                #plt.scatter(centz, el[11], color=col, marker='D')
                    
                if(flag=='True'):
                    plt.errorbar(centz, el[11], yerr=el[7], xerr=centkse, color=col, marker='D', label='Neighbour')
                    flag='False'
                    
                else:
                    plt.errorbar(centz, el[11], yerr=el[7], xerr=centkse, color=col, marker='D')
        
        elif(len(nb)==0): #for loners
             
            plt.scatter(centz, centz, color='magenta', marker='o', s=40)
            #plt.plot(centz, centz, linestyle='None', markeredgecolor='black', markerfacecolor='None', marker='o')
        
        i = i + 1
        
    plt.text(-25.9,-24, 'SF Neighbours', color='blue', fontsize=15)
    plt.text(-25.9,-24.25, 'PE Neighbours', color='red', fontsize=15)
    plt.text(-24.9,-25.3, 'Centrals', color='magenta', fontsize=15)
    
    plt.legend(numpoints=1)
    
############################################################################### 
    
def scatter2(neib, ncentrals):

    centrals = neib[neib[:,12]==1] #ID centrals
    print('CENTRALS', centrals.shape, min(centrals[:,16]), max(centrals[:,16]))
    kscentral =  centrals[:,6]
    
    plt.figure(figsize=(10,8))
    plt.ylim(-0.1,1.2)
    
    #plt.title('Correlation between centrals and their neighbours', size=23)
    plt.xlabel(r'i$_{rest-frame}$ of centrals', size=25)
    plt.ylabel(r'L/L$_c$', size=25)
    plt.xticks(np.arange(-25.8,-24.8,0.2), np.arange(-25.8,-24.8,0.2), fontsize=20)
    plt.yticks(np.arange(0,1.4,0.2), np.arange(0,1.4,0.2), fontsize=20)
    
    
    xx = np.arange(-26,-24, 0.01)
    yy = 0.5*np.ones(len(xx))
    yy2 = 1.5*np.ones(len(xx))
    
    if(work_with=='wide'):
        cl = 20.5
    elif(work_with=='deep'):
        cl = 23.0
        
    yy3 = 10**(-0.4*(cl-kscentral))
    limline  = np.c_[kscentral, centrals[:,16], yy3]
    limline=limline[np.argsort(limline[:,1])]
    
    
    plt.plot(xx, xx/xx, linestyle=':', color='black', linewidth=3, label='1:1')
    plt.plot(xx, yy, linestyle='--', color='gray', linewidth=2, label='50%')
    plt.plot(xx, yy2, linestyle=':', color='gray', linewidth=2)
    plt.fill_between(xx, yy, yy2, color='gray', alpha=0.2)
    
    plt.plot(centrals[:,16], np.ones(len(centrals[:,16])), linestyle='None', marker='o',  mec='magenta', mfc='None', mew=2, ms=8, alpha=1)
    plt.plot(limline[:,1], limline[:,2], linestyle='-.', color='purple', linewidth=2, label='Completeness', alpha=0.6)
    plt.fill_between(limline[:,1], -1*np.ones(len(limline[:,1])), limline[:,2], color='pink', alpha=0.4, hatch='x')
    
    plt.axvspan(min(centrals[:,16]), max(centrals[:,16]), color='y', alpha=0.2, lw=0)
   
    
    if(work_with=='deep'):
        plt.xlim(-25.7, -24.78)
        #plt.xlim(min(centrals[:,11])-0.22, max(centrals[:,11])+0.1)
        
    elif(work_with=='wide'):
        plt.xlim(-25.7, -24.78)
        #plt.xlim(min(centrals[:,11])-0.1, max(centrals[:,11])+0.1)
    
    
    if(work_with=='deep'):
        i = 1 #UMPEG numbers in the Deep fields start with 1 (0 was removed from file because is a double core)
        goto = ncentrals
    else:
        i=0
        goto = ncentrals-1
    
    
    flag = 'True'
    
    plt.errorbar(0,0, xerr=0, yerr=0, mec='black', mfc='None', color='None', label='Neighbours', marker='D', ecolor='black')
    
    while(i<=goto):
        print('Working with group', i)
        
        group = neib[neib[:,13]==i] #Group that holds that same number
        #print('Shape of this group', group.shape)
        nb2 = group[group[:,12]==0]
        central2 = group[group[:,12]==1] #ID central UMPEG
        centm = central2[:,11]
        centme = central2[:,7]
        
        centmi = central2[:,16] #Find rest-frame i mag for central UMPEG
        #print(centm[0], centme[0])
        
        if(len(nb2)>0):
            print('Number nbs', len(nb2))
            for el in nb2:
                rel = 10**(-0.4*(el[16]-centmi)) # xyz el[11] restz based only on mag notredz of central, el[15] restz based on redz of UMPEG, 16 resti based on rest of UMPEG
                if(rel>1):
                    print('Check', rel, el)
                
                ## Option 1
                #rele = np.sqrt(((np.log(10)*10**(0.4*(el[11]-centm)))**2)*(el[7])**2+((np.log(10)*10**(0.4*(el[11]-centm)))**2)*(centme)**2)
                
                ## Option 2
                ffc = 10**(-0.4*centmi)
                ffs = 10**(-0.4*el[15]) ### xyz
                sffc = ffc*centme*(np.log(10)/2.5) 
                sffs = ffs*el[7]*(np.log(10)/2.5)
                
                rele = np.sqrt((((1/ffc)**2)*(sffs**2))+(((ffs/ffc**2)**2)*(sffc**2)))             
                 
                if(el[8]==0):
                    col = 'red'
                elif(el[8]==1):
                    col = 'blue'
                
                plt.scatter(centmi, rel, color=col, marker='D')
                plt.errorbar(centmi, rel, yerr=rele, xerr=centme, color=col, ms=7, marker='D', ecolor=col)
                '''
                if(flag=='True'):
                    plt.errorbar(centm, rel, yerr=rele, xerr=centme, color=col, ms=7, marker='D', ecolor=col, label='Neighbours')
                    flag = 'False'
                else:
                    plt.errorbar(centm, rel, yerr=rele, xerr=centme, color=col, ms=7, marker='D', ecolor=col)
                '''
        elif(len(nb2)==0):
             
             plt.scatter(centmi, np.ones(len(centmi)), color='magenta', marker='o', s=30, lw=2)
        
        i = i + 1
    
    if(work_with=='wide'):
        plt.text(-25.58,0.8, 'SF Neighbours', color='blue', fontsize=15)
        plt.text(-25.58,0.7, 'PE Neighbours', color='red', fontsize=15)
        plt.text(-25.58,1.1, 'UMPEGs Wide Fields', color='magenta', fontsize=15)
        plt.legend(numpoints=1, loc='lower right', fontsize=17)
    elif(work_with=='deep'):
        plt.text(-25.58,0.8, 'SF Neighbours', color='blue', fontsize=15)
        plt.text(-25.58,0.7, 'PE Neighbours', color='red', fontsize=15)
        plt.text(-25.58,1.1, 'UMPEGs Deep Fields', color='magenta', fontsize=15)
        plt.legend(numpoints=1, loc='lower left', fontsize=17)
    
############################################################################### 
    
def cutouts(neib, ncentrals, mrad):
    
    '''
    only needs to be run once, will create the .cl file necessary to create cut-outs
    of your UMPEGs
    '''
    centr = neib[neib[:,12]==1]
    rad = mrad*1000/8.471/0.187 # from Mpc to pix at z~1.6
    outfile = open('test.cl','w')
    
    for item in centr:
        
        fieldt = str(int(item[0]))
        f = fieldt[1:2]
        x = item[2]
        y = item[3]
        group = item[13]
        
        xmin = int(x-rad)
        xmax = int(x+rad+1)
        
        ymin = int(y-rad)
        ymax = int(y+rad+1)
        
        if(xmin<0):
            xmin = 1
        if(ymin<0):
            ymin = 1
            
        if(xmax>19354):
            xmax = 19354
        if(ymax>19354):
            ymax = 19354
        
        
        ksu = item[6]
        
        if(f=='1'):
            patch = str(fieldt[2:8])+'-'+str(fieldt[8:])
        elif(f=='4'):
            if(str(fieldt[8:])=='003100'):
                patch = str(fieldt[2:8])+'-'+str(fieldt[8:])
            else:
                patch = str(fieldt[2:8])+'+'+str(fieldt[8:])
            
        print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+f+'/'+patch+'/CFHTLS_W_Ks_'+patch+'_T0007.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+f+'_'+patch+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_ks.fits', file=outfile)
        print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+f+'/'+patch+'/CFHTLS_W_u_'+patch+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+f+'_'+patch+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_u.fits', file=outfile)
        print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+f+'/'+patch+'/CFHTLS_W_g_'+patch+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+f+'_'+patch+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_g.fits', file=outfile)
        print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+f+'/'+patch+'/CFHTLS_W_r_'+patch+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+f+'_'+patch+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_r.fits', file=outfile)
        print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+f+'/'+patch+'/CFHTLS_W_z_'+patch+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+f+'_'+patch+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_z.fits', file=outfile)
    
    
    outfile.close()
    
############################################################################### 
    
def cutouts_deep(neib, ncentrals, mrad):
    
    '''
    only needs to be run once, will create the .cl file necessary to create cut-outs
    of your UMPEGs
    '''
    neib = neib[neib[:,0]==12] #Just to do it for D2
    centr = neib[neib[:,12]==1]
    rad = mrad*1000/8.471/0.187 # from Mpc to pix at z~1.6
    outfile = open('test.cl','w')
    
    for item in centr:
        
        fieldt = item[0]
        x = item[2]
        y = item[3]
        group = item[13]
        f = 0
        fo = 0
        
        xmin = int(x-rad)
        xmax = int(x+rad+1)
        
        ymin = int(y-rad)
        ymax = int(y+rad+1)
        
        if(xmin<0):
            xmin = 1
        if(ymin<0):
            ymin = 1
            
        if(xmax>19354):
            xmax = 19354
        if(ymax>19354):
            ymax = 19354
        
        
        ksu = item[6]
        
        
        
        if(fieldt==11):
            f = '022559-042940'
            fo = '1'
        elif(fieldt==12):
            f = '100028+021230'
            fo = '2'
            
        elif(fieldt==13):
            f = '141927+524056'
            fo = '3'
            
        elif(fieldt==14):
            f = '221531-174356'
            fo = '4'
            
            
        print('imcopy /Users/osejo/Desktop/Tesis/gzHK/photometry/wmask/WIRDS_Ks_'+f+'_T0002.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+fo+'_'+f+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_ks.fits', file=outfile)
        
        if(fo=='1')or(fo=='3')or(fo=='4')or(fo=='2'):        
            print('imcopy /Volumes/Liz/deep_images/CFHTLS_D-85_g_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+fo+'_'+f+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_g.fits', file=outfile)
            print('imcopy /Volumes/Liz/deep_images/CFHTLS_D-85_z_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+fo+'_'+f+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_z.fits', file=outfile)
            
        #elif(fo=='2'):
            #print('imcopy /Volumes/Liz/deep_images/D2.G.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+fo+'_'+f+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_g.fits', file=outfile)
            #print('imcopy /Volumes/Liz/deep_images/D2.Z.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+fo+'_'+f+'_umpeg_'+str(round(ksu,3))+'_group_'+str(group)+'_z.fits', file=outfile)
    
    
    outfile.close()
    
    
###############################################################################
    
def what_does_back_looklike(bcgo, bcgs, bcgp, neib, binsize):
    
    print('Range of z mags for neib ', min(neib[:,11]), max(neib[:,11]))
    
    testx = np.arange(min(neib[:,11]), max(neib[:,11]), binsize)
    
    bga = func4(testx, bcgo[0][0], bcgo[0][1], bcgo[0][2], bcgo[0][3], bcgo[0][4])
    bgsf = func4(testx, bcgs[0][0], bcgs[0][1], bcgs[0][2], bcgs[0][3], bcgs[0][4])    
    bgpe = func4(testx, bcgp[0][0], bcgp[0][1], bcgp[0][2], bcgp[0][3], bcgp[0][4])
    
    plt.figure()
    plt.yscale('log')
    plt.plot(testx, bga, linestyle='--', color='green')
    plt.plot(testx, bgsf, linestyle='--', color='blue')
    plt.plot(testx, bgpe, linestyle='--', color='red')
    
    
############################################################################### 
    
def restz_umpeg_redz(neibo, ncentralso):
    
    '''
    Define rest-frame Z (or i) magnitdes with respect of redshift of central UMPEG
    instead of muzzin's rainbow plot defined based on the mag of each central
    Created 14 Nov 2017
    '''
    umpeg_catalog_wide = neibo
    outfile = open('test.dat','w')
    ssp = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/sed_fit/make_bbsed/ssp_lr62_bc2003_z001.bbsed')
    agepe = 9.052783
    ebvpe = 0
    
    ssppe = ssp[ssp[:,2]==ebvpe]#Models that more closely reprsent  our PE population
    ssppe = ssppe[ssppe[:,4]==agepe]
    restfzp = ssppe[ssppe[:,0]==0][:,12]
   
    cons = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/sed_fit/make_bbsed/cons_lr62_bc2003_z001.bbsed')
    agesf = 8.324639
    ebvsf = 0.3
    
    conssf = cons[cons[:,2]==ebvsf]#Models that more closely reprsent  our PE population
    conssf = conssf[conssf[:,4]==agesf]
    restfzs = conssf[conssf[:,0]==0][:,12] #Rest-frame Z magnitude for this model
    
    array = []
    abzt = []
    
    
    zmagrest = -99
    
    num_groups = np.unique(umpeg_catalog_wide[:,13])
    
    for gnumber in num_groups:
        group = umpeg_catalog_wide[umpeg_catalog_wide[:,13]==gnumber]
        central_umpeg = group[group[:,12]==1]
        
        central_redz = central_umpeg[0][10]
        #print('For group', group[0][13], 'Central mag and id', central_umpeg[0][6], central_umpeg[0][0], central_umpeg[0][1], 'Reddz', central_redz)
        
        for gal in group:
            
            kind = gal[8]
            ks = gal[6]
            redz = central_redz#gal[10]#central_redz
            cat = 0
            #print(kind, ks, redz)
            
            if(kind==0):
                cat = ssppe
                zmagrest = restfzp
            elif(kind==1):
                cat = conssf
                zmagrest = restfzs
        
            
            for item in cat:
                
                if(item[0]==round(redz, 2)):
                    
                    kmag = item[15] #K magnitude for this model at this peak redshift
                    dmag = kmag - zmagrest #Delta Mag = Kmag model - Restframez
                    abz = ks - dmag # Absolute z mag = ksmag - deltamag
                    #abzt.append(abz)
                    print(redz, kind, item[0], ks, zmagrest, kmag, dmag, abz) ##To check output
                    print(gal[0], gal[1], gal[2], gal[3], gal[4], gal[5], gal[6], gal[7], gal[8], gal[9], gal[10], gal[11], gal[12], gal[13], gal[14], abz[0], file=outfile)
        
    #array = np.c_[umpeg_catalog_wide,  abzt]
    #np.savetxt('test2.dat', array, fmt='%s')
    outfile.close()
        
    
    
    return(0,0)
    
############################################################################### 
    

    
############################################################################### 

  
maglim = 19.5   #Only consider environments of galaxies brighter than this magnitude
artifact = 18.5 #Anything brighter than this will be considered an artifact and will not be included in the analysis
mrad = 0.5 # radius of search (proper radius) around desired object [Mpc] (At z~1.6 0.0166 degrees equals 0.5 proper Mpc ->Also radius used in Tal et al. 2012
work_with = 'wide' #wide, deep or all 
binsize = 0.001


binfrac = 0.1
binz = 0.1
iterations = 10
reg = open('test.reg', 'w')
    
    
print('# Region file format: DS9 version 4.1',file=reg)
print('# Environment around PEGs galaxies (Environment is comprised of SF and PE gzK galaxies', file=reg)
print('global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',file=reg)
print('fk5',file=reg)

#make_cats(work_with, maglim, mrad) #Only needs to be run once, it will determine who fullcat is


if(work_with=='wide'):
    fullcat = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/sfpe_w_cd_nostr_zabs_ks205.dat')
    fullcat = fullcat[fullcat[:,6]>artifact] #Remove possible bright artifacts, for now both SF and PE    
    fullcat = fullcat[(fullcat[:,12]<-24)] # Cut the Rest-frame Z mags to the completeness of the Passive Sample (Due to the diff redshifts btw PEGs and SFGs the Ks<20.5 was not enough, when transformed to Z rest frame, the SF population has a brighter mag completeness at 20.5 than the PEGs->See chapter 3 for details on the redz distribution)
    
    sfcat = fullcat[fullcat[:,8]==1]
    pecat = fullcat[fullcat[:,8]==0]
    
    #print(max(sfcat[:,12]), max(pecat[:,12])) #Different cuts due to different peaks in redz distribution, for testing purposes
    
    binsf = np.arange(-28,-24.0+binsize/2,binsize) #fullcat 
    binsp = np.arange(-28,-24.0+binsize/2,binsize) # PEGs--width of bin must be smaller than Tal12 reported gap btw populations, also it needs to cover the last completeness bin, in this case -24, hence the +half the bin for the max limit
    binss = np.arange(-28,-24.0+binsize/2,binsize) # SFGs--width of bin must be smaller than Tal12 reported gap btw populations, also it needs to cover the last completeness bin, in this case -24, hence the +half the bin for the max limit    
    binnb = np.arange(-28,-24.0+binsize/2,binsize) # neighbours-- width of bin must be smaller than Tal12 reported gap btw populations
    
    
    ### ARE THESE FOLLOWING 3 STILL IN USE?
    #binf = np.arange(0.05,1.15,binfrac) 
    #binbg = np.arange(0.05,1.05,binfrac)
    #bin2 = np.arange(-28,-23,binz)
    
    areaa = 23333.23882974 #Area of this work_with in Mpc2
    
elif(work_with=='deep'):
    
    fullcat = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/sfpe_d_cd_nostr_zabs_ks23.dat')
    fullcat = fullcat[fullcat[:,6]>artifact] #Remove bright artifacts    
    fullcat = fullcat[(fullcat[:,12]<-21.63)] # Cut the Rest-frame Z mags to the completeness of the Passive Sample (Due to the diff redshifts btw PEGs and SFGs the Ks<20.5 was not enough, when transformed to Z rest frame, the SF population has a brighter mag completeness at 20.5 than the PEGs->See chapter 3 for details on the redz distribution)    
    
    sfcat = fullcat[fullcat[:,8]==1]
    pecat = fullcat[fullcat[:,8]==0]
    
    print(max(sfcat[:,12]), max(pecat[:,12]))#Different cuts due to different peaks in redz distribution, for testing purposes
    
    binsf = np.arange(-28,-20.0+binsize/2,binsize) #fullcat 
    binsp = np.arange(-28,-20.0+binsize/2,binsize) # PEGs--width of bin must be smaller than Tal12 reported gap btw populations, also it needs to cover the last completeness bin, in this case -24, hence the +half the bin for the max limit
    binss = np.arange(-28,-20.0+binsize/2,binsize) # SFGs--width of bin must be smaller than Tal12 reported gap btw populations, also it needs to cover the last completeness bin, in this case -24, hence the +half the bin for the max limit    
    binnb = np.arange(-28,-20.0+binsize/2,binsize) # neighbours-- width of bin must be smaller than Tal12 reported gap btw populations
    
   
    areaa = 2334.25386459 #Area of this work_with in Mpc2
    
#elif(work_with=='all'):
#    fullcat = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/sfpe_dw_cd_nostr_zabs_ks205.dat')
#    fullcat = fullcat[(fullcat[:,12]<-24)] # Cut the Rest-frame Z mags to the completeness of the Passive Sample (Due to the diff redshifts btw PEGs and SFGs the Ks<20.5 was not enough, when transformed to Z rest frame, the SF population has a brighter mag completeness at 20.5 than the PEGs->See chapter 3 for details on the redz distribution)   
#    binsp = np.arange(-28,-24.0+binsize/2,binsize) #width of bin must be smaller than Tal12 reported gap btw populations
#    binnb = np.arange(-28,-24.0+binsize/2,binsize) #width of bin must be smaller than Tal12 reported gap btw populations
#    areaa = 25667.49269434 #Area of this work_with in Mpc2


###############################################################################
### LF of SF+PEGs in the Field (Up to Ks < 20.5) or 23 for deep


centerff = (binsf[:-1] + binsf[1:]) / 2 #Use for plotting full cat
centers = (binss[:-1] + binss[1:]) / 2 #Use for plotting Usually same but just in case things change will leave them independent SF
centerp = (binsp[:-1] + binsp[1:]) / 2 #Use for plotting PE


## Combining SF+PE
hall, halledges = np.histogram(fullcat[:,12], binsf, weights=fullcat[:,9]/fullcat[:,10]) #Histogram of background/field counts
halle, halledgese = np.histogram(fullcat[:,12], binsf, weights=fullcat[:,9]) #histogram of background/field not per area, since these will be used for the SQRT(N) errors


## Only SF
hsf, hsfedges = np.histogram(sfcat[:,12], binss, weights=sfcat[:,9]/sfcat[:,10]) #Histogram of background/field counts
hsfe, hsfedgese = np.histogram(sfcat[:,12], binss, weights=sfcat[:,9]) #histogram of background/field not per area, since these will be used for the SQRT(N) errors


## Only PE
hpe, hpeedges = np.histogram(pecat[:,12], binsp, weights=pecat[:,9]/pecat[:,10]) #Histogram of background/field counts
hpee, hpeedgese = np.histogram(pecat[:,12], binsp, weights=pecat[:,9]) #histogram of background/field not per area, since these will be used for the SQRT(N) errors


### Remove empty bins
aall = np.c_[centerff, hall, halle]
asf = np.c_[centers, hsf, hsfe]
ape = np.c_[centerp, hpe, hpee]


aall0 = aall[(aall[:,1]>0) & (aall[:,0]>-25.7)] #remove empty bins and do not include weird point into fit
asf0 = asf[(asf[:,1]>0) & (asf[:,0]>-25.7)] #remove empty bins and do not include weird point into fit
ape0 = ape[(ape[:,1]>0) & (ape[:,0]>-25.7)] #remove empty bins and do not include weird point into fit


### Find a fit for this background correction
bcgo = curve_fit(func4, aall0[:,0], aall0[:,1], sigma=aall0[:,2])
a = bcgo[0][0]
b = bcgo[0][1]
c = bcgo[0][2]
d = bcgo[0][3]
e = bcgo[0][4]

bcgs = curve_fit(func4, asf0[:,0], asf0[:,1], sigma=asf0[:,2])
asff = bcgs[0][0]
bsff = bcgs[0][1]
csff = bcgs[0][2]
dsff = bcgs[0][3]
esff = bcgs[0][4]


bcgp = curve_fit(func4, ape0[:,0], ape0[:,1], sigma=ape0[:,2])
apef = bcgp[0][0]
bpef = bcgp[0][1]
cpef = bcgp[0][2]
dpef = bcgp[0][3]
epef = bcgp[0][4]




if(work_with=='wide'):
    x = np.arange(-27.7,-23.9,0.1) 

elif(work_with=='deep'):
    x = np.arange(-27.7,-23.9,0.1)
    
    
undall = np.sqrt(aall0[:,2])/(np.ones(len(aall0[:,2]))*areaa) #Size of errorbars in poissonian errors SQRT(N)/Area, make sure area is NOT inside of SQRT when doing this
undsf = np.sqrt(asf0[:,2])/(np.ones(len(asf0[:,2]))*areaa) #They all come from the same area
undpe = np.sqrt(ape0[:,2])/(np.ones(len(ape0[:,2]))*areaa)


### Number counts per Mpc2 per mag bin of the background/field (PE+SFGs) in rest-frame Z magnitude bins

'''
plt.figure(figsize=(10,7))
plt.yscale('log') #to make sure the size of the errorbars is correct I will do everything in log space
plt.title('Field surface luminosity fuction of SF and PEGs in '+work_with, fontsize=21)
plt.ylabel('$\Sigma$ (#/$Mpc^2$/'+str(binsize)+' Zmag bin)', fontsize=19)
plt.xlabel('Z [Rest-frame Mag]', fontsize=19)


 
if(work_with=='deep'):
    plt.xlim(-26,-21)
    plt.ylim(10**-4, 1.8*10**0)
    
elif(work_with=='wide'):
    plt.xlim(-26,-24)
    plt.ylim(2*10**-5, 0.2)


##SF+PE
plt.errorbar(aall0[:,0], aall0[:,1], yerr = undall, marker='d', color='green', linestyle='None', markersize=10, alpha=0.5, label='SF+PE')
plt.plot(aall0[:,0], (a*(aall0[:,0]**4)) + (b*aall0[:,0]**3) + (c*aall0[:,0]**2) + (d*aall0[:,0]) + e, color='green', linestyle='--')

##SF
plt.errorbar(asf0[:,0], asf0[:,1], yerr = undsf, marker='d', color='blue', linestyle='None', markersize=10, alpha=0.5, label='SF')
plt.plot(asf0[:,0], (asff*(asf0[:,0]**4)) + (bsff*asf0[:,0]**3) + (csff*asf0[:,0]**2) + (dsff*asf0[:,0]) + esff, color='blue', linestyle='--')

##PE
plt.errorbar(ape0[:,0], ape0[:,1], yerr = undpe, marker='d', color='red', linestyle='None', markersize=10, alpha=0.5, label='PE')
plt.plot(ape0[:,0], (apef*(ape0[:,0]**4)) + (bpef*ape0[:,0]**3) + (cpef*ape0[:,0]**2) + (dpef*ape0[:,0]) + epef, color='red', linestyle='--')

plt.legend(loc='lower right', numpoints=1)
'''


##############################################################################
## Find neighbours around centrals

if(work_with=='deep'):
    neib = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/0.5Mpc/deep_0.5_centrals_18.5-19.5_curve_zirest_umpeg.dat')
    ncentrals = len(neib[neib[:,12]==1])
    
else:
    neib_old, ncentrals = select_nb(fullcat, maglim, mrad, work_with)
    #neib_tosave, ncentrals = restz_umpeg_redz() #At that stage save on to file, only needs to be run once: NO LONGER USED go to /Tal12/1116/extras.py resti_umpeg_redz()
    neib = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/0.5Mpc/wide_0.5_centrals_18.5-19.5_zirest_umpeg.dat')

############################################################################################################
## Check background correction

#what_does_back_looklike(bcgo, bcgs, bcgp, neib, binsize)

##############################################################################

###############################################################################
##### Fractional Graph
### These fractional plots are no longer in use: For Wide use fractional_background.py and for Deep fractional_background_deep.py: Nov 15 2017
#fractional(neib, bcgo, bcgs, bcgp, mrad, binsf, binsp, binss, ncentrals)
#fractional2(neib, bcgo, bcgs, bcgp, mrad, binsf, binsp, binss, ncentrals)

##############################################################################

#### Kind of neighbourhood
#if(work_with=='wide'):
    #kind_of_neighborhood(neib, mrad)
#elif(work_with=='deep'):
    #kind_of_neighborhood_deep(neib, mrad)

##############################################################################


#scatter1(neib, ncentrals)

############################################################################### 

scatter2(neib, ncentrals)

############################################################################### 

#cutouts(neib, ncentrals, mrad)
#cutouts_deep(neib, ncentrals, mrad)

############################################################################### 

#np.savetxt('test.dat', neib, fmt='%s')

###Luminosity function of neighbours around centrals

'''
### LF of SF+PEGs around centrals FULL CAT

centernb = (binnb[:-1] + binnb[1:]) / 2 #Use for plotting


### Do a histogram based on rest-frame Z mags
hnb, hnbedges = np.histogram(neib[:,11], binnb, weights=neib[:,9]/neib[:,14]) #col 11 is where abs mags are, cols # are different thn in the fied
hnbe, hnbedgese = np.histogram(neib[:,11], binnb, weights=neib[:,9])


areasnb = np.unique(neib[:,14]) #All the neighbouing galaxies come from the same areas (which is just adding all the circles)
undnb = np.sqrt(hnbe)/(np.ones(len(hnbe))*areasnb) #poisonian errors, square root of all counts BEFORE they are normalized by area


plt.figure()
plt.yscale('log') #to make sure the errorbars are OK
plt.title('Neighbouring galaxies around PEGs Ks<'+str(maglim), fontsize=21)
plt.ylabel('$\Sigma$ (#/$Mpc^2$/'+str(binsize)+' Zmag bin)', fontsize=19)
plt.xlabel('Z [Rest-frame Mag]', fontsize=19)


plt.ylim(10**-5, 10)
plt.xlim(-28.3,-21.0)


### Not corrected for background:
plt.errorbar(centernb, hnb, yerr = undnb, marker='v', linestyle='None', markersize=10, color='pink', alpha = 0.7, label='Not corrected for background')
#plt.scatter(centernb, hnb, marker='v', color='pink', s=70)


### LF of SF+PEGs around centrals background corrected

hnbc = hnb-hall

undnbc = np.sqrt(hnbc*areasnb)/areasnb

plt.errorbar(centernb, hnbc, yerr = undnbc, marker='o', linestyle='None', markersize=10, color='purple', alpha=0.8, label='Background corrected')
plt.errorbar(centerff, hall, yerr = undall, marker='d', color='green', linestyle='None', markersize=10, alpha = 0.5, label='Background')
#plt.scatter(centernb, hnbc, marker='o', color='green', s=70)
plt.legend(loc='lower right', numpoints=1)

###Do the same as above but dividing the neighbours btw SF and PE

'''

reg.close()
