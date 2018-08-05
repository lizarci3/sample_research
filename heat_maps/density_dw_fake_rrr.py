# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 16:58:45 2017

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
import timeit
import os
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

start=timeit.time.time()


'''
To examine the significance of our heat maps, we try to reproduce high-density environments 
by random-random-random encounters of background galaxies
'''

def make_sure_atleast_three(cont_dense, i):
    
    '''
    Of all the overdense regions that were selected above level l
    Make sure each one has at least 3 fake galaxies inside-just like it was done with the real data
    cont_dense: dictionary
    '''
    
    fakecatn = np.genfromtxt('test_fakecat_dw_it'+str(i)+'.dat')
    
    for item in cont_dense:
        print('Checking which contour has at least 3 obj in ', item)
        
        #ninit = 0
        
        x = cont_dense[item]
        i = 0
        
        if(len(x.allsegs[1])==0):
            print('No group found')
            #ninit = 0 
            
            
        else:
            
            for i in range(len(x.allsegs[1])): # For all 100% contour paths
                j = 0
                
                cont = x.allsegs[1][i] #ith 100% contour level
            
                # Find its nearest contour, but only those that are on x% level (the x is user defined)
                contmin, segmin, indmin, xmin, ymin, dmin = x.find_nearest_contour(cont[0][0], cont[0][1], indices = [0], pixel = False)
                
                bbPath = mplPath.Path(x.allsegs[contmin][segmin])
                
                for fake in fakecatn:
                    
                     xa = fake[2]
                     ya = fake[3]
            
                     xga = (xa/scale)
                     yga = (ya/scale)
                     
                    
                     C = bbPath.contains_point((xga, yga)) # search if the  scaled position is within the x% level
        
                     if(C==True):
                         j = j+1
                if(j>=3):
                    
                    print('Found one optimal')
                    #ninit = ninit + 1
                    
        
                    
    #return(ninit)
###--------------------------------------------------------------------------------------------------------------------


def select_dense_environments(mdi, i):
    
    '''
    Having the final gaussian densiy map, first ID contours at certain, user defined, level
    And then check that they are actually complete (3 or more fake objects inside of limit*level)
    '''

    cont_dense = {}
    
    for each in mdi:
        print('Finding contours in', each)
        plt.figure(num=1)
        cont_dense[each] = plt.contour(mdi[each][int(gsmall.shape[0]/2):(mdi[each].shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(mdi[each].shape[0]-int(gsmall.shape[0]/2))].T, levels=[limgroup*level, level], cmap=plt.get_cmap('Blues'), alpha=1)
    
    make_sure_atleast_three(cont_dense, i)
    
    
    
###-----------------------------------------------------------------------------------------------------------------------
def load_deep_masks():
    
    '''
    Load the masks for the deep fields that will be used to define if a fake position is good or not.
    This is done to ensure the fake data has the same spatial density as the REAL data
    '''
    
    deep = ['d1', 'd2', 'd3', 'd4']
    dimask = {}
    
    for d in deep:
        maskdeep =  np.load('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_'+d+'.npy')
        dimask[d] = maskdeep.astype(bool)
        
    '''
    maskd2 =  np.load('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_d2.npy')
    maskd3 =  np.load('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_d3.npy')
    maskd4 =  np.load('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_d4.npy')
    
    dimask = {'d1':maskd1, 'd2':maskd2, 'd3':maskd3, 'd4':maskd4}
    '''
    return(dimask)
###--------------------------------------------------------------------------------------------------------------------

def load_wide_masks():
    
    '''
    Load the masks for the deep fields that will be used to define if a fake position is good or not.
    This is done to ensure the fake data has the same spatial density as the REAL data
    Laptop does not have enough RAM use Cedar
    '''
    
    wide = ['w1']
    dw1masks = {}
    
    for w in wide:
        names = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_'+w+'.txt', dtype=str)
        
        for each in names:
            tmpMask = np.load('/Volumes/Liz/masks_fake_env/'+w+'/test_'+w+'_'+each+'_gzkmask.npy')
            dw1masks[each] = tmpMask.astype(bool) 
            
    return(dw1masks)
            
###--------------------------------------------------------------------------------------------------------------------
    
def show_density_plot(slzerosd1, slzerosd2, slzerosd3, slzerosd4, slzerosw1, slzerosw4):
    
    fakecatn = np.genfromtxt('test_fakecat_dw_it'+str(i)+'.dat')
    
    '''
    #### Show Gaussian Density Plot
    plt.figure(num=1, figsize=(10,10), dpi=100)
    plt.xlim(0,200)
    plt.ylim(0,200)
    plt.imshow(slzerosd1[int(gsmall.shape[0]/2):(slzerosd1.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosd1.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()
    d1points = fakecatn[fakecatn[:,1]==1]
    plt.scatter(d1points[:,2]/scale, d1points[:,3]/scale, s=50, facecolor='None', edgecolor='red', marker='*')
    plt.scatter(d1points[:,2]/scale, d1points[:,3]/scale, s=25, facecolor='None', edgecolor='gold', marker='*')
    
    
    
    plt.figure(figsize=(10,10), dpi=100)
    plt.xlim(0,200)
    plt.ylim(0,200)
    plt.imshow(slzerosd2[int(gsmall.shape[0]/2):(slzerosd2.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosd2.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()

    plt.figure(figsize=(10,10), dpi=100)
    plt.xlim(0,200)
    plt.ylim(0,200)
    plt.imshow(slzerosd3[int(gsmall.shape[0]/2):(slzerosd3.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosd3.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()
    '''
    plt.figure(figsize=(10,10), dpi=100)
    plt.xlim(0,200)
    plt.ylim(0,200)
    plt.imshow(slzerosd4[int(gsmall.shape[0]/2):(slzerosd4.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosd4.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()
    d1points = fakecatn[fakecatn[:,1]==4]
    plt.scatter(d1points[:,2]/scale, d1points[:,3]/scale, s=50, facecolor='None', edgecolor='red', marker='*')
    plt.scatter(d1points[:,2]/scale, d1points[:,3]/scale, s=25, facecolor='None', edgecolor='gold', marker='*')
   
    
    '''
    plt.figure(figsize=(10,10), dpi=100)
    plt.imshow(slzerosw4[int(gsmall.shape[0]/2):(slzerosw4.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosw4.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()
    
    
    
    plt.figure(figsize=(18,10))
    plt.xlim(1700,0)
    plt.ylim(0,700)
    plt.imshow(slzerosw1[int(gsmall.shape[0]/2):(slzerosw1.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosw1.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()
    '''

###---------------------------------------------------------------------------------------------------------------------------------
      
def make_empty_arrays(scale):
    
    slzerosd1 = np.zeros((30000/scale, 30000/scale)) #The extra pixels secure a 5sigma extra padding for edge effects
    slzerosd2 = np.zeros((30000/scale, 30000/scale))
    slzerosd3 = np.zeros((30000/scale, 30000/scale))
    slzerosd4 = np.zeros((30000/scale, 30000/scale))
    slzerosw1 = np.zeros((300000/scale, 90000/scale))
    slzerosw4 = np.zeros((120000/scale, 93000/scale))
    
    return(slzerosd1, slzerosd2, slzerosd3, slzerosd4, slzerosw1, slzerosw4)
###---------------------------------------------------------------------------------------------------------------------------------
    
def shift(fieldw, patch):
    
    if(fieldw==1):
        field = str(patch[1:6]+patch[7:])
    
    else:
        field = str(patch[:6]+patch[7:])
        
    #print('Which would be', field)
    
    if(fieldw==1):
        xpatch = {'20241041200':8, '20241050800':8,'20241060400':8, '20631041200':7, '20631050800':7, '20631060400':7, '20631070000':7, '21021041200':6, '21021050800':6, '21021060400':6, '21021070000':6, '21410041200':5, '21410050800':5, '21410060400':5, '21800041200':4, '21800050800':4, '21800060400':4, '22150041200':3, '22150050800':3, '22150060400':3, '22539050800':2, '22539060400':2, '22929041200':1, '22929050800':1, '22929060400':1, '23319041200':0, '23319050800':0, '23319060400':0}
        ypatch = {'20241041200':3, '20241050800':2,'20241060400':1, '20631041200':3, '20631050800':2, '20631060400':1, '20631070000':0, '21021041200':3, '21021050800':2, '21021060400':1, '21021070000':0, '21410041200':3, '21410050800':2, '21410060400':1, '21800041200':3, '21800050800':2, '21800060400':1, '22150041200':3, '22150050800':2, '22150060400':1, '22539050800':2, '22539060400':1, '22929041200':3, '22929050800':2, '22929060400':1, '23319041200':3, '23319050800':2, '23319060400':1}

        shiftx = xpatch[(field)]
        shifty = ypatch[(field)]
    
    elif(fieldw==4):
        xpatch = {'220154011900':5, '220154021500':5,'220154031100':5, '220542011900':4, '220542021500':4, '220542031100':4, '220930003100':3, '220930002300':3, '220930011900':3, '220930021500':3, '221318003100':2, '221318002300':2, '221318011900':2, '221318021500':2, '221706002300':1, '221706011900':1, '221706021500':1, '222054002300':0, '222054011900':0}
        ypatch = {'220154011900':2, '220154021500':3,'220154031100':4, '220542011900':2, '220542021500':3, '220542031100':4, '220930003100':0, '220930002300':1, '220930011900':2, '220930021500':3, '221318003100':0, '221318002300':1, '221318011900':2, '221318021500':3, '221706002300':1, '221706011900':2, '221706021500':3, '222054002300':1, '222054011900':2}
    
        shiftx = xpatch[(field)]
        shifty = ypatch[(field)]
        
        
    return(shiftx, shifty)
    
###---------------------------------------------------------------------------------------------------------------------------

def makeGaussian2(size, sigma):
    
    if(size%2 != 0): #make it centered on a whole pixel
        size = size+1
    
    x=np.arange(size)-int(size/2)
    y=np.arange(size)-int(size/2)
    xx, yy = np.meshgrid(x,y)
    
    array = ((1/(2*np.pi*(sigma)**2))*np.exp(-((xx)**2 + (yy)**2)/(2*sigma**2)))
    arraym = (1/np.amax(array))*array
    return arraym, xx, yy
    
###----------------------------------------------------------------------------------------------------------------------
    
def find_fake_positions(minp, maxp, scale, field, f, patch, dimask, dwmask):
    
    '''
    Find the fake positions of Ks<20.5 PEGs taking into account
    that we do not place them in masked regions
    '''
    
    flag = False # As long as we can't be sure they are not in a flagged region repeat procedure
    
    while(flag==False):
        
        mask = []
        fakex = float(np.random.randint(minp, maxp))
        fakey = float(np.random.randint(minp, maxp))
        print('Chose random positions', fakex, fakey)
        
        if(field[:-1]=='d'): ##Check pixels in Deep fields
            
           mask = dimask[field]
           flag = (mask[fakex][fakey]==False) #bool masks have False where there were values of 0 in the pix (i.e., good pixs)
           print('Flag for this position', flag)
            
        else:
            
            #mask = dwmask[patch]
            flag = True#((mask[fakex][fakey]==False))
            print('Flag for this position', flag)
                
        
        if(flag==True):
            print('This is good point')
            return(fakex, fakey)
            
###--------------------------------------------------------------------------------------------------------------------

fields = ['d1','d2','d3','d4', 'w1', 'w4']
scale = 100 #to make things more manageable, scale everything by   
fwhm = 1000/scale # 632 pix or 0.032 degrees represents a box of side 1Mpc in size @ z~1.5 
minp = 1 
maxp = 19354 #For the deep fields
level = 2.7
limgroup = 70/100
its = 1

### Needs to be loaded only once as global variable
dimask = load_deep_masks()
#dwmask = load_wide_masks()
dwmask = 0


### Real catalog
cat = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_dw_cd_nostr_maglim_v_mu_2s01z05k_clean.dat')
#cat = cat[cat[:,6]<20.5] #Only Ks<20.5 define the existence of a dense environment

### Catalogs for testing purposes
cat = cat[0:200][:] #for testing: work with only the first 100 points in your data
#cat = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense_fake/test_full_cat.dat')
#cat = np.arange(1,5,1)


### Define names of patches in w1 and w4
patchesw1 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_w1.txt', dtype=str)
patchesw4 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_w4.txt', dtype=str)


#### Define sigma for each gaussian
sig = (round(fwhm/2.355,0))

#### Make 2D Gaussian
gsmall, xx, yy = makeGaussian2(5*sig, sig)
size = gsmall.shape[0]
print('Made gsmall! It has a shape of ', gsmall.shape)

          
i = 1

while(i<=its):
    print('ITERATION', i)
    slzerosd1, slzerosd2, slzerosd3, slzerosd4, slzerosw1, slzerosw4 = make_empty_arrays(scale) 
    fakecat = open('test_fakecat_dw_it'+str(i)+'.dat', 'w')
    print('# Fake objects placed in both Deep and Wide fields at random positions for iteration '+str(i), file=fakecat)
    print('# (0) Random Chosen Field', file = fakecat)
    print('# (1) Random Chosen Patch (name if on Wide, 0 if on Deep)', file = fakecat)
    print('# (2) X random original', file = fakecat)
    print('# (3) Y random original', file = fakecat)
    print('# (4) X random (scaled, shifted half a gaussian)', file = fakecat)
    print('# (5) Y random (scaled, shifted half a gaussian)', file = fakecat)
    print('# (6) Ks [AB mag]', file = fakecat)
    print('# (7) Ks error', file = fakecat)
    
    for peg in cat:
        
        pixxf = 0
        pixyf = 0
        overlapx = 968
        overlapy = 1290
        x,y,xg,yg = 0,0,0,0
        
        ####Choose a random field
        field = fields[np.random.randint(0,len(fields))]
        f = float(field[-1:])
        print('I chose field', field, f)
        
        if(field[:-1]!='d'): #If it didn't choose Deep Fields
            print('I selected wide', field)
            
            if(f==1):
                patch = patchesw1[np.random.randint(0,len(patchesw1))]
                print('I chose patch', patch)
            elif(f==4):
                patch = patchesw4[np.random.randint(0,len(patchesw4))]
                print('I chose patch', patch)
            
            ### If they are in the wide fields they need to be shifted to their position in the field
            if(f==4)and((patch=='220930-003100')or(patch=='221318+002300')):
                overlapy=1741
                
            shiftx, shifty = shift(f, patch)
                
        else:
            f = float(field[-1:]) # the deep fields don't have smaller patches
            patch = float(field[-1:])
                
        #### Now that it has chosen one of the deep fields or one patch in wide
        #### Find random position in that patch/field        
        
        
        ## The following function also takes into account masks and returns random positions
        ## POSITIONS HAVE NOT BEEN SCALED YET
        pixxf, pixyf = find_fake_positions(minp, maxp, scale, field, f, patch, dimask, dwmask) 
        
        
        ## When working with Wide fields:
        if(field[:-1]!='d'):
            
            # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
            x = pixxf+shiftx*(19354-overlapx)
            y = pixyf+shifty*(19354-overlapy)
            #print('Shifting in wide only', 'Shifted pos', x, y)
             
            xg = round((x/scale),0)+round((gsmall.shape[0]/2),0)
            yg = round((y/scale),0)+round((gsmall.shape[0]/2),0)
            #print('Shifted WIDE by half a gaussian', xg, yg)
            
            ## Save to fake catalog
            print(field, patch, pixxf, pixyf, xg, yg, peg[6], peg[7], file=fakecat)
            
        else: # When working with Deep fields
            
            xg = round((pixxf/scale),0)+int(gsmall.shape[0]/2) #if you don't want to round up x, y positions
            yg = round((pixyf/scale),0)+int(gsmall.shape[0]/2)
            #print('Shifted DEEP by half a gaussian', xg, yg)
            
            ## Save to fake catalog
            print(field, f, pixxf, pixyf, xg, yg, peg[6], peg[7], file=fakecat)
        
        ### Add this gaussian to fake field
        
        if(field=='d1'):
            slzerosd1[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosd1[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
        elif(field=='d2'):
            slzerosd2[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosd2[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
        elif(field=='d3'):
            slzerosd3[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosd3[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
        elif(field=='d4'):
            slzerosd4[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosd4[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
       
       
        elif(field=='w1'):
            slzerosw1[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosw1[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
       
        elif(field=='w4'):
            slzerosw4[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosw4[xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
    
    fakecat.close()
    
    
    ### For each iteration, after it has gone through all the catalog 
    ### Show the output density plots and select dense environments
    
    show_density_plot(slzerosd1, slzerosd2, slzerosd3, slzerosd4, slzerosw1, slzerosw4)
    
    mdi = {'d4':slzerosd4}#, 'd2':slzerosd2, 'd3':slzerosd3, 'd4':slzerosd4, 'w1':slzerosw1, 'w4':slzerosw4}
    select_dense_environments(mdi, i)
    
    
    
    i = i + 1
    
print('t_total=',timeit.time.time()-start, 's')
os.system('say "Liz python is done"')   