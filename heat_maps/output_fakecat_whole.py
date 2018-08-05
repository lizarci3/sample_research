# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:10:58 2017

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
    
def shift(fieldw, patch):
    
    
    if(fieldw=='w1'):
        xpatch = {'20241041200':8, '20241050800':8,'20241060400':8, '20631041200':7, '20631050800':7, '20631060400':7, '20631070000':7, '21021041200':6, '21021050800':6, '21021060400':6, '21021070000':6, '21410041200':5, '21410050800':5, '21410060400':5, '21800041200':4, '21800050800':4, '21800060400':4, '22150041200':3, '22150050800':3, '22150060400':3, '22539050800':2, '22539060400':2, '22929041200':1, '22929050800':1, '22929060400':1, '23319041200':0, '23319050800':0, '23319060400':0}
        ypatch = {'20241041200':3, '20241050800':2,'20241060400':1, '20631041200':3, '20631050800':2, '20631060400':1, '20631070000':0, '21021041200':3, '21021050800':2, '21021060400':1, '21021070000':0, '21410041200':3, '21410050800':2, '21410060400':1, '21800041200':3, '21800050800':2, '21800060400':1, '22150041200':3, '22150050800':2, '22150060400':1, '22539050800':2, '22539060400':1, '22929041200':3, '22929050800':2, '22929060400':1, '23319041200':3, '23319050800':2, '23319060400':1}

        shiftx = xpatch[str(patch)]
        shifty = ypatch[str(patch)]
    
    elif(fieldw=='w4'):
        xpatch = {'220154011900':5, '220154021500':5,'220154031100':5, '220542011900':4, '220542021500':4, '220542031100':4, '220930003100':3, '220930002300':3, '220930011900':3, '220930021500':3, '221318003100':2, '221318002300':2, '221318011900':2, '221318021500':2, '221706002300':1, '221706011900':1, '221706021500':1, '222054002300':0, '222054011900':0}
        ypatch = {'220154011900':2, '220154021500':3,'220154031100':4, '220542011900':2, '220542021500':3, '220542031100':4, '220930003100':0, '220930002300':1, '220930011900':2, '220930021500':3, '221318003100':0, '221318002300':1, '221318011900':2, '221318021500':3, '221706002300':1, '221706011900':2, '221706021500':3, '222054002300':1, '222054011900':2}
    
        shiftx = xpatch[str(patch)]
        shifty = ypatch[str(patch)]
        
        
    return(shiftx, shifty)
    
###---------------------------------------------------------------------------------------------------------------------------

      
def make_empty_arrays(scale):
    
    slzerosd1 = np.zeros((30000/scale, 30000/scale)) #The extra pixels secure a 5sigma extra padding for edge effects
    slzerosd2 = np.zeros((30000/scale, 30000/scale))
    slzerosd3 = np.zeros((30000/scale, 30000/scale))
    slzerosd4 = np.zeros((30000/scale, 30000/scale))
    slzerosw1 = np.zeros((300000/scale, 90000/scale))
    slzerosw4 = np.zeros((150000/scale, 98000/scale))
    
    slzerosfull = {'d1':slzerosd1, 'd2':slzerosd2, 'd3':slzerosd3, 'd4':slzerosd4, 'w1':slzerosw1, 'w4':slzerosw4}
    return(slzerosfull)
    
###---------------------------------------------------------------------------------------------------------------------------

def show_density_plot(slzerosfull, warray, i, field):
    
    fakecatn = np.genfromtxt('/Volumes/Liz/fake_env/output/whole_masks_v'+str(version)+'/'+field+'/test_fakecat_new7_'+field+'_cluster_it'+str(i)+'_v'+str(version)+'.dat')
    
    xmax = {'d1':200, 'd2':200, 'd3':200, 'd4':200, 'w1':1700, 'w4':1100}
    ymax = {'d1':200, 'd2':200, 'd3':200, 'd4':200, 'w1':700, 'w4':800}
    fignum =  {'d1':1, 'd2':2, 'd3':3, 'd4':4, 'w1':10, 'w4':40}
    size1 = {'d1':10, 'd2':10, 'd3':10, 'd4':10, 'w1':18, 'w4':10}
    size2 = {'d1':10, 'd2':10, 'd3':10, 'd4':10, 'w1':10, 'w4':10}
    
    #### Show Gaussian Density Plot
    plt.figure(num=fignum[field], figsize=(size1[field],size2[field]), dpi=100)
    plt.xlim(0,xmax[field])
    plt.ylim(0,ymax[field])
    plt.imshow(slzerosfull[field][int(gsmall.shape[0]/2):(slzerosfull[field].shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosfull[field].shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    plt.tight_layout()
    
    if(field[0:1]=='d'):
        plt.scatter(fakecatn[:,1]/scale, fakecatn[:,2]/scale, s=50, facecolor='None', edgecolor='red', marker='*')
        plt.scatter(fakecatn[:,1]/scale, fakecatn[:,2]/scale, s=25, facecolor='None', edgecolor='gold', marker='*')
        
    else:
        plt.scatter(warray[:,0]/scale, warray[:,1]/scale, s=50, facecolor='None', edgecolor='red', marker='*')
        plt.scatter(warray[:,0]/scale, warray[:,1]/scale, s=25, facecolor='None', edgecolor='gold', marker='*')
    
###---------------------------------------------------------------------------------------------------------------------------

def select_dense_environments(slzerosfull, i, field, warray, nt):
    
    '''
    Having the final gaussian densiy map, first ID contours at certain, user defined, level
    And then check that they are actually complete (3 or more fake objects inside of limit*level)
    '''
    fignum =  {'d1':1, 'd2':2, 'd3':3, 'd4':4, 'w1':10, 'w4':40}
    cont_dense = {}
    
        
    plt.figure(num=fignum[field])
    cont_dense[field] = plt.contour(slzerosfull[field][int(gsmall.shape[0]/2):(slzerosfull[field].shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzerosfull[field].shape[0]-int(gsmall.shape[0]/2))].T, levels=[limgroup*level, level], cmap=plt.get_cmap('Blues'), alpha=1)
    make_sure_atleast_three(cont_dense[field], i, field, warray, nt)
    
    
    
###-----------------------------------------------------------------------------------------------------------------------

def make_sure_atleast_three(cont_dense_field, i, field, warray, nt):
    
    '''
    Of all the overdense regions that were selected above level l
    Make sure each one has at least 3 fake galaxies inside-just like it was done with the real data
    cont_dense: dictionary entrance of just the field that is being worked with
    '''
    
    fakecatn = np.genfromtxt('/Volumes/Liz/fake_env/output/whole_masks_v'+str(version)+'/'+field+'/test_fakecat_new7_'+field+'_cluster_it'+str(i)+'_v'+str(version)+'.dat')
    print('Checking which contour has at least 3 obj in iteration', i)
    
    
    ninit = 0
    x = cont_dense_field
    
    
    if(len(x.allsegs[1])==0):
        print('No group found iteration', i)
        ninit = 0 
        
        
    else:
        
        for i2 in range(len(x.allsegs[1])): # For all 100% contour paths
            j = 0
            m = 0
            
            cont = x.allsegs[1][i2] #ith 100% contour level
        
            # Find its nearest contour, but only those that are on x% level (the x is user defined)
            contmin, segmin, indmin, xmin, ymin, dmin = x.find_nearest_contour(cont[0][0], cont[0][1], indices = [0], pixel = False)
            bbPath = mplPath.Path(x.allsegs[contmin][segmin])
            
            if(field[0:1]=='d'):
                
                for fake in fakecatn:
                    
                     xa = fake[1]
                     ya = fake[2]
            
                     xga = (xa/scale)
                     yga = (ya/scale)
                     
                    
                     C = bbPath.contains_point((xga, yga)) # search if the  scaled position is within the x% level
                     
        
                     if(C==True):
                         j = j+1
                if(j>=3):
                    
                    print('Found one optimal')
                    #ninit = ninit + 1
            else:
                
                for fake in warray:
                          
                    C2 = bbPath.contains_point((fake[0]/scale, fake[1]/scale)) # search if the  scaled position is within the x% level
        
                    if(C2==True):
                         m = m+1
                         plt.scatter(fake[0]/scale, fake[1]/scale, s=60, facecolor='None', edgecolor='magenta', marker='o')
                         
                 
                if(m>=3):
                    
                    print('Found one optimal in wide', field)
                    ninit = ninit + 1
                    
    #print('For iteration', i, 'found', ninit)
    print(i, ninit, file=final)
    nt = nt + ninit
    #print(nt)
                    
         
###-----------------------------------------------------------------------------------------------------------------------
    
fields = ['w1']
scale = 100 #to make things more manageable, scale everything by   
fwhm = 1000/scale # 632 pix or 0.032 degrees represents a box of side 1Mpc in size @ z~1.5 
minp = 1 
maxp = 19354 #For the deep fields
level = 2.7
limgroup = 70/100
its = 100
version = 5
nt = 0

final = open('test_new_5ormore.dat','w')
print('# Number of fake groups found in W1')
print('# 0 Iteration', file=final)
print('# 1 Number of groups found', file=final)

### Define names of patches in w1 and w4
patchesw1 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_w1.txt', dtype=str)
patchesw4 = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_w4.txt', dtype=str)

patcheswide = {'w1':patchesw1, 'w4':patchesw4}

#### Define sigma for each gaussian
sig = (round(fwhm/2.355,0))

#### Make 2D Gaussian
gsmall, xx, yy = makeGaussian2(5*sig, sig)
size = gsmall.shape[0]
print('Made gsmall! It has a shape of ', gsmall.shape)



i = 1

while(i<=its):
    warray = []
    slzerosfull = [] 
    
    ### Make Empty arrays for each field
    slzerosfull = make_empty_arrays(scale)    

    
    for field in fields:
        
        print('Field', field, 'iteration', i)
        cat = np.genfromtxt('/Volumes/Liz/fake_env/output/whole_masks_v'+str(version)+'/'+field+'/test_fakecat_new7_'+field+'_cluster_it'+str(i)+'_v'+str(version)+'.dat')
        warray = []
        
        for peg in cat:
        
            # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
            x = peg[1]
            y = peg[2]
            
             
            xg = round((x/scale),0)+round((gsmall.shape[0]/2),0)
            yg = round((y/scale),0)+round((gsmall.shape[0]/2),0)
            
            
            vec = np.hstack((float(x), float(y)))
            
            if(len(warray)==0):
                warray = vec
            else:
                warray = np.vstack((warray, vec))
                
        
            ### Add this gaussian to fake field
            slzerosfull[field][xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)] = slzerosfull[field][xg-round((size/2),0):xg+round((size/2),0), yg-round((size/2),0):yg+round((size/2),0)]+gsmall
        
        #show_density_plot(slzerosfull, warray, i, field)
        select_dense_environments(slzerosfull, i, field, warray, nt)
        
        
        
    i = i + 1
    
final.close()
print('t_total=',timeit.time.time()-start, 's')
os.system('say "Liz python is done"')   