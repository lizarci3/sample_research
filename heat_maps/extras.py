# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 12:42:04 2017

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt


'''
Module to take into account masked regions when positioning random galaxies
'''

def create_bool_flags(field):
    
    '''
    Create a boolean array that contains the flags of Deep fields, not needed apparently they already exist
    '''
    filters = ['gc']#, 'zc', 'kc']
    total = np.zeros((19354,19354))
    
    for filt in filters:
        stimage = fits.open('/Users/osejo/Desktop/Tesis/gzHK/photometry/wmask/d'+str(field)+'/111106/ascii/stmask'+str(field)+'.'+filt+'.fits', memmap=True)
        stdata = stimage[0].data
        total = total + stdata
        
        
        stimage.close()
    total[total>0] = False
    total[total==0] = True
    return(total)

def test_existing_masks():
    test = np.load('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_d1.npy')

    return(test)
    
def transfer_masks():
    
    '''
    Transfer all masks to my external drive or to mahone
    '''
    print('Starting')
    field = 'w4' 
    
    names = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_'+field+'.txt', dtype=str)
    
    for patch in names:
        
        print(patch)
        #print('scp osejo@mars.smu.ca:/disc/disc61/tmp/osejo/mahone/background/test_'+field+'_'+patch+'_gzkmask.npy /Volumes/Liz/masks_fake_env/')
        os.system('scp osejo@mars.smu.ca:/disc/disc61/tmp/osejo/mahone/background/test_'+field+'_'+patch+'_gzkmask.npy /Volumes/Liz/masks_fake_env/'+field)

def transfer_masks_deep():
    
    '''
    Transfer all masks to my external drive or to mahone
    '''
    print('Starting')
    fields = ['d2', 'd3', 'd4'] 
     
    for f in fields:
        
        print(f)
        #print('scp /Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_'+f+'.npy  larcilao@mahone.ace-net.ca::/home/larcilao/nqs/fake_groups/masks')
        os.system('scp /Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/background/mask_gzk_'+f+'.npy  larcilao@mahone.ace-net.ca:/home/larcilao/nqs/fake_groups/masks')

def transfer_masks_wide():
    
    '''
    Transfer all masks to nqs
    '''
    
    print('Starting')
    fields = ['w1', 'w4'] 
    patcheswide = {}
    
    for field in fields:
        
        patcheswide[field] = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_'+field+'.txt', dtype=str)
        print(field)
        
        for patch in patcheswide[field]:
            print(patch)
            #print('scp /Volumes/Liz/masks_fake_env/'+field+'/test_'+field+'_'+patch+'_gzkmask.npy  larcilao@mahone.ace-net.ca:/home/larcilao/nqs/fake_groups/masks')
            os.system('scp /Volumes/Liz/masks_fake_env/'+field+'/test_'+field+'_'+patch+'_gzkmask.npy  larcilao@mahone.ace-net.ca:/home/larcilao/nqs/fake_groups/masks')


def efingmasks():
    
    '''
    Look at the distribution of the original data, then try to create a 
    fake catalog using the mask and see if it agrees. I chose a reprentative 
    shape in W1 to make it easy
    '''
    
    fullcat = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_dw_cd_nostr_maglim_v_mu_2s01z05k_clean.dat')
    
    #chosenone = '022150-050800'
    #chosenone = '022150-041200'
    #chosenone = '022929-041200'
    totest = ['023319-041200', '022150-050800', '022150-041200', '022929-041200']
    
    for chosenone in totest:
        chosenonei = int(chosenone[1:6]+chosenone[7:])
        minp = 1
        maxp = 19354
        
        masktmp = np.load('/Volumes/Liz/masks_fake_env/w1/test_w1_'+chosenone+'_gzkmask.npy')
        w1mask = masktmp.astype(bool)     
        
        real = fullcat[fullcat[:,0]==chosenonei]
        
        plt.figure()
        plt.xlim(0,maxp)
        plt.ylim(0,maxp)
        plt.title(chosenone+' yx pos')
        plt.scatter(real[:,2], real[:,3], marker='o', color='blue')
        
        for gal in real:
            
            flag = False # As long as we can't be sure they are not in a flagged region repeat procedure
            
            while(flag==False):
                
            
               fakex = int(np.random.randint(minp, maxp))
               fakey = int(np.random.randint(minp, maxp))
               
               flag = ((w1mask[fakey][fakex]==False))
               print('Flag for this position', flag)
               
               if(flag==True):
                   print('This pos is good')
                   plt.scatter(fakex, fakey, marker='s', color='red')
        

#field = 1
#totala = create_bool_flags(field)
#test = test_existing_masks()
#np.savetxt('test.dat', totala)
#transfer_masks()
#transfer_masks_wide()
efingmasks()