# -*- coding: utf-8 -*-
"""
Created on Tue May  3 13:30:06 2016

@author: osejo
"""
from __future__ import print_function
from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
import timeit
import pyregion
from scipy.spatial import cKDTree
import matplotlib.path as mplPath
from ellipse_fit import fit_ellipse
import matplotlib.patches as patchesbi
from extras_groups import plot_groups_deep, color_plots_deep, plot_groups_deep_wsf
import datetime

start=timeit.time.time()
now = datetime.datetime.now()

def compare_kdtree_dense_environments(dist, nb):

    #dist = 0.016 #radius of search
    #nb = 0 #number of neighbours within this radius
    
    xtop, ytop = kdtree_matching(nb, dist)
    plt.scatter(xtop, ytop, marker='h', edgecolor='green', facecolor='None', s=5)

def create_patches(field):

    '''
    In the current figure, plot patches
    '''
    
    region_name = '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/masks/deep/d'+str(field)+'.gzhk_xy_scaled.reg'
    
    r = pyregion.open(region_name)
    
    patch_list, artist_list = r.get_mpl_patches_texts(origin=0)
    
    for p in patch_list:
        #print(p)
            
        p.set_facecolor('gray')
        p.set_edgecolor('gray')
        #ax.add_patch(p)
        plt.gca().add_patch(p)
        
 

    for t in artist_list:
        #ax.add_patch(t)    
        plt.gca().add_patch(t)
        
    if(field==4):
        extra_region = '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/d4_ext_scaled.reg'
        
        re = pyregion.open(extra_region)
    
        patch_liste, artist_liste = re.get_mpl_patches_texts(origin=0)
        
        for pe in patch_liste:
           
                
            pe.set_facecolor('gray')
            pe.set_edgecolor('gray')
            #ax.add_patch(p)
            plt.gca().add_patch(pe)
            
     
    
        for te in artist_liste:
            #ax.add_patch(t)    
            plt.gca().add_patch(te)
            
    if(field==3):
        extra_region = '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/d3_ext_scaled.reg'
        
        re = pyregion.open(extra_region)
    
        patch_liste, artist_liste = re.get_mpl_patches_texts(origin=0)
        
        for pe in patch_liste:
           
                
            pe.set_facecolor('gray')
            pe.set_edgecolor('gray')
            #ax.add_patch(p)
            plt.gca().add_patch(pe)
            
     
    
        for te in artist_liste:
            #ax.add_patch(t)    
            plt.gca().add_patch(te)
            
    if(field==1):
        extra_region = '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/d1_ext_scaled.reg'
        
        re = pyregion.open(extra_region)
    
        patch_liste, artist_liste = re.get_mpl_patches_texts(origin=0)
        
        for pe in patch_liste:
           
                
            pe.set_facecolor('gray')
            pe.set_edgecolor('gray')
            #ax.add_patch(p)
            plt.gca().add_patch(pe)
            
     
    
        for te in artist_liste:
            #ax.add_patch(t)    
            plt.gca().add_patch(te)
            
def kdtree_matching(nb, dist):
    
    '''
    How does this compare with kdtree-closest nb matching?
    I will only return objects who have nb number of neighbours within a distance dist, 
    And this wil be two free parameters that will help me reproduce the output after a certain
    contour overdensity: Say every yellow on the hot cmap
    '''
    array = []
    
    rae = dpl[:,4]
    dece = dpl[:,5]
    
    tree = cKDTree(zip(rae, dece))
    points =zip(rae, dece)   
    indices = tree.query_ball_point(points, r=dist)
    
    for i in xrange(len(indices)):
        p=[]
        
        if (len(indices[i])>nb):
            #print (dpl[i][4], dpl[i][5])
            p = np.hstack((dpl[i][4], dpl[i][5]))
            
            if(len(array)==0):
                array = p
            else:
                array = np.vstack((array, p))
    
    x = (array[:,0]-38.8298006033)/(-0.0052057334996)
    y = (array[:,1]+7.49747268452)/(0.00516703956132)        
    return(x,y)        
    
    #print(array)


def makeGaussian2(size, sigma):
    
    if(size%2 != 0): #make it centered on a whole pixel
        size = size+1
    
    
    #x = np.arange(-size/2, size/2) #centered mesh
    #y = np.arange(-size/2, size/2)
    
    x=np.arange(size)-int(size/2)
    y=np.arange(size)-int(size/2)
    xx, yy = np.meshgrid(x,y)
    
    array = ((1/(2*np.pi*(sigma)**2))*np.exp(-((xx)**2 + (yy)**2)/(2*sigma**2)))
    arraym = (1/np.amax(array))*array
    return arraym, xx, yy


def plot_shifted_full_matrix():
    
    plt.figure(figsize=(10,4), dpi=100)
    plt.xlim(0,1700)
    plt.ylim(0,800)
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.axis('equal')
    plt.title('PEGs W1, gaussian density plot, Ks<20.5 fwhm= '+str(fwhm)+' pix') 
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')     
    plt.imshow(slzeros.T, cmap=plt.get_cmap('hot'), origin='lower')
    plt.colorbar() 

def select_overdense_regions(level, scale, slzeros):
    
    '''
    In this subroutine we select objects that live in overdense Ks < limoffield regions
    BUT once you have the groups up to the limit of the corresponding field you can stop at 
    20.5 if only want a sub sample. Each dp catalog is complete up until the completeness limit
    of each field
    '''
    
    lizcrazy= np.where(slzeros>=level)
    print('Will be selecting objects with values above '+str(level))
    xylizcrazy = np.c_[lizcrazy[0], lizcrazy[1]]

    xylizcrazy[:,0]=xylizcrazy[:,0]-int(gsmall.shape[0]/2)
    xylizcrazy[:,1]=xylizcrazy[:,1]-int(gsmall.shape[0]/2)

    #plt.scatter(xylizcrazy[:,0], xylizcrazy[:,1], marker='d', edgecolor='blue', facecolor='None', s=5)


    for item in dp:
        
        x = item[2]
        y = item[3]
        
        xg = int(x/scale)
        yg = int(y/scale)
    
        for dense in xylizcrazy:
            if(int(xg)==int(dense[0]))and(int(yg)==int(dense[1])):
                #print(field, item[1], xg, yg, item[2], item[3], item[4], item[5], item[7], item[8])#, file=dense_file)
                plt.scatter(item[2]/scale, item[3]/scale, marker='d', edgecolor='Blue', facecolor='None', s=10)
            
def xyradec(field):
    
    '''
    Find a suitable conversion x,y ->ra, dec
    
    '''
    
    data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_d'+str(field)+'.dat')
    
    x = data[:,2]/scale
    y = data[:,3]/scale
    
    ra = data[:,4]
    dec = data[:,5]
    
    mx, bx = np.polyfit(x, ra, 1)
    my, by = np.polyfit(y, dec, 1)
    
    
    return(mx, bx, my, by)
    
def transform_axis(field, mx, bx, my, by):
    
    if(field==1):
        
        tick_lblsx = [36.2, 36.4, 36.6, 36.8]
        #tick_lblsx = [36.1, 36.2, 36.3, 36.4, 36.5, 36.6, 36.7, 36.8]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [-4.9, -4.7, -4.5, -4.3, -4.1]
        #tick_lblsy = [-4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1]
        convy = (tick_lblsy-by)/(my)   
        
    elif(field==2):
        
        tick_lblsx = [149.8, 150, 150.2, 150.4, 150.6]
        #tick_lblsx = [149.8, 149.9, 150, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [1.8, 2.0, 2.2, 2.4, 2.6]
        #tick_lblsy = [1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8]
        convy = (tick_lblsy-by)/(my)
        
    elif(field==3):
        
        tick_lblsx = [214.4, 214.8, 215.2,]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [52.4, 52.6,  52.8, 53]
        convy = (tick_lblsy-by)/(my)
        
    elif(field==4):
        
        tick_lblsx = [333.4, 333.8,  334.2]
        #tick_lblsx = [333.4, 333.5, 333.6, 333.7, 333.8, 333.9, 334, 334.1, 334.2]

        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [-18.2, -18, -17.8, -17.6, -17.4]
        #tick_lblsy = [-18.2, -18.1, -18, -17.9, -17.8, -17.7, -17.6, -17.5, -17.4]
        
        convy = (tick_lblsy-by)/(my)

    tick_lblsx = np.array(tick_lblsx)
    tick_locsx = convx.tolist()
    tick_lblsy = np.array(tick_lblsy)
    tick_locsy = convy.tolist()

    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
    
def select_groups(cont_dense, dp, fwhm, level, limgroup):
    
    x = cont_dense
    
    i = 0

    for i in range(len(x.allsegs[1])): # For all 100% contour paths
    
       
    
        mass = 0
        volume = 0
        area = 0
        
        array = []
        ggroup = []
        
        cont = x.allsegs[1][i] #ith 100% contour level
    
        # Find its nearest contour, but only those that are on x% level (the x is user defined)
        contmin, segmin, indmin, xmin, ymin, dmin = x.find_nearest_contour(cont[0][0], cont[0][1], indices = [0], pixel = False)
        
        bbPath = mplPath.Path(x.allsegs[contmin][segmin])
    
        # Find the volume of this contour by modeling like an ellipse
        axis, phi, center = fit_ellipse(x.allsegs[contmin][segmin])    
    
        a = axis[0]*scale*0.187*8.462/1000 #axis in Mpc  
        b = axis[1]*scale*0.187*8.462/1000
    
        #a = axis[0]*scale*0.128*8.462/1000 #axis in Mpc
        #b = axis[1]*scale*0.128*8.462/1000
        
        c = (a+b)/2
        print(axis[0],axis[1])
        
        ad = axis[0]*scale*0.187/3600 #axis in Deg
        bd = axis[1]*scale*0.187/3600
        
        #ad = axis[0]*scale*0.128/3600 #axis in Deg
        #bd = axis[1]*scale*0.128/3600
        
        print(ad,bd)
        
        volume = 4*np.pi*a*b*c/3 #in Mpc^3
        area = np.pi*a*b # in Mpc^2
        aread = np.pi*ad*bd # in deg^2
        print('Area in deg', aread)
        
        for gal in dp:
            
             vec = []
             xa = gal[2]
             ya = gal[3]
    
             xga = (xa/scale)#.astype(int)
             yga = (ya/scale)#.astype(int)
            
             C = bbPath.contains_point((xga, yga)) # search if the  scaled position is within the x% level

             if(C==True):
                 m = 10**(-0.348*gal[7]+18.284)
                 mass = mass+m
                 
                 #vec = np.hstack((i, gal))
                 
                 if(len(array)==0):
                     array = gal
                 else:
                     array = np.vstack((array, gal))
                     
        # How many of the galaxies in the group is massive? (at least 3 in group to classify)
                     
        ggroup = array[array[:,7]<20.5]
        
        if(len(ggroup)>=3): # If there are at least 3 massive galaxies in group
        
            # ID the approx center of the group              
            medx = np.mean(array[:,2])
            medy = np.mean(array[:,3])
         
            # Find the number density of this group
            nd = len(array[:,0])/volume
        
            # Find the mass density for this group
            md = np.log10(mass)/volume
        
            # Find the number area density
            nad = len(array[:,0])/area
        
        
            print(field, i, medx, medy, np.log10(mass), volume, nd, md, nad, aread, file=charac) 
            plt.figure(num=field)
            plt.scatter(array[:,2]/100, array[:,3]/100, marker='d', edgecolor='Blue', facecolor='None', s=10)
            
    
            #np.savetxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*100))+'/test_new_group_'+str(i)+'_'+str(limgroup*100)+'%_d'+str(field)+'_fwhm_'+str(int(fwhm*100))+'_l'+str(level)+'.dat', array, fmt='%s')


def find_sf(bbPath,i,aread,afield,ngroup):
    
        print('Hola, in find_sf')
        sfarray = []
        sdbri = sd[sd[:,7]<20.5]
        
        backs = len(sd[:,0])/afield # Number of SFGs in this field per deg^2 (up to the completeness of the field)
        #print('Background correction SFGs per deg^2', backs, 'Max mag SFGs', max(sd[:,7]))
        
        backsb = len(sdbri[:,0])/afield # Number of SFGs in this field per deg^2 (up to the completeness of the field)
        #print('Background correction SFGs per deg^2', backs, 'Max mag SFGs', max(sd[:,7]))
        
        
        for sfgal in sd:
            
            S = bbPath.contains_point((sfgal[2]/scale, sfgal[3]/scale))
            
            if(S==True):
                
                if(len(sfarray)==0):
                    sfarray = sfgal
                else:
                    sfarray = np.vstack((sfarray, sfgal))
        
        
        #np.savetxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*100))+'/updated_positions/test_sfin_group_'+str(ngroup)+'_'+str(limgroup*100)+'%_d'+str(field)+'_fwhm_'+str(int(fwhm*100))+'_l'+str(level)+'.dat', sfarray, fmt='%s')
        
        if(len(sfarray)!=0): #if you found anyone in group
            print(sfarray.shape)
            
            if(sfarray.shape[0]==88):
                sfarray = sfarray.reshape(1,88)
                
            nggds = len(sfarray[:,0])/aread #Number of elements in group per deg^2
            
            print('#SF GALS')
            print('Original #', len(sfarray), 'Original #/deg^2', nggds)
            
            nggdcs = nggds - backs # Number of elements per deg^2 corrected for foreground/background contamination in group
            print('Corrected elements in group per deg^2', nggdcs, 'Group', i)
            nggs = nggdcs*aread  # Number of corrected elements in group
            print('Corrected elements in group', nggs, 'Group', i)
        
        ### For Ks < 20.5
        
        brisf = sfarray[sfarray[:,7]<20.5]
        
        if(len(brisf)!=0): #if you found anyone in group
            print('Bright', brisf.shape)
            
            nggdsb = len(brisf[:,0])/aread #Number of elements in group per deg^2
            print('OriginalB #', len(brisf), 'OriginalB #/deg^2', nggdsb)
            
            nggdcsb = nggdsb - backsb # Number of elements per deg^2 corrected for foreground/background contamination in group
            print('CorrectedB elements in group per deg^2', nggdcsb, 'Group', i)
            nggsb = nggdcsb*aread  # Number of corrected elements in group
            print('CorrectedB elements in group', nggsb, 'Group', i)
     
        else:
            
            nggsb = 0 #Did not find any bright SFG in group
            
        return(len(sfarray[:,0]), nggs, len(brisf[:,0]), nggsb)
        
        
        ### If you want to see them
        #print(sfarray.shape)
        #plt.figure(num=field)
        #plt.scatter(sfarray[:,2]/scale, sfarray[:,3]/scale, marker='p', edgecolor='blue', facecolor='None', s=40)
        #plt.scatter(sfarray[:,2]/scale, sfarray[:,3]/scale, marker='p', edgecolor='cyan', facecolor='None', s=0.5*40)
        
        
        
def select_groups_corrected(cont_dense, dp, sd, fwhm, level, limgroup, field, xall, yall, msall, xmass, ymass, msmass, tempx, tempy, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb):
        
    '''
    The difference between this one and the previous function is that
    in this one I will attempt to correct for background numbers
    '''
    areas_deep = {'1':0.692839,'2':0.911096,'3':0.452202,'4':0.459256}
    #print(field, 'Field')
    
    ###To correct for background numbers of PEGs: How many in this field per deg2
    afield = areas_deep[str(field)]
    back = len(dp[:,0])/afield # Number of PEGs in this field per deg^2 (up to the completeness of the field)
    backbright = len(dpl[:,0])/afield # Number of PEGs in this field per deg^2 (up to Ks < 20.5)
    
    
    ###To correct for background numbers of SFGS+PEGs: How many in this field per deg2
    backsfpe = (len(dp[:,0])+len(sd[:,0]))/afield # Number of SFGS+PEGs in this field per deg^2 (up to the completeness of the field)
    #sdbright = sd[sd[:,7]<20.5]    
    backbrightsfpe = (len(dpl[:,0])+len(sd[sd[:,7]<20.5]))/afield # Number of PEGs in this field per deg^2 (up to Ks < 20.5)
    
    ### To correct for background in mass estimates down to Ks < 23: total mass per deg^2 of field
    mfield = np.sum(10**(-0.348*dp[:,7]+18.284))
    backm = mfield/afield  #Total mass PEGs per deg^2 in the field
    #print('M/deg^2 field', backm)
    
    
    ### To correct for background in mass estimates down to Ks<20.5: total mass per deg^2 of field
    mfieldbri = np.sum(10**(-0.348*dpl[:,7]+18.284))
    backmbri = mfieldbri/afield  #Total mass PEGs per deg^2 in the field
    #print('M/deg^2 field', backm)
    
    
    
    
    ## Search for groups inside these contours
    x = cont_dense
    ngroup = 1 #number assigned to complete group 
    
    i = 0

    for i in range(len(x.allsegs[1])): # For all 100% contour paths
        print('Contour', i)
        mass = 0
        volume = 0
        area = 0
        
        array = []
        ggroup = []
        
        cont = x.allsegs[1][i] #ith 100% contour level
    
        # Find its nearest contour, but only those that are on x% level (the x is user defined)
        contmin, segmin, indmin, xmin, ymin, dmin = x.find_nearest_contour(cont[0][0], cont[0][1], indices = [0], pixel = False)
        
        bbPath = mplPath.Path(x.allsegs[contmin][segmin])
    
        # Find the volume of this contour by modeling like an ellipse
        axis, phi, center = fit_ellipse(x.allsegs[contmin][segmin]) 
        
        print('Starting with', axis[0], axis[1])
        
        a = axis[0]*scale*0.187*8.462/1000 #axis in Mpc  
        b = axis[1]*scale*0.187*8.462/1000
        c = (a+b)/2
        print('In Mpc, for testing', a, b, c)
        
        ad = axis[0]*scale*0.187/3600 #axis in Deg
        bd = axis[1]*scale*0.187/3600
        cd = (ad+bd)/2
        
        print(ad, bd, cd)
    
        aread = np.pi*ad*bd # in deg^2
        volume = 4/3*np.pi*a*b*c # volume in cubic megaparsecs
        
        #print('Area in deg', aread)
        
        for gal in dp:
            
             xa = gal[2]
             ya = gal[3]
    
             xga = (xa/scale)
             yga = (ya/scale)
            
             C = bbPath.contains_point((xga, yga)) # search if the  scaled position is within the x% level

             if(C==True):
                 m = 10**(-0.348*gal[7]+18.284)
                 mass = mass+m
                 
                 
                 if(len(array)==0):
                     array = gal
                 else:
                     array = np.vstack((array, gal))
                     
        # How many of the galaxies in the group is massive? (at least 3 in group to classify)
                     
        ggroup = array[array[:,7]<20.5]
        print('Created group', len(ggroup))
        np.savetxt('test.dat', ggroup, fmt='%s')
        
        if(len(ggroup)>=3): # If there are at least 3 massive galaxies in group
        
            print('Group has morethan 3 bright!')
            
            # Find how many SFGs are in the group
            nsf, nsfc, nbsf, nbsfc = find_sf(bbPath,i,aread,afield,ngroup)
            
            # ID the approx center of the group              
            medx = np.mean(array[:,2])
            medy = np.mean(array[:,3])
            
            '''
            ########## Do Background corrections down to the completeness of field
            '''
            # Find corrected number of elements in group
            nggd = len(array[:,0])/aread #Number of elements in group per deg^2
            nggdc = nggd - back # Number of elements per deg^2 corrected in group
            ngg = nggdc*aread  # Number of corrected elements in group
            #print('Corrected elements Ks 23 in group', ngg, 'Group', i)
            
            ## Find corrected number of sf+pegs in group
            #nsfpe = (len(array[:,0])+len(nsf))/aread
            #nsfpec = nsfpe-backsfpe
            #nsfpecf = nsfpec*aread #number of sf+pegs in group corected
            
            # Find corrected number of elements per Mpc^2
            nareacm = nggdc/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find the CORRECTED number for background number of elements per Mpc^3
            nvcm = nareacm*(3/(4*c))
            
            # Find the mass per squared degree for this group
            #print('Mass of group Ks < 23', mass, np.log10(mass))
            md = mass/aread #Total mass of group/area of group in deg^2
            #print('Group (M)/deg2', md)
            
            # Find CORRECTED mass per squared degree for this group
            mdc = md - backm #mgroup/deg^2 - mfield/deg^2 =  mgroup/deg^2 corrected
            
            # Find CORRECTED total mass of group
            mmmgc = mdc*aread 
        
            # Find CORRECTED mass per Mpc^2
            mdcm = mdc/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find CORRECTED mass per cubic megaparsecs
            mvcm = mdcm*(3/(4*c))
            
            
            '''
            ########## Do Background corrections down Ks < 20.5
            '''
            # Find corrected number of elements in group
            nggdbri = len(ggroup[:,0])/aread #Number of Ks < 20.5 elements in group per deg^2
            nggdcbri = nggdbri - backbright # Number of elements per deg^2 corrected in group
            nggbri = nggdcbri*aread  # Number of corrected elements in group
            #print('Corrected elements Ks 20.5 in group', nggbri, 'Group', i)
            
            # Find corrected number of elements per Mpc^2
            nareacmbri = nggdcbri/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find the CORRECTED number for background number of elements per Mpc^3
            nvcmbri = nareacmbri*(3/(4*c))
            
            # Find the mass for Ks < 20.5 galaxies per squared degree for this group
            massbri = np.sum(10**(-0.348*ggroup[:,7]+18.284))
            
            #print('Mass of group Ks< 20.5', massbri, np.log10(massbri))
            mdbri = massbri/aread #Total mass of group/area of group in deg^2
            #print('Group (M)/deg2', md)
            
            # Find CORRECTED mass per squared degree for this group
            mdcbri = mdbri - backmbri #mgroup/deg^2 - mfield/deg^2 =  mgroup/deg^2 corrected
            
            # Find CORRECTED total mass of group
            mmmgcbri = mdcbri*aread 
        
            # Find CORRECTED mass per Mpc^2
            mdcmbri = mdcbri/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find CORRECTED mass per cubic megaparsecs
            mvcmbri = mdcmbri*(3/(4*c))
            
            
            ####FROM HERE
            
            '''
            # Find the number of objects in this group per deg2
            narea = len(array[:,0])/aread
            
            # Find the CORRECTED for background number of objects in this group per deg2 and Mpc^2
            nareac = narea - back #in degrees squared
            nareacm = nareac/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            
            # Find the CORRECTED number for background per Mpc^3
            nvcm = nareacm*(3/(4*c))
            
        
            # Find the mass per squared degree for this group
            md = mass/aread
            
            # Find CORRECTED mass per squared degree for this group
            mdc = md - backm
            mdcm = mdc/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find CORRECTED TOTAL MASS
            mmmgc = mdc*aread  
            
            # Find corrected mass per cubic megaparsecs
            mvcm = mdcm*(3/(4*c))
            '''
            
        
            #print(field, i, medx, medy, np.log10(mmmgc), volume, nggdc, nareacm, nvcm, np.log10(mdc), np.log10(mdcm), np.log10(mvcm), aread, len(array), ngg, nsf, nsfc, array[0][4], array[0][5], file=charac) 
            
            plt.figure(num=field)
            
            fig = plt.figure(num=field)
            ax = fig.add_subplot(111)
            pat = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='#3c0ff0', lw=2) 
            #ax.add_patch(pat)
            
            
            #### Make markings where the groups are
            
            centpx = np.mean(bbPath.vertices[:,0])
            centpy = np.mean(bbPath.vertices[:,1])
            
            print('centpx', centpx, 'centpy', centpy)
            
            if(field!=4):
                cir = patchesbi.Circle((centpx, centpy), 7, facecolor='None', edgecolor='white', lw=2)
            else:
                cir = patchesbi.Circle((75, 113), 7, facecolor='None', edgecolor='white', lw=2)
                
            ax.add_patch(cir)
            
            if(field==2):
                xsh = 7
                ysh = -20
                       
            else:              
                xsh = 7
                ysh = 10
            
            ax.text(centpx+xsh, centpy+ysh, 'D'+str(field)+'_'+str(ngroup), color='#3c0ff0', fontsize=12, fontweight='bold')
            ax.text(centpx+xsh, centpy+ysh, 'D'+str(field)+'_'+str(ngroup), color='white', fontsize=12, fontweight='bold')
            
            tempx, tempy = plot_groups_deep(field, bbPath, xall, yall, msall, xmass, ymass, msmass, i, tempx, tempy, ngroup)
            tempx, tempy = plot_groups_deep_wsf(field, bbPath, xall, yall, msall, xmass, ymass, msmass, i, tempx, tempy, xsfs, ysfs, xsfsb, ysfsb, ngroup)
            
            
            ## Make colorplots of group
            color_plots_deep(field, i, centpx, centpy, scale, bbPath, ggroup, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb, ngroup)
             
            ### If you want to know the positions (in degrees), of the brightest cluster member
            brightest = array[array[:,7]==np.min(array[:,7])]
            #print('For field', field, 'Group', i, 'Brightest gal has Ks', brightest[:,7])
            #print('And position', brightest[:,4][0], brightest[:,5][0])
            
            print(field, 0, ngroup, medx, medy, np.log10(mmmgc), volume, nggdc, nareacm, nvcm, np.log10(mdc), np.log10(mdcm), np.log10(mvcm), aread, len(array[:,0]), ngg, nsf, nsfc, brightest[:,4][0], brightest[:,5][0], len(ggroup[:,0]), nggbri, nvcmbri, np.log10(mmmgcbri), np.log10(mvcmbri), nbsf, nbsfc, file=charac) 
            ####    0    1    2      3     4        
            
            #plt.scatter(array[:,2]/100, array[:,3]/100, marker='d', edgecolor='Blue', facecolor='None', s=10)
            #plt.text(medx/scale+(0.05*medx/scale), medy/scale+(0.03*medy/scale), 'Log[M*/Msol]='+str(round(np.log10(mass),2)), color='blue', fontsize=10, bbox=dict(facecolor='red', alpha=0.3))
            #plt.text(medx/scale+(0.08*medx/scale), medy/scale-(0.03*medy/scale), r'Log[M*/Msol]/Mpc$^3$='+str(round(np.log10(mvcm),2)), color='blue', fontsize=10, bbox=dict(facecolor='red', alpha=0.3))
            #plt.text(medx/scale-(0.05*medx/scale), medy/scale-(0.03*medy/scale), r'#/Mpc$^3$='+str(round(nvcm,2)), color='blue', fontsize=10, bbox=dict(facecolor='red', alpha=0.3))
          
    
            #np.savetxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*100))+'/updated_positions/test_group_'+str(ngroup)+'_'+str(limgroup*100)+'%_d'+str(field)+'_fwhm_'+str(int(fwhm*100))+'_l'+str(level)+'.dat', array, fmt='%s')
            ngroup = ngroup + 1
            
        else:
            tempx, tempy, centpx, centpy,  bbPath, ggroup, array = 0, 0, 0, 0, 0, 0, 0
        #return(tempx, tempy, bbPath)
        
        
    return(tempx, tempy, field, i, centpx, centpy, scale, bbPath, ggroup, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb)
    

### REAL Data  

print('Starting!')  


fields = [2]
scale = 100 #to make things more manageable, scale everything by   
fwhm = 1000/scale # 632 pix or 0.032 degrees represents a box of side 1Mpc in size @ z~1.5 
tempx = 0
tempy = 0

show_groups = 'no'

print('Do you want to have points of ALL PEGs in your plots? yes/no:')
points = raw_input()

ksmassive = 20.5
print('Do you want to have points of all massive '+str(ksmassive)+' PEGS in your plots? yes/no:')
massive = raw_input() 

print('Do you want to have points of ultra-massive  PEGS in your plots? yes/no:')
ultramassive = raw_input()

if(ultramassive=='yes'):
    print('Ultra-Massive galaxies brighter than Ks?:')
    ultraksmassive = float(raw_input())

#if(massive=='yes'):
    #print('Massive galaxies brighter than Ks?:')
    #ksmassive = float(raw_input())

#if(massive!='yes'):
    #ksmassive = 20.5#19.8 # the idea with the stars is to show which objects make up the gaussian density field and distinguish them from the rest (represented with circles)
    
print('Do you want to add SF Galaxies to density plots? yes/no')
sfin = raw_input()
if (sfin=='yes'):
    print('Magnitude limit for the SF Population?')
    limsf = float(raw_input())
    
print('Do you want to select overdense regions? yes/no:')
overdense=raw_input()

if(overdense=='yes'):
    print('Overdense above?: (e.g. 2.7)')
    level = float(raw_input())
    
    print('Groups will be completed to what percentage of level? (e.g. 70)')
    limgroup = float(raw_input())/100
    
    print('Do you want to see complete groups? yes/no:')
    show_groups = raw_input()


    charac = open('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/smf/updated_positions/test_characteristics_up_'+str(now.year)+str(now.month)+str(now.day)+'_l'+str(level)+'_f'+str(fwhm*100)+'_c_'+str(limgroup)+'_deep.dat','w')
    print('# Groups in the Deep fields selected with a fwhm of '+str(fwhm*100)+' and  a level of '+str(level), file = charac)
    print('# (0) Field', file = charac)
    print('# (1) Place-holder to match Wide files-ignore-', file = charac)
    print('# (2) Group', file = charac)
    print('# (3) X Center of group', file = charac)
    print('# (4) Y Center of group', file = charac)
    print('# (5) Corrected Total Mass [Log(M*/Msol)]', file=charac)
    print('# (6) Volume [Mpc^3] # Modeled using an ellipse', file = charac)
    print('# (7) Number Density [#/deg^2]', file = charac)
    print('# (8) Number Density [#/Mpc^2]', file = charac)
    print('# (9) Number Density [#/Mpc^3]', file = charac)
    print('# (10) Mass Area Density [Log(M*/Msol/deg^2)]', file = charac)
    print('# (11) Mass Area Density [Log(M*/Msol/Mpc^2)]', file = charac)
    print('# (12) Mass Density [Log(M*/Msol/Mpc^3)]', file = charac)
    print('# (13) Area of group [deg^2]', file = charac)
    print('# (14) Number of PEGs in group -NOT corrected for foreground/background contamination-', file = charac)
    print('# (15) Number of PEGs -corrected for foreground/background contaminatio- in group', file = charac)
    print('# (16) Number of SFGs in group', file = charac)
    print('# (17) Number of SFGs -corrected for foreground/background contaminatio- in group', file = charac)
    print('# (18) RA [deg] for brightest item in group', file = charac)
    print('# (19) DEC [deg] for brightest item in group', file = charac)
    print('# (20) Number of PEGs with Ks< 20.5 in group', file = charac)
    print('# (21) Number of PEGs with Ks<20.5 -corrected for foreground/background contaminatio- in group', file = charac)
    print('# (22) Number Density of PEGs with Ks<20.5 -corrected b/f- [#/Mpc^3]', file = charac)
    print('# (23) Corrected Total Mass PEGs Ks<20.5 [Log(M*/Msol)]', file=charac)
    print('# (24) Mass Density PEGs Ks<20.5 -corrected b/f- [Log(M*/Msol/Mpc^3)]', file = charac)
    print('# (25) Number of SFGs down to Ks<20.5 in group', file = charac)
    print('# (26) Number of SFGs Ks<20.5 -corrected for foreground/background contaminatio- in group', file = charac)
    
   
        
for field in fields:
    
    if(field==4):
        dp = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_d4_cd_clean.dat') #artifacts removed
        
    else:
        dp = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_d'+str(field)+'_cd.dat') # These catalogs have been stopped at the completeness of the sample. Double cores and artifacts removed
    
    sd = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/sf_d'+str(field)+'_cd.dat')
   
    sd = sd[sd[:,7]<23.5]
    sdlim = sd[sd[:,7]<20.5]
    
    
    
    dpl = dp[dp[:,7]<20.5]# the contour plots will always stop at the completeness of the wide survey      
    dpm = dp[dp[:,7]<ksmassive]
    
    if(ultramassive=='yes'):
        dpum = dp[dp[:,7]<ultraksmassive]
        mum = 10**((-0.348)*dpum[:,7]+18.284)
        msum = 100*((mum)/max(mum))    
    
    ## The completeness of the SF population in the Deep fields is OK due to gzHKs unlike the Wide fields 
    

    ma = 10**((-0.348)*dp[:,7]+18.284)
    ml = 10**((-0.348)*dpl[:,7]+18.284)
    mm = 10**((-0.348)*dpm[:,7]+18.284)
    
    
    
    msa = 100*((ma)/max(ma))
    msl = 100*((ml)/max(ml))
    msm = 100*((mm)/max(mm))
    
    
    xytoradec = open('test'+str(field)+'.dat','w')
    dense_file = open('test_dense_deep'+str(field)+'.dat','w')

    ra = dpl[:,4]
    dec = dpl[:,5]

    pixx = dpl[:,2]
    pixy = dpl[:,3]


    ###### The extra pixels secure a 5sigma extra padding for edge effects
        
    slzeros = np.zeros((30000/scale, 30000/scale)) #d1


    print(slzeros.shape)

     
    sig = (round(fwhm/2.355,0)) #sigma = 2.355*fwhm
    print('I am working with a sigma of ', sig)


    gsmall, xx, yy = makeGaussian2(5*sig, sig)
    size = gsmall.shape[0]
    print('Made gsmall! It has a shape of ', gsmall.shape)


    plt.figure(num=field, figsize=(10,10), dpi=100) #d1

    ax=plt.subplot(111)
    plt.axis('scaled')
    plt.xlim(max(pixx)/scale, min(pixx)/scale)
    plt.ylim(min(pixy)/scale, max(pixy)/scale)

    
    #plt.title('Gaussian density plot, D'+str(field), fontsize=25) 
    #plt.xlabel('x [pix]')
    #plt.ylabel('y [pix]')
    plt.xlabel('RA [deg]', fontsize=22)
    plt.ylabel('DEC [deg]', fontsize=22)


    for gal in dpl: #construct the gaussian density plot which will always go to Ks < 20.5
    
        x = gal[2]
        y = gal[3]
    
        # shift every point by at least 1/2 gaussian to avoid edge effects
    
        #xg = int(x/scale)+int(gsmall.shape[0]/2) # if you WANT to round up x, y positions
        #yg = int(y/scale)+int(gsmall.shape[0]/2)
        
        xg = x/scale+(round(size/2,0)) #if you don't want to round up x, y positions
        yg = y/scale+(round(size/2,0))
   
        # Write  the xy coordinates and corresponding RA, DEC of the object to make the unit conversion from pix to deg    
        #print(field, gal[1], x, y, gal[4], gal[5], file=xytoradec)    
    
        slzeros[xg-round(size/2,0):xg+round(size/2,0), yg-round(size/2,0):yg+round(size/2,0)] = slzeros[xg-round(size/2,0):xg+round(size/2,0), yg-round(size/2,0):yg+round(size/2,0)]+gsmall
        
    ########## In pixels
    #draw points where galaxies up to a certain Ks reside in the density plot
    
    ## Draw masked regions
    #create_patches(field)
    
  
    ### Draw points of all PEGs down to the actual completeness of your field
    xa = dp[:,2]
    ya = dp[:,3]
    
    xga = (xa/scale).astype(int)
    yga = (ya/scale).astype(int)
    
    
    
    if(points=='yes'):
        
        #plt.scatter(xga, yga, marker='o', edgecolor='darkred', facecolor='None', s=msa) # Plot rounded positions
        #plt.scatter(xga, yga, marker='o', edgecolor='gold', facecolor='None', s=0.5*msa) # Plot rounded positions
        
        plt.scatter(xa/scale, ya/scale, marker='o', edgecolor='darkred', facecolor='None', s=msa) # Plot real positions
        plt.scatter(xa/scale, ya/scale, marker='o', edgecolor='gold', facecolor='None', s=0.5*msa) # Plot real positions
        
        
    #### Plot the most massive objects
        
    xm = dpm[:,2]
    ym = dpm[:,3]
    
    if(ultramassive=='yes'):
        xum = dpum[:,2]
        yum = dpum[:,3]
    
    #xgm = (xm/scale).astype(int)
    #ygm = (ym/scale).astype(int)
    
    
    if(massive=='yes')and(ksmassive==20.5):
        
        #plt.scatter(xgm, ygm, marker='*', edgecolor='darkred', facecolor='None', s=msm) # Plot rounded positions
        #plt.scatter(xgm, ygm, marker='*', edgecolor='gold', facecolor='None', s=0.5*msm) # Plot rounded positions
        
        plt.scatter((xm/scale), (ym/scale), marker='*', edgecolor='darkred', facecolor='None', s=msm)
        plt.scatter((xm/scale), (ym/scale), marker='*', edgecolor='gold', facecolor='None', s=0.5*msm)
   
    if(ultramassive=='yes'):
        
      
        plt.scatter((xum/scale), (yum/scale), marker='*', edgecolor='magenta', facecolor='None', s=2*msum, linewidth=2)
        plt.scatter((xum/scale), (yum/scale), marker='*', edgecolor='None', facecolor='white', s=2*0.5*msum)
        
  
    if(sfin=='yes'):
        sd = sd[sd[:,7]<limsf]
        xsf = sd[:,2]
        ysf = sd[:,3]
        
        msf = 40*((sd[:,7])/max(sd[:,7])) #scale in terms of their Ks magnitudes, since we do not have the masses
        
        marsf = 40*((sd[:,7])/max(sd[:,7]))
        marsfb = 40*((sdlim[:,7])/max(sdlim[:,7]))
    
        plt.scatter(xsf/scale, ysf/scale, marker='*', edgecolor='blue', facecolor='None', s=msf)
        plt.scatter(xsf/scale, ysf/scale, marker='*', edgecolor='cyan', facecolor='None', s=0.5*msf)
        
        
    ### Plot also masked regions
    #create_patches(field)
    
    ### Change axis from x,y to RA,DEC
    
    mx, bx, my, by = xyradec(field) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis(field, mx, bx, my, by)    
    plt.xticks(tick_locsx, tick_lblsx, fontsize=20)
    plt.yticks(tick_locsy, tick_lblsy, fontsize=20)
    
    #### Show Gaussian Density Plot
    plt.imshow(slzeros[(round(gsmall.shape[0]/2,0)):(slzeros.shape[0]-(round(gsmall.shape[0]/2,0))), (round(gsmall.shape[0]/2,0)):(slzeros.shape[0]-(round(gsmall.shape[0]/2,0)))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=20) 
    #plt.contour(slzeros[int(gsmall.shape[0]/2):(slzeros.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzeros.shape[0]-int(gsmall.shape[0]/2))].T, cmap=plt.get_cmap('hot'))
    plt.tight_layout()
    #np.savetxt('gaussian_density_d'+str(field)+'.dat', slzeros[int(gsmall.shape[0]/2):(slzeros.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzeros.shape[0]-int(gsmall.shape[0]/2))].T, fmt='%s')
   
    #### Select overdense regions
 
    if(overdense=='yes'):
        
        #cont_dense = plt.contour(slzeros[int(gsmall.shape[0]/2):(slzeros.shape[0]-int(gsmall.shape[0]/2)), int(gsmall.shape[0]/2):(slzeros.shape[0]-int(gsmall.shape[0]/2))].T, levels=[limgroup*level, level], cmap=plt.get_cmap('Blues'), alpha=0.01)
        cont_dense = plt.contour(slzeros[(round(gsmall.shape[0]/2,0)):(slzeros.shape[0]-(round(gsmall.shape[0]/2,0))), (round(gsmall.shape[0]/2,0)):(slzeros.shape[0]-(round(gsmall.shape[0]/2,0)))].T, levels=[limgroup*level, level], cmap=plt.get_cmap('Blues'), alpha=0.7)        
        plt.clabel(cont_dense, fontsize=50)
        
        
    if(show_groups=='yes')and(field!=3): 
        
        if (sfin!='yes'):
            print('No SF info')
            xsf = 0
            ysf = 0
            msf = 0 
            marsf = 0
            marsfb = 0
            
        # If you want to save the contour path of your groups, use next line
        #select_groups(cont_dense, dp, fwhm, level, limgroup)
        #tempx, tempy, bbPath = select_groups_corrected(cont_dense, dp, fwhm, level, limgroup, field, xa/scale, ya/scale, msa, xm/scale, ym/scale, msm, tempx, tempy, xsf/scale, ysf/scale, sdlim[:,2]/scale, sdlim[:,3]/scale)
        
        ### To make faster color groups
        tempx, tempy, field, i, centpx, centpy, scale, bbPath, ggroup, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb = select_groups_corrected(cont_dense, dp, sd, fwhm, level, limgroup, field, xa/scale, ya/scale, msa, xm/scale, ym/scale, msm, tempx, tempy, xsf/scale, ysf/scale, sdlim[:,2]/scale, sdlim[:,3]/scale, marsf, marsfb)
        # In the image show the information for the FULL (i.e. complete) group
        #plot_info_group(level, scale, fwhm, field, datacharac)
    

    
    
    
    
    
    #### To compare with KD-Tree output
    #compare_kdtree_dense_environments(0.016,0)
    
    #### Save current figure      
    plt.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*scale))+'/d'+str(field)+'_f'+str(int(fwhm*scale))+'_l'+str(level)+'_superplot_test2.pdf', dpi=100)
    
    ## If you want to see the SHIFTED full matrix
    #plot_shifted_full_matrix()
    
    
    xytoradec.close()
    dense_file.close()
    
charac.close()  

print(tempx, tempy) 
    
print('t_total=',timeit.time.time()-start, 's')
os.system('say "Liz python is done"')

       
