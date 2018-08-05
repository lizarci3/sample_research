# -*- coding: utf-8 -*-
"""
Created on Tue May  3 13:30:06 2016

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.image as mpimg
import pylab
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import timeit
import pyregion
from scipy.spatial import cKDTree
import matplotlib.path as mplPath
from ellipse_fit import fit_ellipse
from matplotlib.path import Path
import matplotlib.patches as patchesbi
from extras_groups import plot_groups_wide, plot_groups_wide_wsf, color_plots_wide


start=timeit.time.time()

def compare_kdtree_dense_environments(dist, nb):

    #dist = 0.016 #radius of search
    #nb = 0 #number of neighbours within this radius
    
    xtop, ytop = kdtree_matching(nb, dist)
    plt.scatter(xtop, ytop, marker='h', edgecolor='green', facecolor='None', s=5)

def create_patches():
    
    '''
    In the current figure, plot patches
    '''
    
    patches = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/names_w'+str(fieldw)+'.txt', dtype=str)
    
    for patch in patches:
        region_name = '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/masks/w'+str(fieldw)+'/w'+str(fieldw)+'_'+patch+'_gzk_pix.reg'
    
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

def kdtree_matching(nb, dist):
    
    '''
    How does this compare with kdtree-closest nb matching?
    I will only return objects who have nb number of neighbours within a distance dist, 
    And this wil be two free parameters that will help me reproduce the output after a certain
    contour overdensity: Say every yellow on the hot cmap
    '''
    array = []
    
    rae = pw[:,4]
    dece = pw[:,5]
    
    tree = cKDTree(zip(rae, dece))
    points =zip(rae, dece)   
    indices = tree.query_ball_point(points, r=dist)
    
    for i in xrange(len(indices)):
        p=[]
        
        if (len(indices[i])>nb):
            #print (pw[i][4], pw[i][5])
            p = np.hstack((pw[i][4], pw[i][5]))
            
            if(len(array)==0):
                array = p
            else:
                array = np.vstack((array, p))
    
    x = (array[:,0]-38.8298006033)/(-0.0052057334996)
    y = (array[:,1]+7.49747268452)/(0.00516703956132)        
    return(x,y)        
    
    #print(array)

def makeGaussian(size, fwhm, sigma, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
        
    #return (np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)) #approximation using exponenial functions
    return ((1/(2*np.pi*sigma**2))*np.exp(-((xx)**2 + (yy)**2)/(2*sigma**2))) # symmetric 2D Gaussian distribution
    

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

def shift(fieldw, field):
    
    if(fieldw==1):
        xpatch = {'20241041200':8, '20241050800':8,'20241060400':8, '20631041200':7, '20631050800':7, '20631060400':7, '20631070000':7, '21021041200':6, '21021050800':6, '21021060400':6, '21021070000':6, '21410041200':5, '21410050800':5, '21410060400':5, '21800041200':4, '21800050800':4, '21800060400':4, '22150041200':3, '22150050800':3, '22150060400':3, '22539050800':2, '22539060400':2, '22929041200':1, '22929050800':1, '22929060400':1, '23319041200':0, '23319050800':0, '23319060400':0}
        ypatch = {'20241041200':3, '20241050800':2,'20241060400':1, '20631041200':3, '20631050800':2, '20631060400':1, '20631070000':0, '21021041200':3, '21021050800':2, '21021060400':1, '21021070000':0, '21410041200':3, '21410050800':2, '21410060400':1, '21800041200':3, '21800050800':2, '21800060400':1, '22150041200':3, '22150050800':2, '22150060400':1, '22539050800':2, '22539060400':1, '22929041200':3, '22929050800':2, '22929060400':1, '23319041200':3, '23319050800':2, '23319060400':1}

        shiftx = xpatch[field]
        shifty = ypatch[field]
    
    elif(fieldw==4):
        xpatch = {'220154011900':5, '220154021500':5,'220154031100':5, '220542011900':4, '220542021500':4, '220542031100':4, '220930003100':3, '220930002300':3, '220930011900':3, '220930021500':3, '221318003100':2, '221318002300':2, '221318011900':2, '221318021500':2, '221706002300':1, '221706011900':1, '221706021500':1, '222054002300':0, '222054011900':0}
        ypatch = {'220154011900':2, '220154021500':3,'220154031100':4, '220542011900':2, '220542021500':3, '220542031100':4, '220930003100':0, '220930002300':1, '220930011900':2, '220930021500':3, '221318003100':0, '221318002300':1, '221318011900':2, '221318021500':3, '221706002300':1, '221706011900':2, '221706021500':3, '222054002300':1, '222054011900':2}
    
        shiftx = xpatch[field]
        shifty = ypatch[field]
        
        
    return(shiftx, shifty)


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
    No longer used, now we use the contour regions
    and the peak and x% completeness of that contour to 
    ID  a group
    '''
    
    lizcrazy= np.where(slzeros>=level) #top or level
    print('Will be selecting objects with values above '+str(level))
    xylizcrazy = np.c_[lizcrazy[0], lizcrazy[1]]

    xylizcrazy[:,0]=xylizcrazy[:,0]-int(gsmall.shape[0]/2)
    xylizcrazy[:,1]=xylizcrazy[:,1]-int(gsmall.shape[0]/2)

    #plt.scatter(xylizcrazy[:,0], xylizcrazy[:,1], marker='d', edgecolor='blue', facecolor='None', s=5)


    for item in pw:
    
        field = str(int(item[43]))
        shiftx, shifty = shift(fieldw, field)
    
        x = item[2]+shiftx*(19354-968)
        y = item[3]+shifty*(19354-1290)
    
        xg = int(x/scale)
        yg = int(y/scale)
    
        for dense in xylizcrazy:
            if(int(xg)==int(dense[0]))and(int(yg)==int(dense[1]))and(item[17]<select):
                print(item[43], item[1], xg, yg, item[4], item[5], item[17], item[23], file=dense_file)
                #plt.scatter(xg, yg, marker='d', edgecolor='Blue', facecolor='None', s=10)

def xyradec(fieldw):
    
    '''
    Find a suitable conversion x,y ->ra, dec
    
    '''
    
    data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_w'+str(fieldw)+'.dat')
    
    x = data[:,2]
    y = data[:,3]
    
    ra = data[:,4]
    dec = data[:,5]
    
    mx, bx = np.polyfit(x, ra, 1)
    my, by = np.polyfit(y, dec, 1)
    
    return(mx, bx, my, by)
 
def make_empty_array(fieldw):
    
    if(fieldw==1):
        #slzeros = np.zeros((174000/scale, 70000/scale)) # shape of original mesh, with some padding for W1
        slzeros = np.zeros((300000/scale, 90000/scale))
    elif(fieldw==4):
        slzeros = np.zeros((120000/scale, 90000/scale)) # shape of original mesh, with some padding for W4

    return(slzeros)


def transform_axis(fieldw, mx, bx, my, by):
    
    if(fieldw==1):
        
        tick_lblsx = [31,32,33,34,35,36,37,38]
        #tick_lblsx = np.arange(31,39,0.02) #for cut-outs
        
        tick_lblsy = [-6.5, -6,-5.5, -5, -4.5]
        #tick_lblsy = np.arange(-6.5,-4.4,0.01) # for cut-outs
        
    elif(fieldw==4):
        
        tick_lblsx = [335,334,333,332,331]
        ##tick_lblsx = np.arange(331,336,1)
        #tick_lblsx = np.arange(331,337,0.01)
        
        tick_lblsy = [-0.5, 0, 0.5, 1, 1.5, 2.0, 2.5]
        #tick_lblsy = np.arange(-0.5,2.6,0.01)
        


    convx = (tick_lblsx-bx)/(mx)
    convy = (tick_lblsy-by)/(my)
    
    tick_lblsx = np.array(tick_lblsx)
    tick_locsx = convx.tolist()
    tick_lblsy = np.array(tick_lblsy)
    tick_locsy = convy.tolist()

    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
    

def select_groupsw(cont_densew, pw, fwhm, level, limgroupw, characw): 
    
    x = cont_densew # Comes from the contour plots

    # How many 100% level contours do we have? (most likely this will represent the number of groups you have in that field)
    # some may be fake check at the end

    i = 0

    for i in range(len(x.allsegs[1])): # for all 100% level contours find its closest x% (user defined) contour that corresponds to that 100%
        
        mass = 0
        volume = 0
        area = 0
        
        array = []
        ggroup = []
    
        cont = x.allsegs[1][i] #ith 100% contour level
        
        # Find its nearest x% neighbouring contour
        contmin, segmin, indmin, xmin, ymin, dmin = x.find_nearest_contour(cont[0][0], cont[0][1], indices = [0], pixel = False)
        
        bbPath = mplPath.Path(x.allsegs[contmin][segmin])  
        
        # Find the volume of this contour by modeling like an ellipse
        axis, phi, center = fit_ellipse(x.allsegs[contmin][segmin])    
    
        a = axis[0]*100*0.187*8.462/1000  
        b = axis[1]*100*0.187*8.462/1000
        c = (a+b)/2
        
        ad = axis[0]*100*0.187/3600 #axis in Deg
        bd = axis[1]*100*0.187/3600
    
        volume = 4*np.pi*a*b*c/3
        area = np.pi*a*b
        aread = np.pi*ad*bd # in deg^2
        
        for gal in pw:
            
            field = str(int(gal[43]))
            shiftx, shifty = shift(fieldw, field)
    
            # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
            xa = gal[2]+shiftx*(19354-968)
            ya = gal[3]+shifty*(19354-1290)
    
            
            xag = (xa/scale) #float or int
            yag = (ya/scale)
            
             
            C = bbPath.contains_point((xag, yag)) # search if the  scaled position is within the x% level

            if(C==True):
                m = 10**(-0.348*gal[17]+18.284)
                mass = mass+m
                if(len(array)==0):
                    array = gal
                else:
                    array = np.vstack((array, gal))
                     
        # How many of the galaxies in the group is massive? (at least 3 in group to classify)
                     
        ggroup = array[array[:,17]<20.5]
        
        if(len(ggroup)>=3): # If there are at least 3 massive galaxies in group
        
            # ID the approx center of the group              
            medx = np.mean(array[:,2])
            medy = np.mean(array[:,3])
            #fieldg = np.unique(array[:,43])#np.mean(array[:,43]) #or numpy unique
         
            # Find the number density of this group
            nd = len(array[:,0])/volume
        
            # Find the mass density for this group
            md = np.log10(mass)/volume
        
            # Find the number area density
            nad = len(array[:,0])/area
        
        
            print(fieldw, i, medx, medy, np.log10(mass), volume, nd, md, nad, aread, file=characw)
            
            plt.figure(num=fieldw*10)
            
            
            for element in array:
                
                fieldarray = str(int(element[43]))
                shiftx, shifty = shift(fieldw, fieldarray)
    
                # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
                xarray = (element[2]+shiftx*(19354-968))/scale
                yarray = (element[3]+shifty*(19354-1290))/scale
    
                #plt.scatter(xarray, yarray, marker='d', edgecolor='Blue', facecolor='None', s=10)
            
          
    
            #np.savetxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*100))+'/test_group_'+str(i)+'_'+str(limgroupw*100)+'%_w'+str(fieldw)+'_fwhm_'+str(int(fwhm*100))+'_l'+str(level)+'.dat', array, fmt='%s')

def find_sf(bbPath,i,aread, afield,ngroup):
    
        #print('Hola')
        sfarray = []
        
        ###To correct for background numbers: How many in this field per deg2
       
        backs = len(sw[:,0])/afield # Number of SFGs in this field per deg^2 (up to the completeness of the field)
        #print('Background correction SFGs per deg^2', backs, 'Max mag SFGs', max(sw[:,17]))
        
            
        for sfgal in sw: #for each sf in catalog with Ks < 20.5
             
            overlapx = 968
            overlapy = 1290
             
            fieldsf = str(int(sfgal[43]))
             
            if(fieldw==4)and((fieldsf==220930003100)or(fieldsf==221318002300)):
                overlapy=1741
             
            shiftx, shifty = shift(fieldw, fieldsf)
            
            xasf = sfgal[2]+shiftx*(19354-overlapx)
            yasf = sfgal[3]+shifty*(19354-overlapy)
    
            xgasf = (xasf/scale)
            ygasf = (yasf/scale)
            
              
            S = bbPath.contains_point((xgasf, ygasf))
            
            if(S==True): #find if it is inside boundaries of group
                
                if(len(sfarray)==0):
                    sfarray = sfgal
                else:
                    sfarray = np.vstack((sfarray, sfgal))
        
        # Find corrected number of elements in group
        
        sfadata = np.array(sfarray)
        #print(i, sfadata.shape, sfadata.shape[0], len(sfarray), len(sfarray))
        
        ## I removed this one for testing, put back if you want catalog
        #np.savetxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*100))+'/updated_positions/test_sfin_group_'+str(ngroup)+'_'+str(limgroupw*100)+'%_w'+str(fieldw)+'_fwhm_'+str(int(fwhm*100))+'_l'+str(level)+'.dat', sfarray, fmt='%s')
        
        nggs = 0
        
        if(sfadata.shape[0]==44):
            sfadata = sfadata.reshape(1,44)
       
        if(len(sfadata)!=0): #if you found anyone in group
            #print(sfadata.shape)
            
            nggds = len(sfadata[:,0])/aread #Number of elements in group per deg^2
            #print('Original #', len(sfadata), 'Original #/deg^2', nggds)
            
            nggdcs = nggds - backs # Number of elements per deg^2 corrected for foreground/background contamination in group
            #print('Corrected elements in group per deg^2', nggdcs, 'Group', i)
            nggs = nggdcs*aread  # Number of corrected elements in group
            #print('Corrected elements in group', nggs, 'Group', i)
            
        
        return(len(sfadata), nggs)
    
        
        
         
        ### IF you want to plot them
        '''
        plt.figure(num=fieldw*10)
        overlapx = 968
        overlapy = 1290
        
        
        if(sfadata.shape[0]!=0)and(sfadata.shape[0]!=44): #IF they are empty or have just one element, the following script won't work
            
            for itemsf in sfadata:
                        
                fieldasf = str(int(itemsf[43]))
                
                if(fieldw==4)and((fieldasf==220930003100)or(fieldasf==221318002300)):
                    overlapy=1741
                        
                shiftx, shifty = shift(fieldw, fieldasf)
                    
                # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
                sfx = itemsf[2]+shiftx*(19354-overlapx)
                sfy = itemsf[3]+shifty*(19354-overlapy)
            
                        
                plt.scatter(sfx/scale, sfy/scale, marker='p', edgecolor='blue', facecolor='None', s=40)
                plt.scatter(sfx/scale, sfy/scale, marker='p', edgecolor='cyan', facecolor='None', s=0.5*40)
        
        elif(sfadata.shape[0]==44): #when there is only 1 element
            
            fieldasf = str(int(sfadata[43]))
            print(fieldasf)
                
            if(fieldw==4)and((fieldasf==220930003100)or(fieldasf==221318002300)):
                overlapy=1741
                    
            shiftx, shifty = shift(fieldw, fieldasf)
            
                
            # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
            sfx = sfadata[2]+shiftx*(19354-overlapx)
            sfy = sfadata[3]+shifty*(19354-overlapy)
        
                    
            plt.scatter(sfx/scale, sfy/scale, marker='p', edgecolor='blue', facecolor='None', s=40)
            plt.scatter(sfx/scale, sfy/scale, marker='p', edgecolor='cyan', facecolor='None', s=0.5*40)
        '''
        
            
def select_groups_correctedw(cont_densew, pw, fwhm, level, limgroupw, fieldw, scale, sw, gsmall):
    
    '''
    The difference between this one and the previous function is that
    in this one I will attempt to correct for background numbers
    '''
    areas_wide = {'1':15.53,'4':9.56}
    #print('Mag lim of selected objects', max(pw[:,17]))
    
    ###To correct for background numbers: How many in this field per deg2
    afield = areas_wide[str(fieldw)]
    back = len(pw[:,0])/afield # Number of PEGs in this field per deg^2 (up to the completeness of the field)
    #print('Background correction per deg^2', back)
    
    ### To correct for background in mass estimates: total mass per deg^2 of field
    mfield = np.sum(10**(-0.348*pw[:,17]+18.284))
    backm = mfield/afield
    
    
    
    ## Search for groups inside these contours
    x = cont_densew
    ngroup = 1 #number assigned to complete group 
    i = 0

    for i in range(len(x.allsegs[1])): # For all 100% contour paths
        
        #print('I am in contour', i)
        
        massg = 0
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
        c = (a+b)/2
        
        
        ad = axis[0]*scale*0.187/3600 #axis in Deg
        bd = axis[1]*scale*0.187/3600
        cd = (ad+bd)/2
        
    
        aread = np.pi*ad*bd # in deg^2
        volume = 4/3*np.pi*a*b*c # volume in cubic megaparsecs
        
        #print('Area in deg', aread)
        
        for gal in pw:
                
             overlapx = 968
             overlapy = 1290
             
             field = str(int(gal[43]))
             
             if(fieldw==4)and((field==220930003100)or(field==221318002300)):
                 overlapy=1741
             
             shiftx, shifty = shift(fieldw, field)
            
             xa = gal[2]+shiftx*(19354-overlapx)
             ya = gal[3]+shifty*(19354-overlapy)
    
             xga = (xa/scale)
             yga = (ya/scale)
            
             C = bbPath.contains_point((xga, yga)) # search if the  scaled position is within the x% level

             if(C==True):
                 m = 10**(-0.348*gal[17]+18.284)
                 massg = massg+m
                 
                 
                 if(len(array)==0):
                     array = gal
                 else:
                     array = np.vstack((array, gal))
                     
        # How many of the galaxies in the group is massive? (at least 3 in 100% equi-level path to classify)
                     
        ggroup = array[array[:,17]<20.5]
        
        if(len(ggroup)>=3): # If there are at least 3 massive galaxies in group
            
            
            #print('Contour', i, 'Did make the cut') 
            #print(bbPath, file=test_contours2)
            #from catalog_groups import combinedcat
            #combinedcat(bbPath, i, contmin, segmin, indmin, xmin, ymin, dmin, fieldw, scale, mastercat)
            
            # Find how many SFGs are in the group
            nsf, ncsf = find_sf(bbPath,i,aread,afield,ngroup)  #it returns the original and corrected number of SFGs in ith group          
            
            allfields = array[:,43]
            field = allfields[1]
        
            # ID the approx center of the group              
            medx = np.mean(array[:,2])
            medy = np.mean(array[:,3])
            
            # Find corrected number of elements in group
            nggd = len(array[:,0])/aread #Number of elements in group per deg^2
            nggdc = nggd - back # Number of elements per deg^2 corrected in group
            #print('Corrected elements in group per deg^2', nggdc, 'Group', i)
            ngg = nggdc*aread  # Number of corrected elements in group
            #print('Corrected elements in group', ngg, 'Group', i)
            
            # Find corrected number of elements per Mpc^2
            nareacm = nggdc/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find the CORRECTED number for background number of elements per Mpc^3
            nvcm = nareacm*(3/(4*c))
            
            # Find the mass per squared degree for this group
            md = massg/aread #Total mass of group/area of group in deg^2
            
            # Find CORRECTED mass per squared degree for this group
            mdc = md - backm #mgroup/deg^2 - mfield/deg^2 =  mgroup/deg^2 corrected
            
            # Find CORRECTED total mass of group
            mmmgc = mdc*aread 
        
            # Find CORRECTED mass per Mpc^2
            mdcm = mdc/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find CORRECTED mass per cubic megaparsecs
            mvcm = mdcm*(3/(4*c))
            
            '''
            # Find the number of objects in this group per deg2
            narea = len(array[:,0])/aread
            
            # Find the CORRECTED for background number of objects in this group per deg2 and Mpc^2
            nareac = narea - back #in degrees squared
            nareacm = nareac/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            
            # Find the CORRECTED number for background per Mpc^3
            nvcm = nareacm/(4/3*c)
        
            # Find the mass per squared degree for this group
            md = np.log10(mass)/aread
            
            # Find CORRECTED mass per squared degree for this group
            mdc = md - backm
            mdcm = mdc/(60**2)/(60**2)/(8.462**2)*(1000**2)  # in squared megaparsecs
            
            # Find corrected mass per cubic megaparsecs
            mvcm = mdcm/(4/3*c)
            '''
        
            fig = plt.figure(num=fieldw*20)
            ax = fig.add_subplot(111)
            pat = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='None', lw=1.3) 
            ax.add_patch(pat)
            
            centpx = np.mean(bbPath.vertices[:,0])
            centpy = np.mean(bbPath.vertices[:,1])

            print(centpx, centpy, i, file=clusters)
            
            cir = patchesbi.Circle((centpx, centpy), 17, facecolor='None', edgecolor='white')
            ax.add_patch(cir)
            
            if(ngroup==19)and(fieldw==1):
                ax.text(centpx+32, centpy-46, 'W'+str(fieldw)+'_'+str(ngroup), color='#3c0ff0', fontsize=10, fontweight='bold')
                ax.text(centpx+32, centpy-46, 'W'+str(fieldw)+'_'+str(ngroup), color='white', fontsize=10, fontweight='bold')
            
            elif(ngroup==14)and(fieldw==1):
                ax.text(centpx+35, centpy+26, 'W'+str(fieldw)+'_'+str(ngroup), color='#3c0ff0', fontsize=10, fontweight='bold')
                ax.text(centpx+35, centpy+26, 'W'+str(fieldw)+'_'+str(ngroup), color='white', fontsize=10, fontweight='bold')
            
            elif(fieldw==4)and(ngroup==7):
                
                ax.text(centpx+26, centpy-37, 'W'+str(fieldw)+'_'+str(ngroup), color='#3c0ff0', fontsize=10, fontweight='bold')
                ax.text(centpx+26, centpy-37, 'W'+str(fieldw)+'_'+str(ngroup), color='white', fontsize=10, fontweight='bold')
                
            else:
                ax.text(centpx+26, centpy+26, 'W'+str(fieldw)+'_'+str(ngroup), color='#3c0ff0', fontsize=10, fontweight='bold')
                ax.text(centpx+26, centpy+26, 'W'+str(fieldw)+'_'+str(ngroup), color='white', fontsize=10, fontweight='bold')
            
            
            '''
            ## If you want to plot individually each group
            '''            
            if(ngroup==14):#and(ngroup<=15): #for testing
                plot_groups_wide_wsf(fieldw, bbPath, ngroup, pw, scale, sw)
                tempx, tempy, tempng, tpatch, tbbpath, tggroup = color_plots_wide(fieldw, field, ngroup, centpx, centpy, scale, bbPath, ggroup, sw)#, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb)
                #plot_groups_wide(fieldw, bbPath, i, pw, scale, ngroup, gsmall)
           
            ## If you want to know the position of the brightest source in group
            brightest = array[array[:,17]==np.min(array[:,17])]
            #print('For field', fieldw, 'Group', i, 'Brightest gal has Ks', brightest[:,17])
            #print('Group', i, 'And position', brightest[:,4][0], brightest[:,5][0])
            #print('W'+str(fieldw)+'_'+str(i)+' '+str(brightest[:,4][0])+' '+str(brightest[:,5][0]))
            #print(brightest[:,4][0], brightest[:,5][0])
            
            print(fieldw, field, ngroup, medx, medy, np.log10(mmmgc), volume, nggdc, nareacm, nvcm, np.log10(mdc), np.log10(mdcm), np.log10(mvcm), aread, len(array[:,0]), ngg, nsf, ncsf, brightest[:,4][0], brightest[:,5][0], file=characw) 
           
            for element in array:
                
                overlapx = 968
                overlapy = 1290
             
                field = str(int(element[43]))
             
                if(fieldw==4)and((field==220930003100)or(field==221318002300)):
                    overlapy=1741
             
                shiftx, shifty = shift(fieldw, field)
    
                # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
                xarray = (element[2]+shiftx*(19354-overlapx))/scale
                yarray = (element[3]+shifty*(19354-overlapy))/scale
                #plt.figure(num=fieldw*10)
                #plt.scatter(xarray, yarray, marker='d', edgecolor='Blue', facecolor='None', s=10)
                
            
            #np.savetxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_'+str(int(fwhm*100))+'/updated_positions/test_group_'+str(ngroup)+'_'+str(limgroupw*100)+'%_w'+str(fieldw)+'_fwhm_'+str(int(fwhm*100))+'_l'+str(level)+'.dat', array, fmt='%s')
            
            ## New group search            
            ngroup = ngroup + 1
            
    return(tempx, tempy, tempng, tpatch, tbbpath, tggroup)
    #return(0,0,0,0,0,0)
    
#################################################################################
            
def plot_all_sf(sw, gsmall):
    
    for gal in sw: #All sfg with Ks < 20.5
    
    
        overlapx = 968
        overlapy = 1290
    
     
        field = str(int(gal[43]))
        
        if(fieldw==4)and((field==220930003100)or(field==221318002300)):
            overlapy=1741
            
        shiftx, shifty = shift(fieldw, field)
        
        # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
        xs = gal[2]+shiftx*(19354-overlapx)
        ys = gal[3]+shifty*(19354-overlapy)
        
        plt.scatter(xs/scale, ys/scale, marker='*', edgecolor='green', facecolor='None', s=100) #Real positions
        plt.scatter(xs/scale, ys/scale, marker='*', edgecolor='cyan', facecolor='None', s=50)
    

### REAL Data  

print('Starting!')  
fieldw = 4 #Choose the field you want to work with 
select = 20.5 #To show PEGs up to this mag. INDEPENDENT of massive PEGs selection

tempx = 0
tempy = 0
   
pw = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_w'+str(fieldw)+'_fullcat_nostr.dat') #artifacts, double cores removed from catalog
pw = pw[(pw[:,17]<20.5) & (pw[:,17]>18.5)] # second condition does not change number, or elements of groups
sw = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/sf_w'+str(fieldw)+'_fullcat_nostr.dat') 
sw = sw[sw[:,17]<20.5] #due to completeness in g
xytoradec = open('test.dat','w')
dense_file = open('test_170704_dense_wide'+str(fieldw)+'.dat','w')
test_contours2 = open('test_contours2.dat','w')
clusters = open('test_clusters.dat', 'w')

mastercat = open('test_mastercat2.dat','w')


ra = pw[:,4]
dec = pw[:,5]
ks = pw[:,17]
mass = 10**((-0.348)*ks+18.284)



# Find a suitable conversion between pix -> degrees. This transformation was calculated including the 100 scaling
mx, bx, my, by = xyradec(fieldw)

scale = 100 #to make things more manageable, scale everything by 

# The extra pixels secure a 5sigma extra padding for edge effects
slzeros = make_empty_array(fieldw)
print(slzeros.shape)

fwhm = 1000/scale # 632 pix or 0.032 degrees represents a box of side 1Mpc in size @ z~1.5  
sig = int(round(fwhm/2.355,0)) # fwhm = 2.355*sigma
print('I am working with a sigma of ', sig)


gsmall, xx, yy = makeGaussian2(5*sig, sig)
size = gsmall.shape[0]
print('Made gsmall! It has a shape of ', gsmall.shape)

if(fieldw==1):
    shape_fig = (14,5)
    edgy = 0
    
elif(fieldw==4):
    shape_fig = (12,8)
    edgy = 0.1
    
plt.figure(num=fieldw*20, figsize=shape_fig, dpi=100)
ax=plt.subplot(111)


plt.axis('scaled')
plt.xlim(((min(ra-edgy)-bx)/(mx)),((max(ra)-bx)/(mx)))
plt.ylim((min(dec)-by)/(my),(max(dec)-by)/(my))


#plt.title('PEGs W'+str(fieldw)+', gaussian density plot, Ks<20.5 fwhm= '+str(fwhm*100)+' pix') 
#plt.xlabel('x [pix]')
#plt.ylabel('y [pix]')


plt.xlabel('RA [deg]', fontsize=22)
plt.ylabel('DEC [deg]', fontsize=22)

print('Do you want to have points of all PEGs in your plots? yes/no:')
points = raw_input()

print('Do you want to have points of PEGs with Ks< '+str(select)+' in your plots? yes/no:')
selectp = raw_input()

print('Do you want to have points of all massive PEGS in your plots? yes/no:')
massive = raw_input()



if(massive=='yes'):
    print('Massive galaxies brighter than Ks?:')
    ksmassive = float(raw_input())
    
print('Do you want to select overdense regions? yes/no:')
overdense=raw_input()

#print('Do you want to show the position of groups? yes/no:')
#pgroups = raw_input()

if(overdense=='yes'):
    print('Overdense above? : (e.g. 2.7)')
    level = float(raw_input())
    
    print('Groups will be completed to what percentage of level? (e.g. 70)')
    limgroupw = float(raw_input())/100
    
    characw = open('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/smf/updated_positions/test_up_171016_characteristics_l'+str(level)+'_f'+str(fwhm*100)+'_c_'+str(limgroupw)+'_w'+str(fieldw)+'.dat','w')
    '''    
    print('# Groups in the Wide '+str(fieldw)+' fields selected with a fwhm of '+str(fwhm*100)+' and  a level of '+str(level), file = characw)
    print('# (0) Wide Field', file = characw)
    print('# (1) Tile/Patch', file = characw)
    print('# (2) Group', file = characw)
    print('# (3) X Center of group', file = characw)
    print('# (4) Y Center of group', file = characw)
    print('# (5) Mass [Log(M*/Msol)]', file = characw)
    print('# (6) Volume [Mpc^3] # Modeled using an ellipse', file = characw)
    print('# (7) Number Density [#/deg^2]', file = characw)
    print('# (8) Number Density [#/Mpc^2]', file = characw)
    print('# (9) Number Density [#/Mpc^3]', file = characw)
    print('# (10) Mass Area Density [Log(M*/Msol/deg^2)]', file = characw)
    print('# (11) Mass Area Density [Log(M*/Msol/Mpc^2)]', file = characw)
    print('# (12) Mass Density [Log(M*/Msol/Mpc^3)]', file = characw)
    print('# (13) Area of group [deg^2]', file = characw)
    print('# (14) Number of PEGs in group', file = characw)
    print('# (15) Number of PEGs -corrected for foreground/background contaminatio- in group', file = characw)
    print('# (16) Number of SFGs in group', file = characw)
    print('# (17) Number of SFGs -corrected for foreground/background contaminatio- in group', file = characw)
    print('# (18) RA [deg] for brightest item in group', file = characw)
    print('# (19) DEC [deg] for brightest item in group', file = characw)
    '''
    
    print('# Groups in the Wide '+str(fieldw)+' fields selected with a fwhm of '+str(fwhm*100)+' and  a level of '+str(level), file = characw)
    print('# (0) Field', file = characw)
    print('# (1) Tile/Patch', file = characw)
    print('# (2) Group', file = characw)
    print('# (3) X Center of group', file = characw)
    print('# (4) Y Center of group', file = characw)
    print('# (5) Corrected Total Mass [Log(M*/Msol)]', file=characw)
    print('# (6) Volume [Mpc^3] # Modeled using an ellipse', file = characw)
    print('# (7) Number Density [#/deg^2]', file = characw)
    print('# (8) Number Density [#/Mpc^2]', file = characw)
    print('# (9) Number Density [#/Mpc^3]', file = characw)
    print('# (10) Mass Area Density [Log(M*/Msol/deg^2)]', file = characw)
    print('# (11) Mass Area Density [Log(M*/Msol/Mpc^2)]', file = characw)
    print('# (12) Mass Density [Log(M*/Msol/Mpc^3)]', file = characw)
    print('# (13) Area of group [deg^2]', file = characw)
    print('# (14) Number of PEGs in group -NOT corrected for foreground/background contamination-', file = characw)
    print('# (15) Number of PEGs -corrected for foreground/background contaminatio- in group', file = characw)
    print('# (16) Number of SFGs in group', file = characw)
    print('# (17) Number of SFGs -corrected for foreground/background contaminatio- in group', file = characw)
    print('# (18) RA [deg] for brightest item in group', file = characw)
    print('# (19) DEC [deg] for brightest item in group', file = characw)
      
    print('Do you want complete groups? yes/no:')
    show_groups = raw_input()

    

for gal in pw: #All passive galaxies with Ks < 20.5
    
    
    overlapx = 968
    overlapy = 1290

    massgal = 10**((-0.348)*gal[17]+18.284)
    msy = 150*((massgal)/max(mass))
    
    
    field = str(int(gal[43]))
    if(fieldw==4)and((field==220930003100)or(field==221318002300)):
        overlapy=1741
        
    shiftx, shifty = shift(fieldw, field)
    
    # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
    x = gal[2]+shiftx*(19354-overlapx)
    y = gal[3]+shifty*(19354-overlapy)
    
    
    # Next shift every point by at least 1/2 gaussian so that you won't have edge effects
    xg = round(x/scale, 0)+round(gsmall.shape[0]/2, 0)
    yg = round(y/scale, 0)+round(gsmall.shape[0]/2, 0)
   
    # Plot all massive
    if(massive=='yes'):
        if(gal[17]<ksmassive):
            #plt.scatter(xg-int(gsmall.shape[0]/2), yg-int(gsmall.shape[0]/2), marker='s', edgecolor='darkred', facecolor='None', s=msy) #5 w1, 8 w4
            #plt.scatter(xg-int(gsmall.shape[0]/2), yg-int(gsmall.shape[0]/2), marker='s', edgecolor='gold', facecolor='None', s=0.5*msy) # 3 w1, 5 w4
            
            plt.scatter(x/scale, y/scale, marker='*', edgecolor='magenta', facecolor='None', s=2*msy) #Real positions
            plt.scatter(x/scale, y/scale, marker='*', edgecolor='None', facecolor='white', s=2*0.5*msy) # 
    
    # Plot ALL PEGs to plot on the pixel graph
    if(points=='yes'):
        #plt.scatter(xg-int(gsmall.shape[0]/2), yg-int(gsmall.shape[0]/2), marker='o', edgecolor='darkred', facecolor='None', s=msy)
        #plt.scatter(xg-int(gsmall.shape[0]/2), yg-int(gsmall.shape[0]/2), marker='o', edgecolor='gold', facecolor='None', s=0.5*msy)
    
        plt.scatter(x/scale, y/scale, marker='o', edgecolor='darkred', facecolor='None', s=msy) #Real positions
        plt.scatter(x/scale, y/scale, marker='o', edgecolor='gold', facecolor='None', s=0.5*msy)
    
    # Plot PEGs down to the selected Ks mag (which is normally a bit past the completeness of the sample)
    if(selectp=='yes'):
        if(gal[17]<select):
            
            plt.scatter(x/scale, y/scale, marker='*', edgecolor='darkred', facecolor='None', s=2*msy)
            plt.scatter(x/scale, y/scale, marker='*', edgecolor='gold', facecolor='None', s=2*0.5*msy)
    
            #plt.scatter(xg-int(gsmall.shape[0]/2), yg-int(gsmall.shape[0]/2), marker='*', edgecolor='darkred', facecolor='None', s=msy)
            #plt.scatter(xg-int(gsmall.shape[0]/2), yg-int(gsmall.shape[0]/2), marker='*', edgecolor='gold', facecolor='None', s=0.5*msy)
    
    # Write  the xy coordinates and corresponding RA, DEC of the object to make the unit conversion from pix to deg    
    print(gal[43], gal[1], xg-round(gsmall.shape[0]/2,0), yg-round(gsmall.shape[0]/2,0), gal[4], gal[5], file=xytoradec)    
    
    '''
    ### Gaussian Density Plot
    '''    
    if(gal[17]<20.5): #only make the gaussian density plots if Ks < 20.5   
        slzeros[xg-size/2:xg+size/2, yg-size/2:yg+size/2] = slzeros[xg-size/2:xg+size/2, yg-size/2:yg+size/2]+gsmall
   
## Plot ALL SFGs Ks<20.5 to make sure you are not missing anyone, visually inspect -> for testing only
#plot_all_sf(sw, gsmall)

####################### In pixels

### Plot also masked regions
#create_patches()


### Change axis from x,y to RA,DEC
mx, bx, my, by = xyradec(fieldw) #Find a suitable conversion between pix->deg both in RA and DEC    
tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis(fieldw, mx, bx, my, by)    
plt.xticks(tick_locsx, tick_lblsx, fontsize=18)
plt.yticks(tick_locsy, tick_lblsy, fontsize=18)
    

## Plot for comparison the Physical Scale of the Virgo Cluster
if(fieldw==1):
    plt.plot(np.arange(500, 525,0.01), 500*np.ones(len(np.arange(500, 525, 0.01))), linestyle='-', color='yellow')
    plt.text(600,530, 'Virgo Cluster', color='yellow', fontsize=12)
    
    plt.plot(np.arange(400, 496,0.01), 200*np.ones(len(np.arange(400, 496, 0.01))), linestyle='-', color='white')
    plt.text(500,220, '15 Mpc', color='white', fontsize=12)

if(fieldw==4):
    plt.plot(np.arange(850, 875,0.01), 300*np.ones(len(np.arange(850, 875, 0.01))), linestyle='-', color='yellow')
    plt.text(950,320, 'Virgo Cluster', color='yellow', fontsize=12)
    plt.plot(np.arange(400, 496,0.01), 200*np.ones(len(np.arange(400, 496, 0.01))), linestyle='-', color='white')
    plt.text(500,220, '15 Mpc', color='white', fontsize=12)

         
### Show gaussian density plot
plt.imshow(slzeros[round(gsmall.shape[0]/2,0):(slzeros.shape[0]-round(gsmall.shape[0]/2,0)), round(gsmall.shape[0]/2,0):(slzeros.shape[0]-round(gsmall.shape[0]/2,0))].T, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')
plt.colorbar()

#np.savetxt('gaussian_density2_w'+str(fieldw)+'.dat', slzeros[round(gsmall.shape[0]/2,0):(slzeros.shape[0]-round(gsmall.shape[0]/2,0)), round(gsmall.shape[0]/2,0):(slzeros.shape[0]-round(gsmall.shape[0]/2,0))].T, fmt='%s')


#### Select overdense regions 
if(overdense=='yes'): #select overdense regions
    
    cont_densew = plt.contour(slzeros[round(gsmall.shape[0]/2,0):(slzeros.shape[0]-round(gsmall.shape[0]/2,0)), round(gsmall.shape[0]/2,0):(slzeros.shape[0]-round(gsmall.shape[0]/2,0))].T, levels=[limgroupw*level, level], cmap=plt.get_cmap('Blues'), alpha=0.01)
    plt.clabel(cont_densew, fontsize=50)

    if(show_groups=='yes'):
    
        #select_groupsw(cont_densew, pw, fwhm, level, limgroupw, characw)
        tempx, tempy, tempng, tpatch, tbbpath, tggroup = select_groups_correctedw(cont_densew, pw, fwhm, level, limgroupw, fieldw, scale, sw, gsmall)
        




xytoradec.close()
dense_file.close()
characw.close()
mastercat.close()
clusters.close()

'''
if(show_groups=='yes'):
    data_groups = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/test_new_170714_characteristics_l'+str(level)+'_f'+str(fwhm*100)+'_wide'+str(fieldw)+'_c'+str(limgroupw)+'.dat')
    
    for iji in data_groups:
        
        overlapx = 968
        overlapy = 1290
    
        field = str(int(iji[1]))
        #print(field, type(field))
        #print('Info groups')
        #print('field', field)
        
        if(fieldw==4)and((field==220930003100)or(field==221318002300)):
            overlapy=1741
        
        shiftx, shifty = shift(fieldw, field)
        #print(overlapx, overlapy)
    
        # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
        xiji = (iji[3]+shiftx*(19354-overlapx))/scale
        yiji = (iji[4]+shifty*(19354-overlapy))/scale
        #print(xiji, yiji)
        
        plt.figure(num=fieldw*10)
        
        plt.text(xiji+0.004*xiji, yiji+0.007*yiji, 'Log[M*/Msol]='+str(round(iji[5],2)), color='blue', fontsize=12, bbox=dict(facecolor='yellow', alpha=0.5))
        ##plt.text(xiji+0.004*xiji, yiji-0.01*yiji, r'Log[M*/Msol]/Mpc$^3$='+str(round(iji[6],2)), color='blue', fontsize=10, bbox=dict(facecolor='yellow', alpha=0.5))
        ##plt.text(xiji-0.004*xiji, yiji-0.01*yiji, r'#/Mpc$^3$='+str(round(iji[12],2)), color='blue', fontsize=10, bbox=dict(facecolor='yellow', alpha=0.5))
            
        plt.text(xiji+0.005*xiji, yiji-0.01*yiji, r'Log[M*/Msol]/Mpc$^3$='+str(round(iji[12],2)), color='blue', fontsize=12, bbox=dict(facecolor='yellow', alpha=0.5))
        plt.text(xiji-0.005*xiji, yiji-0.01*yiji, r'#/Mpc$^3$='+str(round(iji[9],2)), color='blue', fontsize=12, bbox=dict(facecolor='yellow', alpha=0.5))
        
'''   
       

print('t_total=',timeit.time.time()-start, 's')
os.system('say "Liz python is done"')

       
