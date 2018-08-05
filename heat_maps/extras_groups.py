# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 17:06:47 2017

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import os
import sys
import numpy as np
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as patchesbi
from astropy.io import fits
from PyAstronomy import pyasl
import matplotlib.path as path



def cutouts_wide(field, i, centpx, centpy, scale, f, fint):
    
    '''
    only needs to be run once, will create the .cl file necessary to create cut-outs
    of your UMPEGs
    '''
    
    overlapx = 968
    overlapy = 1290
             
    outfile = open('test_'+str(i)+'_w'+str(field)+'.cl','w')
    
    if(field==4)and((fint==220930003100)or(fint==221318002300)):
        overlapy=1741
                    
    shiftx, shifty = shift(field, fint)
    print('Shifts', shiftx, shifty)
    
    # de-shift each position so that they correspond to their positios in each patch
    xpatch = ((centpx*scale)-shiftx*(19354-overlapx))
    ypatch = ((centpy*scale)-shifty*(19354-overlapy))
    
    print('xpatch, ypatch', xpatch, ypatch)
    
    xmin = int(xpatch-1000) #10*scale was decided when cutting around the groups
    xmax = int(xpatch+1000)
    
    ymin = int(ypatch-1000)
    ymax = int(ypatch+1000)
    
    if(xmin<0):
        xmin = 1
    if(ymin<0):
        ymin = 1
            
    if(xmax>19354):
        xmax = 19354
    if(ymax>19354):
        ymax = 19354
            
    
    print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+str(field)+'/'+f+'/CFHTLS_W_Ks_'+f+'_T0007.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+str(field)+'_group_'+str(i)+'_ks.fits')
    print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+str(field)+'/'+f+'/CFHTLS_W_g_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+str(field)+'_group_'+str(i)+'_g.fits')
    print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+str(field)+'/'+f+'/CFHTLS_W_z_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+str(field)+'_group_'+str(i)+'_z.fits')
    
    print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+str(field)+'/'+f+'/CFHTLS_W_Ks_'+f+'_T0007.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+str(field)+'_group_'+str(i)+'_ks.fits', file=outfile)
    print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+str(field)+'/'+f+'/CFHTLS_W_g_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+str(field)+'_group_'+str(i)+'_g.fits', file=outfile)
    print('imcopy /home/larcilao/nqs/nqs_k_catalog/w'+str(field)+'/'+f+'/CFHTLS_W_z_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] w'+str(field)+'_group_'+str(i)+'_z.fits', file=outfile)
       
    

    
    outfile.close()
    return(i, xmin, xmax, ymin, ymax)
    
def makeFalseColor_wide(field, patch,  spatch, i, bandIndex=[0,1,2], limits=None):
    
    
    convFiles = ['/Volumes/Liz/color_groups/w'+field+'/w'+field+'_group_'+i+'_g.fits', '/Volumes/Liz/color_groups/w'+field+'/w'+field+'_group_'+i+'_z.fits', '/Volumes/Liz/color_groups/w'+field+'/w'+field+'_group_'+i+'_ks.fits']
 
    
    blue = fits.getdata(convFiles[bandIndex[0]])
    green = fits.getdata(convFiles[bandIndex[1]])
    red = fits.getdata(convFiles[bandIndex[2]])
    
    
    if limits is None:
        limits = [np.max(blue), np.max(green), np.max(red)]
        print ('Limits blue green red', limits)
        
       
    if(field=='1'):
        
        #red_dic = {'20241041200':[0,7000], '20241050800':[0,7000],'20241060400':[0,7000], '20631041200':[0,7000], '20631050800':[0,7000], '20631060400':[0,7000], '20631070000':[0,7000], '21021041200':[0,7000], '21021050800':[0,7000], '21021060400':[1000,7000], '21021070000':[100,130000], '21410041200':[0,7000], '21410050800':[0,7000], '21410060400':[0,7000], '21800041200':[0,7000], '21800050800':[0,7000], '21800060400':[0,7000], '22150041200':[0,7000], '22150050800':[0,7000], '22150060400':[0,7000], '22539050800':[0,7000], '22539060400':[0,7000], '22929041200':[0,7000], '22929050800':[0,7000], '22929060400':[0,7000], '23319041200':[0,7000], '23319050800':[0,7000], '23319060400':[1000,7000]}
        #green_dic = {'20241041200':[0,10000], '20241050800':[0,10000],'20241060400':[0,10000], '20631041200':[0,10000], '20631050800':[0,10000], '20631060400':[0,10000], '20631070000':[0,10000], '21021041200':[0,10000], '21021050800':[0,10000], '21021060400':[0,10000], '21021070000':[7,10000], '21410041200':[0,10000], '21410050800':[0,10000], '21410060400':[0,10000], '21800041200':[0,10000], '21800050800':[0,10000], '21800060400':[500,1000], '22150041200':[0,10000], '22150050800':[0,10000], '22150060400':[0,10000], '22539050800':[0,10000], '22539060400':[0,10000], '22929041200':[0,10000], '22929050800':[0,10000], '22929060400':[0,10000], '23319041200':[0,10000], '23319050800':[0,10000], '23319060400':[0,10000]}
        #blue_dic = {'20241041200':[0,8000], '20241050800':[0,8000],'20241060400':[0,8000], '20631041200':[0,8000], '20631050800':[0,8000], '20631060400':[0,8000], '20631070000':[0,8000], '21021041200':[0,8000], '21021050800':[0,8000], '21021060400':[0,8000], '21021070000':[0,8000], '21410041200':[0,8000], '21410050800':[0,8000], '21410060400':[0,8000], '21800041200':[0,8000], '21800050800':[0,8000], '21800060400':[0,8000], '22150041200':[0,8000], '22150050800':[0,8000], '22150060400':[0,8000], '22539050800':[0,8000], '22539060400':[0,8000], '22929041200':[0,8000], '22929050800':[0,8000], '22929060400':[0,8000], '23319041200':[0,8000], '23319050800':[0,8000], '23319060400':[0,8000]}
        
        red_dic = {'1':[550,7000], '2':[600,7000], '3':[600,7000], '4':[600,7000],'5':[600,7000],'6':[600,7000],'7':[600,7000],'8':[600,7000],'9':[600,7000],'10':[600,7000],'11':[600,7000],'12':[600,7000],'13':[600,7000],'14':[600,7000],'15':[600,7000],'16':[600,7000],'17':[600,7000],'18':[600,7000],'19':[600,7000]}
        green_dic = {'1':[5,10000], '2':[10, 10000], '3':[10,10000], '4':[10,10000], '5':[10,10000],'6':[10,10000],'7':[10,10000],'8':[10,10000],'9':[10,10000],'10':[10,10000],'11':[10,10000],'12':[10,10000],'13':[10,10000],'14':[10,10000],'15':[10,10000],'16':[10,10000],'17':[10,10000],'18':[10,10000],'19':[10,10000]}
        blue_dic = {'1':[0,8000], '2':[0, 8000], '3':[0,8000], '4':[0,8000], '5':[0,8000],'6':[0, 8000],'6':[0, 8000],'7':[0, 8000],'8':[0, 8000],'9':[0, 8000],'10':[0, 8000],'11':[0, 8000],'12':[0, 8000],'13':[0, 8000],'14':[0, 8000],'15':[0, 8000],'16':[0, 8000],'17':[0, 8000],'18':[0, 8000],'19':[0, 8000]}
        
    elif(field=='4'):
             
        red_dic = {'1':[550,7000], '2':[600,7000], '3':[600,7000], '4':[600,7000],'5':[600,7000],'6':[600,7000],'7':[600,7000],'8':[600,7000],'9':[600,7000],'10':[600,7000],'11':[600,7000],'12':[600,7000],'13':[600,7000],'14':[600,7000],'15':[600,7000]}
        green_dic = {'1':[5,10000], '2':[40, 10000], '3':[10,10000], '4':[10,10000], '5':[10,10000],'6':[10,10000],'7':[10,10000],'8':[10,10000],'9':[10,10000],'10':[10,10000],'11':[10,10000],'12':[10,10000],'13':[10,10000],'14':[10,10000],'15':[10,10000]}
        blue_dic = {'1':[0,8000], '2':[0, 8000], '3':[0,8000], '4':[0,8000], '5':[0,8000],'6':[0, 8000],'6':[0, 8000],'7':[0, 8000],'8':[0, 8000],'9':[0, 8000],'10':[0, 8000],'11':[0, 8000],'12':[0, 8000],'13':[0, 8000],'14':[0, 8000],'15':[0, 8000]}
       
    
    fieldstr = str(int(patch))
    
    
    #minr = red_dic[fieldstr][0]
    #maxr = red_dic[fieldstr][1]
    
    #ming = green_dic[fieldstr][0]
    #maxg = green_dic[fieldstr][1]
    
    #minb = blue_dic[fieldstr][0]
    #maxb = blue_dic[fieldstr][1]
    
    minr = red_dic[i][0]
    maxr = red_dic[i][1]
    
    ming = green_dic[i][0]
    maxg = green_dic[i][1]
    
    minb = blue_dic[i][0]
    maxb = blue_dic[i][1]
    
    falsecolor = np.zeros((blue.shape[0], blue.shape[1], 3), dtype=float)
    falsecolor2 = np.zeros((red.shape[0], red.shape[1], 3), dtype=float)
    
    falsecolor[:,:,0] = linear(red, minr, maxr)
    falsecolor[:,:,1] = linear(green, ming, maxg)
    falsecolor[:,:,2] = linear(blue, minb, maxb)

    falsecolor2[:,:,0] = linear(red, 0, limits[2]*0.7)
    #falsecolor2 = linear(red, 0, limits[2])
    
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/color_plots/deep/testing/6_r'+str(minr)+'_'+str(maxr)+'_g'+str(ming)+'_'+str(maxg)+'_b'+str(minb)+'_b'+str(maxb)+'.png')
    return (falsecolor, falsecolor2, minr, maxr, ming, maxg, minb, maxb)
    
    
def xyradecwide_cutouts_indgroup(fieldw, xmin, ymin, overlapx, shiftx, overlapy, shifty, patch):
    
    '''
    Find a suitable conversion x,y ->ra, dec for the color plots of wide
    
    '''
    print(xmin, ymin)
    
    data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/w'+str(fieldw)+'/'+str(patch)+'/k-k_'+str(patch)+'.cat')
    
    x = (data[:,1])
    y = (data[:,2])
    
    ra = data[:,25]
    dec = data[:,26]
    
    
    mx, bx = np.polyfit(x-xmin, ra, 1)
    my, by = np.polyfit(y-ymin, dec, 1)
    
    
    return(mx, bx, my, by)
    
def transform_axis_cutouts_wide_indgroup(fieldw, mx, bx, my, by, i):
    
    if(fieldw==1):
        
        tick_lblsx_dic = {1:[32.4, 32.45], 2:[32.75, 32.8], 3:[32.15,32.2], 4:[31.0], 5:[34.85, 34.9], 6:[32.55], 7:[35.4, 35.45], 8:[32.95, 33.0], 9:[35.0], 10:[32.2, 32.25], 11:[38.25], 12:[37.7], 13:[33.75], 14:[38.35], 15:[33.6], 16:[31.95], 17:[30.7], 18:[38.1], 19:[33.15, 33.2]}
        tick_lblsy_dic = {1:[-6.7, -6.65, -6.6], 2:[-6.55], 3:[-6.4], 4:[-6.35], 5:[-6.2, -6.15], 6:[-6.05, -6.0], 7:[-6.0, -5.95], 8:[-5.9], 9:[-5.9, -5.85], 10:[-5.7], 11:[-5.35], 12:[-5.3, -5.25], 13:[-5.2, -5.15], 14:[-5.0, -4.95], 15:[-4.85], 16:[-4.6], 17:[-4.55, -4.5], 18:[-4.4], 19:[-4.25]}
        
         
    elif(fieldw==4):
        
        tick_lblsx_dic = {1:[332.4, 332.45], 2:[332.85], 3:[330.2], 4:[331.85], 5:[330.1], 6:[334.0, 334.05], 7:[333.4, 333.45], 8:[333.6], 9:[330.7, 330.75], 10:[330.1, 330.15], 11:[331.75], 12:[333.9], 13:[332.4, 332.45], 14:[330.6, 330.65], 15:[330.15]}
        tick_lblsy_dic = {1:[1.0], 2:[1.0], 3:[1.25, 1.30], 4:[1.30], 5:[1.65], 6:[1.9], 7:[1.9, 1.95], 8:[2.1, 2.15], 9:[2.1, 2.15], 10:[2.2, 2.25], 11:[2.25, 2.3], 12:[2.3], 13:[2.3, 2.35], 14:[2.55], 15:[2.55, 2.6]}
          
    tick_lblsx = tick_lblsx_dic[i]
    tick_lblsy = tick_lblsy_dic[i]

    convx = (tick_lblsx-bx)/(mx)
    convy = (tick_lblsy-by)/(my)
    
    tick_lblsx = np.array(tick_lblsx)
    tick_locsx = convx.tolist()
    tick_lblsy = np.array(tick_lblsy)
    tick_locsy = convy.tolist()
    
    print((tick_locsx, tick_lblsx))
    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
        
def color_plots_wide(field, patch, i, centpx, centpy, scale, bbPath, ggroup, sw):#, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb):
        
    '''
    Make a gzK color plot of the Wide protoclusters
    '''
    patch = str(int(patch))
    
    if(field==1):
        
        spatch = '0'+patch[:5]+'-'+patch[5:]
        xlimits = {1:[2000, 0], 4:[1660, 400], 7:[2000, 10]}
        ylimits = {1:[0, 2000], 4:[400, 1600], 7:[10, 2000]}

            
    elif(field==4):
        
        xlimits = {1:[1700, 500], 2:[1197, 0], 4:[1660, 400], 7:[2000, 10]}
        ylimits = {1:[400, 1700], 2:[400,1580], 4:[400, 1600], 7:[10, 2000]}
        
        mark = '+'
        if(patch[6:]=='003100'):
            mark = '-'
            
        spatch = patch[:6]+mark+patch[6:]
       
        
         
    ### Make a cut-out of the group, only once
    print('Will make a cut-out of field', field, 'Group', i, 'Patch', spatch, patch)
        
    ## Run independently to get the cuts
    i, xmin, xmax, ymin, ymax = cutouts_wide(field, i, centpx, centpy, scale, spatch, patch)
    
    
    falsecolor, falsecolor2, minr, maxr, ming, maxg, minb, maxb = makeFalseColor_wide(str(field), patch,  spatch, str(i), bandIndex=[0,1,2], limits=None)
    
    fig, ax = plt.subplots(figsize=(10,10), facecolor='w')
    
    
    if(i in xlimits):
        minx = xlimits[i][0]
        maxx = xlimits[i][1]
        
        miny = ylimits[i][0]
        maxy = ylimits[i][1]

    else:
        minx = 1600
        maxx = 400
    
        miny = 400
        maxy = 1600
        
        
    if(i==3)or(i==4):
        ax.set_title(r'Group W'+str(field)+' '+str(i), size=27, y=1.01)
    else:
        ax.set_title(r'Group W'+str(field)+' '+str(i), size=27, y=1.03)
        
    ax.set_xlabel('RA [deg]', size=25)
    ax.set_ylabel('DEC [deg]', size=25)
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
    
    overlapx = 968
    overlapy = 1290
             
             
    
    if(field==4)and((patch==220930003100)or(patch==221318002300)):
        overlapy=1741
                    
    shiftx, shifty = shift(field, patch)
    print('Shifts', shiftx, shifty)
    
    # de-shift each position so that they correspond to their positios in each patch
    
    vpath = scale*(bbPath.vertices)
    vpathf = np.c_[(vpath[:,0]-shiftx*(19354-overlapx))-xmin, (vpath[:,1]-shifty*(19354-overlapy))-ymin]
    
   
    patn = patchesbi.PathPatch(path.Path(vpathf), facecolor='None', edgecolor='white', lw=2, alpha=0.5) 
    ax.add_patch(patn)
    
     
    ### X axis flipped-to match contour-images
    ax.scatter((ggroup[:,2])-xmin, (ggroup[:,3])-ymin, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
    
    for sfgal in sw:
        
        sfpatch = str(int(sfgal[43]))
        
        
        if(sfpatch==patch):
            ax.scatter(sfgal[2]-xmin, sfgal[3]-ymin, marker='*', edgecolor='cyan', facecolor='None', s=170, alpha=0.5)
        
        
       
    
    
    #ax.scatter(xall*scale-(xmin), yall*scale-(ymin), marker='o', edgecolor='gold', facecolor='None', s=3.7*msall, alpha=0.7) # Plot all PEGs
    #ax.scatter(xsfs*scale-(xmin), ysfs*scale-(ymin), marker='o', edgecolor='cyan', facecolor='None', s=marsf, alpha=0.5) # Plot all SFGs
    #ax.scatter(xsfsb*scale-(xmin), ysfsb*scale-(ymin), marker='*', edgecolor='cyan', facecolor='None', s=4*marsfb, alpha=0.5) # Plot Ks<20.5 SFGs
    
    
    # Test colors
    
    if(i==2)or(i==3):
        x_pos = 1400
    else:
        x_pos = 300
        
    #ax.text(x_pos, 600, 'Red: '+str(minr)+' '+str(maxr), color='white')
    #ax.text(x_pos, 700, 'Blue: '+str(minb)+' '+str(maxb), color='white')
    #ax.text(x_pos, 800, 'Green: '+str(ming)+' '+str(maxg), color='white')
    
    ## Thest these positions using a region file on the cut-out regions
    #region_files_totest(np.c_[xall*scale-(xmin+500), yall*scale-(ymin+500)])
    #print(np.c_[xall*scale-xmin, yall*scale-ymin])
    
    
    ## Complement objects that are in the adjacent tiles
    
    
    if(field==1)and(i==2):
        ax.scatter(898.191, 1068.68, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        ax.scatter(1291.89, 1109.09, marker='*', edgecolor='cyan', facecolor='None', s=170, alpha=0.5)
        ax.scatter(513.972, 1541.72, marker='*', edgecolor='cyan', facecolor='None', s=170, alpha=0.5)
    
    if(field==1)and(i==7):
        ax.scatter(1650.11, 201.384, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        ax.scatter(188.331, 629.947, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==1)and(i==16):
        ax.scatter(1593.51, 608.037, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
        
    if(field==1)and(i==19):
        ax.scatter(487.677, 1179.74, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==4)and(i==2):
        ax.scatter(112.471, 1021.07, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==4)and(i==3):
        ax.scatter(962.635, 420.729, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==4)and(i==4):
        ax.scatter(485.142, 533.506, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==4)and(i==7):
        ax.scatter(1784.84, 925.84, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        ax.scatter(1791.26, 1870.62, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==4)and(i==6):
        ax.scatter(1273.58, 1588.55, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
    
    if(field==4)and(i==10):
        ax.scatter(1432.02, 1451.7, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
    if(field==4)and(i==14):
        ax.scatter(1660.96, 1126.42, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.5)
        
        
    
        
        
          
    mx, bx, my, by = xyradecwide_cutouts_indgroup(field, xmin, ymin, overlapx, shiftx, overlapy, shifty, spatch) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis_cutouts_wide_indgroup(field, mx, bx, my, by, i)    
    plt.xticks(tick_locsx, tick_lblsx, fontsize=22)
    plt.yticks(tick_locsy, tick_lblsy, fontsize=22)
    

    ax.imshow(falsecolor, vmin=0, origin='lower')
    plt.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/color_groups/deep/color_d'+str(field)+'_group_'+str(i)+'.pdf', dpi=100)
    
    #return(0,0,0,0,0,0)
    return(centpx, centpy,i,patch,bbPath, ggroup)
    
    #return(0,0,0,0,0,0)
    
def region_files_totest(datacutout):
    
    '''
    #Test the position of objects when doing the color figures
    '''    
    
    reg=open('test_region.reg','w')
    print('# Region file format: DS9 version 4.1',file=reg)
    print('# Test position of objects in color images', file=reg)
    print('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',file=reg)
    print('image',file=reg)
    
    for n in datacutout: 
        
        if(n[0]>=0)and(n[1]>=0):
            print ('point('+str(n[0])+','+str(n[1])+') #point= circle color= green', file=reg)
            
         
    reg.close() 
    
def color_plots_deep(field, i, centpx, centpy, scale, bbPath, ggroup, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb, ngroup):
    
    '''
    Make a gzK color plot of the Deep protoclusters
    '''
    
    if(field==1):
        f = '022559-042940'
        maxxlim = 500
        minxlim = 1500
        minylim = 500
        maxylim = 1500

    elif(field==2):
        f = '100028+021230'
        maxxlim = 300
        minxlim = 1700
        minylim = 300
        maxylim = 1600
            
            
    elif(field==3):
        f = '141927+524056'
            
            
    elif(field==4):
        f = '221531-174356'
        maxxlim = 200
        minxlim = 1790
        minylim = 210
        maxylim = 1800
        
        
    
        
    ### Make a cut-out of the group, only once
    #print('Will make a cut-out of field', field, 'Group', i)
        
    ## Run independently to get the cuts
    i, xmin, xmax, ymin, ymax = cutouts_deep(field, i, centpx, centpy, scale, f)
    
    falsecolor, falsecolor2, minr, maxr, ming, maxg, minb, maxb = makeFalseColor(str(field), str(f),  str(i), bandIndex=[0,1,2], limits=None)
    
    fig, ax = plt.subplots(figsize=(10,10), facecolor='w')
    

    ax.set_title(r'Group D'+str(field)+'_'+str(ngroup), size=27)
    ax.set_xlabel('RA [deg]', size=25)
    ax.set_ylabel('DEC [deg]', size=25)
    ax.set_xlim(minxlim, maxxlim)
    ax.set_ylim(minylim, maxylim)
    
    vpath = scale*(bbPath.vertices)
    vpathf = np.c_[vpath[:,0]-xmin, vpath[:,1]-ymin]
    
    patn = patchesbi.PathPatch(path.Path(vpathf), facecolor='None', edgecolor='white', lw=2, alpha=0.5) 
    ax.add_patch(patn)
    
    #### Original
    #ax.scatter(ggroup[:,2]-xmin, ggroup[:,3]-ymin, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.7)
    
    #ax.scatter(xall*scale-(xmin), yall*scale-(ymin), marker='o', edgecolor='gold', facecolor='None', s=3.7*msall, alpha=0.7) # Plot all PEGs
    #ax.scatter(xsfs*scale-(xmin), ysfs*scale-(ymin), marker='o', edgecolor='cyan', facecolor='None', s=4.7*marsf, alpha=0.7) # Plot all SFGs
    #ax.scatter(xsfsb*scale-(xmin), ysfsb*scale-(ymin), marker='*', edgecolor='cyan', facecolor='None', s=4.7*marsfb, alpha=0.7) # Plot Ks<20.5 SFGs
    
    ### X axis flipped-to match contour-images
    ax.scatter(ggroup[:,2]-xmin, ggroup[:,3]-ymin, marker='*', edgecolor='gold', facecolor='None', s=230, alpha=0.7, lw=0.6)
    
    ax.scatter(xall*scale-(xmin), yall*scale-(ymin), marker='o', edgecolor='gold', facecolor='None', s=3.7*msall, alpha=0.7, lw=0.6) # Plot all PEGs
    ax.scatter(xsfs*scale-(xmin), ysfs*scale-(ymin), marker='o', edgecolor='cyan', facecolor='None', s=marsf, alpha=0.7, lw=0.4) # Plot all SFGs
    ax.scatter(xsfsb*scale-(xmin), ysfsb*scale-(ymin), marker='*', edgecolor='cyan', facecolor='None', s=4*marsfb, alpha=0.7, lw=0.4) # Plot Ks<20.5 SFGs
   
    # Test colors
    #ax.text(1000, 600, 'Red: '+str(minr)+' '+str(maxr), color='white')
    #ax.text(1000, 700, 'Blue: '+str(minb)+' '+str(maxb), color='white')
    #ax.text(1000, 800, 'Green: '+str(ming)+' '+str(maxg), color='white')
    
    ## Thest these positions using a region file on the cut-out regions
    #region_files_totest(np.c_[xall*scale-(xmin+500), yall*scale-(ymin+500)])
    #print(np.c_[xall*scale-xmin, yall*scale-ymin])
    
     
    mx, bx, my, by = xyradecdeep_cutouts(field, xmin, ymin) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis_cutouts(field, mx, bx, my, by)    
    plt.xticks(tick_locsx, tick_lblsx, fontsize=20)
    plt.yticks(tick_locsy, tick_lblsy, fontsize=20)
   

    ax.imshow(falsecolor, vmin=0, origin='lower')
    plt.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/color_groups/updated_positions/deep/color_d'+str(field)+'_group_'+str(ngroup)+'.pdf', dpi=100)
    
def cutouts_deep(field, i, centpx, centpy, scale, f):
    
    '''
    only needs to be run once, will create the .cl file necessary to create cut-outs
    of your UMPEGs
    '''
    
    outfile = open('test.cl','w')
    xmin = int(centpx*scale-10*scale) #10 was decided when cutting around the groups
    xmax = int(centpx*scale+10*scale)
    
    ymin = int(centpy*scale-10*scale)
    ymax = int(centpy*scale+10*scale)
    
    if(xmin<0):
        xmin = 1
    if(ymin<0):
        ymin = 1
            
    if(xmax>19354):
        xmax = 19354
    if(ymax>19354):
        ymax = 19354
            
    
   
    print('imcopy /Users/osejo/Desktop/Tesis/gzHK/photometry/wmask/WIRDS_Ks_'+str(f)+'_T0002.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+str(field)+'_group_'+str(i)+'_ks.fits', file=outfile)
    print('imcopy /Volumes/Liz/deep_images/CFHTLS_D-85_g_'+str(f)+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+str(field)+'_group_'+str(i)+'_g.fits', file=outfile)
    print('imcopy /Volumes/Liz/deep_images/CFHTLS_D-85_z_'+f+'_T0007_MEDIAN.fits['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+'] d'+str(field)+'_group_'+str(i)+'_z.fits', file=outfile)
        
    

    outfile.close()
    return(i, xmin, xmax, ymin, ymax)

def showgroups(fieldw):
    
    '''
    Make circles that show an approximate of the groups positions
    '''
    group = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/group_catalog_w'+str(fieldw)+'.dat')
    ngroups = np.unique(group[:,38])
    
    i = 0
    
    
def test():
    
    
    #paths = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/test_contours2,dat')
    gauss = np.genfromtxt('gaussian_density_d1.dat')
    plt.figure()
    plt.imshow(gauss[51:69, 49:67], cmap=plt.get_cmap('hot'), origin='lower', aspect='equal' )
    #plt.imshow(gauss[49:67, 51:69], cmap=plt.get_cmap('hot'), origin='lower', aspect='equal' )
    
    plt.figure()
    plt.imshow(gauss, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal' )

def xyradecdeep_cutouts(field, xmin, ymin):
    
    '''
    Find a suitable conversion x,y ->ra, dec
    
    '''
    scale = 100
    
    data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_d'+str(field)+'.dat')
    
    x = (data[:,2])-xmin
    y = (data[:,3])-ymin
    
    ra = data[:,4]
    dec = data[:,5]
    
    mx, bx = np.polyfit(x, ra, 1)
    my, by = np.polyfit(y, dec, 1)
    
    
    return(mx, bx, my, by)
    
def xyradecdeep(field):
    
    '''
    Find a suitable conversion x,y ->ra, dec
    
    '''
    scale = 100
    
    data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_d'+str(field)+'.dat')
    
    x = (data[:,2])/scale
    y = (data[:,3])/scale
    
    ra = data[:,4]
    dec = data[:,5]
    
    mx, bx = np.polyfit(x, ra, 1)
    my, by = np.polyfit(y, dec, 1)
    
    
    return(mx, bx, my, by)
    
def xyradecwide(fieldw):
    
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
    
def xyradecwide_cutouts(fieldw, xmin, ymin, overlapx, shiftx, overlapy, shifty):
    
    '''
    Find a suitable conversion x,y ->ra, dec
    
    '''
    print(xmin, ymin)
    scale = 100
    
    data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_w'+str(fieldw)+'.dat')
    
    x = (data[:,2])-xmin/scale #because cut-outs start at pix 0
    y = (data[:,3])-ymin/scale
    
    ra = data[:,4]
    dec = data[:,5]
    
    
    mx, bx = np.polyfit(x, ra, 1)
    my, by = np.polyfit(y, dec, 1)
    
    
    return(mx, bx, my, by)

    
def plot_groups_wide_old(fieldw, bbPath):
    
    '''
    Do a zoom-in figure of the groups found 
    in the wide fields
    '''
    
    
    gauss = np.genfromtxt('gaussian_density_w'+str(fieldw)+'.dat')
    
    #xcent = (np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))/2
    #ycent = (np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))/2
    
    fig = plt.figure(figsize=(12,12), dpi=100)
    plt.axis('scaled')
    ax2 = fig.add_subplot(111)
    plt.xlim(np.min(bbPath.vertices[:,0])-10, np.max(bbPath.vertices[:,0])+10)
    plt.ylim(np.min(bbPath.vertices[:,1])-10, np.max(bbPath.vertices[:,1])+10)
    #pat2 = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='#3c0ff0', lw=1.3) 
    #ax2.add_patch(pat2)
    ax2.imshow(gauss[int(np.min(bbPath.vertices[:,0])-10):int(np.max(bbPath.vertices[:,0])+10), int(np.min(bbPath.vertices[:,1])-10):int(np.max(bbPath.vertices[:,1])+10)], cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')

def transform_axis_cutouts(field, mx, bx, my, by):
    
    if(field==1):
        
        tick_lblsx = [36.7, 36.68, 36.66]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [-4.7, -4.68]
        convy = (tick_lblsy-by)/(my)   
        
    elif(field==2):
        
        
        tick_lblsx = [150]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [2.7]
        convy = (tick_lblsy-by)/(my)
        
    elif(field==3):
        
        tick_lblsx = [214.2, 214.4, 214.6, 214.8, 215, 215.2, 215.4]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [52.4, 52.5, 52.6, 52.7, 52.8, 52.9, 53, 53.1]
        convy = (tick_lblsy-by)/(my)
        
    elif(field==4):
        
        tick_lblsx = [333.97, 334, 334.03]
        #tick_lblsx = [333.4, 333.5, 333.6, 333.7, 333.8, 333.9, 334, 334.1, 334.2]

        convx = (tick_lblsx-bx)/(mx)
        #tick_lblsy = [-18.2, -18, -17.8, -17.6, -17.4]
        tick_lblsy = [-17.62, -17.65, -17.68]
        
        convy = (tick_lblsy-by)/(my)

    tick_lblsx = np.array(tick_lblsx)
    tick_locsx = convx.tolist()
    tick_lblsy = np.array(tick_lblsy)
    tick_locsy = convy.tolist()

    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
    
def transform_axis_cutouts_wide(fieldw, mx, bx, my, by):
    
    if(fieldw==1):
        
        tick_lblsx = [32.35, 32.4, 32.45]
        #tick_lblsx = np.arange(31,39,0.02) #for cut-outs
        
        tick_lblsy = [-6.7, -6.65, -6.6]
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
    
    print((tick_locsx, tick_lblsx))
    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
    
def transform_axis(field, mx, bx, my, by):
    
    if(field==1):
        
        tick_lblsx = [36.66, 36.68, 36.7]
        #tick_lblsx = np.arange(36.1, 36.8, 0.05)#[36.1, 36.2, 36.3, 36.4, 36.5, 36.6, 36.7, 36.8]
        convx = (tick_lblsx-bx)/(mx)
        #tick_lblsy = np.arange(-4.9, -4.1, 0.05)#[-4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1]
        tick_lblsy = [-4.72, -4.7, -4.68]
        convy = (tick_lblsy-by)/(my)   
        
    elif(field==2):
        
        #tick_lblsx = [149.8, 150, 150.2, 150.4, 150.6]
        tick_lblsx = [149.8, 149.9, 150, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6]
        convx = (tick_lblsx-bx)/(mx)
        #tick_lblsy = [1.8, 2.0, 2.2, 2.4, 2.6]
        tick_lblsy = [1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8]
        convy = (tick_lblsy-by)/(my)
        
    elif(field==3):
        
        tick_lblsx = [214.2, 214.4, 214.6, 214.8, 215, 215.2, 215.4]
        convx = (tick_lblsx-bx)/(mx)
        tick_lblsy = [52.4, 52.5, 52.6, 52.7, 52.8, 52.9, 53, 53.1]
        convy = (tick_lblsy-by)/(my)
        
    elif(field==4):
        
        tick_lblsx = [333.97, 334, 334.03]
        #tick_lblsx = [333.4, 333.5, 333.6, 333.7, 333.8, 333.9, 334, 334.1, 334.2]

        convx = (tick_lblsx-bx)/(mx)
        #tick_lblsy = [-18.2, -18, -17.8, -17.6, -17.4]
        tick_lblsy = [-17.62, -17.65, -17.68]
        
        convy = (tick_lblsy-by)/(my)

    tick_lblsx = np.array(tick_lblsx)
    tick_locsx = convx.tolist()
    tick_lblsy = np.array(tick_lblsy)
    tick_locsy = convy.tolist()

    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
    
def transform_axis_wide(fieldw, mx, bx, my, by):
    
    if(fieldw==1):
    
        tick_lblsx = np.arange(30,39,0.05) #for cut-outs
        tick_lblsy = np.arange(-7.0,-3.5,0.05) # for cut-outs
        
    elif(fieldw==4):
    
        tick_lblsx = np.arange(330,337,0.05)
        tick_lblsy = np.arange(-0.5,2.7,0.05)
        


    convx = (tick_lblsx-bx)/(mx)
    convy = (tick_lblsy-by)/(my)
    
    tick_lblsx = np.array(tick_lblsx)
    tick_locsx = convx.tolist()
    tick_lblsy = np.array(tick_lblsy)
    tick_locsy = convy.tolist()

    return (tick_locsx, tick_lblsx, tick_locsy, tick_lblsy)
    
def plot_groups_deep(field, bbPath, xall, yall, msall, xmass, ymass, msmass, i, tempx, tempy, ngroup):
    
    '''
    Do a zoom-in figure of the groups found 
    in the deep fields
    '''
    
    gauss = np.genfromtxt('gaussian_density_d'+str(field)+'.dat')
    mx, bx, my, by = xyradecdeep(field)
    
    fig, ax2 = plt.subplots(figsize=(10,10), facecolor='w')
    
    
    ax2.imshow(gauss, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')     
    pat2 = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='white', lw=2) #3c0ff0
    ax2.add_patch(pat2)
    
    
    ax2.set_xlabel('RA [deg]', fontsize=25)
    ax2.set_ylabel('DEC [deg]', fontsize=25)
    #ax2.set_title('D'+str(field)+'_'+str(round(np.mean(bbPath.vertices[:,0])*mx+bx,3))+'_'+str(round(np.mean(bbPath.vertices[:,1])*my+by,3)), fontsize=22)
    ax2.set_title('Group D'+str(field)+'_'+str(ngroup), fontsize=27)    
    
    ax2.scatter(xall, yall, marker='o', edgecolor='darkred', facecolor='None', s=5*msall) # Plot all PEGs
    ax2.scatter(xall, yall, marker='o', edgecolor='gold', facecolor='None', s=0.5*5*msall) 
    
    ax2.scatter(xmass, ymass, marker='*', edgecolor='darkred', facecolor='None', s=2.3*5*msmass) # Plot PEGs with Ks < 20.5 (default but is user-assigned)
    ax2.scatter(xmass, ymass, marker='*', edgecolor='gold', facecolor='None', s=2.3*0.5*5*msmass)
    
    mx, bx, my, by = xyradecdeep(field) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis(field, mx, bx, my, by)    
    plt.xticks(tick_locsx, tick_lblsx, fontsize=20)
    plt.yticks(tick_locsy, tick_lblsy, fontsize=20)
    
    ##ax2.set_xlim(np.max(bbPath.vertices[:,0])+7, np.min(bbPath.vertices[:,0])-7)
    ##ax2.set_ylim(np.min(bbPath.vertices[:,1])-7, np.max(bbPath.vertices[:,1])+7)
    
    centx = np.mean(bbPath.vertices[:,0])
    centy = np.mean(bbPath.vertices[:,1])
    
    cutouts = {'1':[5,5,5,5], '2':[7,7,7,6], '4':[8,8,8,8]}
    fieldstr = str(field)
    ax2.set_xlim(centx+cutouts[fieldstr][0], centx-cutouts[fieldstr][1])
    ax2.set_ylim(centy-cutouts[fieldstr][2], centy+cutouts[fieldstr][3])
    
    
    if(abs(np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))>tempx):
        tempx = (np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))
        
    if((np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))>tempy):
        tempy = (np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))
        
    #print(np.min(bbPath.vertices[:,0]), np.max(bbPath.vertices[:,0]), np.mean(bbPath.vertices[:,0]))
    #print(np.min(bbPath.vertices[:,1]), np.max(bbPath.vertices[:,1]), np.mean(bbPath.vertices[:,1]))
    print(tempx, tempy)
    plt.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/updated_positions/group_d'+str(field)+'_group_'+str(ngroup)+'.pdf', dpi=100)
    return(tempx, tempy)
    
def plot_groups_deep_wsf(field, bbPath, xall, yall, msall, xmass, ymass, msmass, i, tempx, tempy, xsfs, ysfs, xsfsb, ysfsb, ngroup):
    
    '''
    Do a zoom-in figure of the groups found 
    in the deep fields
    '''
    
    gauss = np.genfromtxt('gaussian_density_d'+str(field)+'.dat')
    mx, bx, my, by = xyradecdeep(field)
    
    fig, ax2 = plt.subplots(figsize=(10,10), facecolor='w')
    #fig, ax2 = plt.subplots(1, figsize=(10,10), dpi=100)
    
    
    ax2.imshow(gauss, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')     
    pat2 = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='white', lw=2) #3c0ff0
    ax2.add_patch(pat2)
    
    
    ax2.set_xlabel('RA [deg]', fontsize=25)
    ax2.set_ylabel('DEC [deg]', fontsize=25)
    #ax2.set_title('D'+str(field)+'_'+str(round(np.mean(bbPath.vertices[:,0])*mx+bx,3))+'_'+str(round(np.mean(bbPath.vertices[:,1])*my+by,3)), fontsize=22)
    ax2.set_title('Group D'+str(field)+'_'+str(ngroup), fontsize=27)    
    
    ax2.scatter(xall, yall, marker='o', edgecolor='darkred', facecolor='None', s=5*msall, linewidths=2) # Plot all PEGs
    ax2.scatter(xall, yall, marker='o', edgecolor='gold', facecolor='None', s=0.5*5*msall, linewidths=2) 
    
    ax2.scatter(xmass, ymass, marker='*', edgecolor='darkred', facecolor='None', s=2.3*5*msmass, linewidths=2) # Plot PEGs with Ks < 20.5 (default but is user-assigned)
    ax2.scatter(xmass, ymass, marker='*', edgecolor='gold', facecolor='None', s=2.3*0.5*5*msmass, linewidths=2)
    
    ax2.scatter(xsfs, ysfs, marker='o', edgecolor='blue', facecolor='None', alpha=0.6, s=1.5*40, linewidths=2) # Plot all SFGs
    ax2.scatter(xsfs, ysfs, marker='o', edgecolor='cyan', facecolor='None', alpha=0.6, s=1.5*20, linewidths=2)
    
    ax2.scatter(xsfsb, ysfsb, marker='*', edgecolor='blue', facecolor='None', alpha=0.6, s=3*90, linewidths=2) # Plot SFGs with Ks < 20.5 (default but is user-assigned)
    ax2.scatter(xsfsb, ysfsb, marker='*', edgecolor='cyan', facecolor='None', alpha=0.6, s=3*45, linewidths=2)
    
    mx, bx, my, by = xyradecdeep(field) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis(field, mx, bx, my, by)    
    plt.xticks(tick_locsx, tick_lblsx, fontsize=20)
    plt.yticks(tick_locsy, tick_lblsy, fontsize=20)
    
    #ax2.set_xlim(np.max(bbPath.vertices[:,0])+7, np.min(bbPath.vertices[:,0])-7)
    #ax2.set_ylim(np.min(bbPath.vertices[:,1])-7, np.max(bbPath.vertices[:,1])+7)
    
    centx = np.mean(bbPath.vertices[:,0])
    centy = np.mean(bbPath.vertices[:,1])
    
    
    cutouts = {'1':[7,7,7,7], '2':[8,8,8,7], '4':[9,9,9,9]}
    fieldstr = str(field)
    ax2.set_xlim(centx+cutouts[fieldstr][0], centx-cutouts[fieldstr][1])
    ax2.set_ylim(centy-cutouts[fieldstr][2], centy+cutouts[fieldstr][3])
    #ax2.set_xlim(centx+7, centx-7)
    #ax2.set_ylim(centy-7, centy+6)
    
    
    if(abs(np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))>tempx):
        tempx = (np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))
        
    if((np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))>tempy):
        tempy = (np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))
        
    #print(np.min(bbPath.vertices[:,0]), np.max(bbPath.vertices[:,0]), np.mean(bbPath.vertices[:,0]))
    #print(np.min(bbPath.vertices[:,1]), np.max(bbPath.vertices[:,1]), np.mean(bbPath.vertices[:,1]))
    print(tempx, tempy)
    plt.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/updated_positions/group_sfpe_d'+str(field)+'_group_'+str(ngroup)+'.pdf', dpi=100)
    return(tempx, tempy)

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


def plot_groups_wide(fieldw, bbPath, i, pw, scale, ngroup, gsmall):
    
    '''
    Do a zoom-in figure of the groups found 
    in the wide fields
    '''
    print('Inside plot groups!')
    print('Group i', i)
    
    if(fieldw==1):
        figsizes = {1:10, 2:7, 3:7, 4:7, 5:7, 6:7, 7:10, 8:7, 9:7, 10:7, 11:7, 12:7, 13:7, 14:7, 15:7, 16:7, 17:7, 18:7, 19:7}
        ssymbol = 2
    else: 
        figsizes = {1:7, 2:7, 3:7, 4:7, 5:7, 6:7, 7:10, 8:7, 9:7, 10:7, 11:6, 12:7, 13:7, 14:7, 15:7, 16:7, 17:7, 18:7, 19:7}
        ssymbol = 7
        
    gauss = np.genfromtxt('gaussian_density2_w'+str(fieldw)+'.dat')
    
    fig, ax2 = plt.subplots(1, figsize=(10,10), dpi=100)
    
    
    ax2.imshow(gauss, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')     
    pat2 = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='white', lw=2) 
    ax2.add_patch(pat2)
    
    
    ax2.set_xlabel('RA [deg]', fontsize=27)
    ax2.set_ylabel('DEC [deg]', fontsize=27) 
    
    ### Change axis from x,y to RA,DEC
    
    mx, bx, my, by = xyradecwide(fieldw) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis_wide(fieldw, mx, bx, my, by)    
    plt.xticks(tick_locsx, tick_lblsx, fontsize=25)
    plt.yticks(tick_locsy, tick_lblsy, fontsize=25)
    
    ##ax2.set_xlim(np.max(bbPath.vertices[:,0])+7, np.min(bbPath.vertices[:,0])-7)
    ##ax2.set_ylim(np.min(bbPath.vertices[:,1])-7, np.max(bbPath.vertices[:,1])+7)
    
    ax2.set_title('W'+str(fieldw)+'_'+str(ngroup), fontsize=33)
    
    centx = np.mean(bbPath.vertices[:,0])
    centy = np.mean(bbPath.vertices[:,1])
    
    ax2.set_xlim(centx+figsizes[ngroup], centx-figsizes[ngroup])
    ax2.set_ylim(centy-figsizes[ngroup], centy+figsizes[ngroup])
    
    mass = 10**((-0.348)*pw[:,17]+18.284)
    #mass2 = np.sort(mass)
    #secmass = mass2[-2] # Second largest mass, not used at the moment
    
    if(fieldw==1):
        smarker = 870
    elif(fieldw==4):
        smarker = 400
    
    for gal in pw: #All passive galaxies with Ks < 20.5
    
        fieldg = str(int(gal[43]))
        overlapx = 968
        overlapy = 1290

        massgal = 10**((-0.348)*gal[17]+18.284)
        msy = smarker*((massgal)/np.max(mass)) #230 for w4 700 w1
    
        
        if(fieldw==4)and((fieldg==220930003100)or(fieldg==221318002300)):
            overlapy=1741
        
        shiftx, shifty = shift(fieldw, str(int(fieldg)))
    
        # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
        x = gal[2]+shiftx*(19354-overlapx)
        y = gal[3]+shifty*(19354-overlapy)
        #print('Pongo', x/scale, y/scale)
        plt.scatter(x/scale, y/scale, marker='*', edgecolor='darkred', facecolor='None', s=ssymbol*msy) #Real positions
        plt.scatter(x/scale, y/scale, marker='*', edgecolor='gold', facecolor='None', s=ssymbol*0.5*msy)
    
    fig.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/updated_positions/group_w'+str(fieldw)+'_group_'+str(ngroup)+'.pdf', dpi=100, bbox_inches='tight')
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/test_'+str(i)+'.pdf', dpi=100, bbox_inches='tight')
    
    
    '''
    
    ## To ensure all images will have the same size
    if(abs(np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))>tempx):
        tempx = (np.max(bbPath.vertices[:,0])-np.min(bbPath.vertices[:,0]))
        
    if((np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))>tempy):
        tempy = (np.max(bbPath.vertices[:,1])-np.min(bbPath.vertices[:,1]))
        
    #print(np.min(bbPath.vertices[:,0]), np.max(bbPath.vertices[:,0]), np.mean(bbPath.vertices[:,0]))
    #print(np.min(bbPath.vertices[:,1]), np.max(bbPath.vertices[:,1]), np.mean(bbPath.vertices[:,1]))
    print(tempx, tempy)
    plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/group_d'+str(field)+'_group_'+str(i)+'.pdf', dpi=100)
    return(tempx, tempy)
    
    '''
    
    #fig = plt.figure(num=fieldw*10)
    #plt.text('W'+str(fieldw)+'_', centx, centy+10)
    #return(0,0)
    
def plot_groups_wide_wsf(fieldw, bbPath, i, pw, scale, sw):
    
    '''
    Do a zoom-in figure of the groups found 
    in the wide fields
    '''
    print('Inside plot groups!')
    print('Group i', i)
    
    if(fieldw==1):
        figsizes = {1:10, 2:7, 3:7, 4:7, 5:7, 6:7, 7:10, 8:7, 9:7, 10:7, 11:7, 12:7, 13:7, 14:7, 15:7, 16:7, 17:7, 18:7, 19:7}
        ssymbol = 2
    else: 
        figsizes = {1:7, 2:7, 3:7, 4:7, 5:7, 6:7, 7:10, 8:7, 9:7, 10:7, 11:6, 12:7, 13:7, 14:7, 15:7, 16:7, 17:7, 18:7, 19:7}
        ssymbol = 7
  
    
    gauss = np.genfromtxt('gaussian_density2_w'+str(fieldw)+'.dat')
    
    fig, ax2 = plt.subplots(1, figsize=(10,10), dpi=100)
    
    
    ax2.imshow(gauss, cmap=plt.get_cmap('hot'), origin='lower', aspect='equal')     
    pat2 = patchesbi.PathPatch(bbPath, facecolor='None', edgecolor='white', lw=2) 
    ax2.add_patch(pat2)
    
    
    ax2.set_xlabel('RA [deg]', fontsize=27)
    ax2.set_ylabel('DEC [deg]', fontsize=27) 
    
    ### Change axis from x,y to RA,DEC
    
    mx, bx, my, by = xyradecwide(fieldw) #Find a suitable conversion between pix->deg both in RA and DEC    
    tick_locsx, tick_lblsx, tick_locsy, tick_lblsy = transform_axis_wide(fieldw, mx, bx, my, by)    
    #plt.xticks(tick_locsx, tick_lblsx, fontsize=25)
    #plt.yticks(tick_locsy, tick_lblsy, fontsize=25)
    
    ##ax2.set_xlim(np.max(bbPath.vertices[:,0])+7, np.min(bbPath.vertices[:,0])-7)
    ##ax2.set_ylim(np.min(bbPath.vertices[:,1])-7, np.max(bbPath.vertices[:,1])+7)
    
    ax2.set_title('W'+str(fieldw)+'_'+str(i), fontsize=33)
    
    centx = np.mean(bbPath.vertices[:,0])
    centy = np.mean(bbPath.vertices[:,1])
    
    ax2.set_xlim(centx+figsizes[i], centx-figsizes[i])
    ax2.set_ylim(centy-figsizes[i], centy+figsizes[i])
    
    mass = 10**((-0.348)*pw[:,17]+18.284)
    #mass2 = np.sort(mass)
    #secmass = mass2[-2] # Second largest mass, not used at the moment
    
    if(fieldw==1):
        smarker = 870
    elif(fieldw==4):
        smarker = 400
    
    for gal in pw: #All passive galaxies with Ks < 20.5
    
        fieldg = str(int(gal[43]))
        overlapx = 968
        overlapy = 1290

        massgal = 10**((-0.348)*gal[17]+18.284)
        msy = smarker*((massgal)/np.max(mass)) #230 for w4 700 w1
    
        
        if(fieldw==4)and((fieldg==220930003100)or(fieldg==221318002300)):
            overlapy=1741
        
        shiftx, shifty = shift(fieldw, str(int(fieldg)))
    
        # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
        x = gal[2]+shiftx*(19354-overlapx)
        y = gal[3]+shifty*(19354-overlapy)
    
        plt.scatter(x/scale, y/scale, marker='*', edgecolor='darkred', facecolor='None', s=ssymbol*msy) #Real positions
        plt.scatter(x/scale, y/scale, marker='*', edgecolor='gold', facecolor='None', s=ssymbol*0.5*msy)

    
    for gals in sw: #All sf galaxies with Ks < 20.5
    
        fieldgs = str(int(gals[43]))
        overlapx = 968
        overlapy = 1290


        if(fieldw==4)and((fieldgs==220930003100)or(fieldgs==221318002300)):
            overlapy=1741
        
        shiftx, shifty = shift(fieldw, str(int(fieldgs)))
    
        # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
        xs = gals[2]+shiftx*(19354-overlapx)
        ys = gals[3]+shifty*(19354-overlapy)
    
        plt.scatter(xs/scale, ys/scale, marker='*', edgecolor='blue', facecolor='None', s=200) #Real positions
        plt.scatter(xs/scale, ys/scale, marker='*', edgecolor='cyan', facecolor='None', s=150)

    
    fig.tight_layout()
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/group_w'+str(fieldw)+'_group_sfpe_'+str(i)+'.pdf', dpi=100, bbox_inches='tight')
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/fwhm_1000/updated_positions/group_w'+str(fieldw)+'_group_sfpe_'+str(i)+'.pdf', dpi=100, bbox_inches='tight')
    
    
def makeFalseColor(df, fieldi, groupi, bandIndex=[0,1,2], limits=None):
    
    
    
    convFiles = ['/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/color_groups/deep/d'+df+'_group_'+groupi+'_g.fits',
                 '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/color_groups/deep/d'+df+'_group_'+groupi+'_z.fits',
                 '/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense/140916/color_groups/deep/d'+df+'_group_'+groupi+'_ks.fits']
   
    
    blue = fits.getdata(convFiles[bandIndex[0]])
    green = fits.getdata(convFiles[bandIndex[1]])
    red = fits.getdata(convFiles[bandIndex[2]])
    
    
    if limits is None:
        limits = [np.max(blue), np.max(green), np.max(red)]
        print ('Limits blue green red', limits)
        
    red_dic = {'1':[5,5000], '2':[12,5000], '4':[0,7000]}
    green_dic = {'1':[0, 10000], '2':[0, 10000], '4':[0, 10000]}
    blue_dic = {'1':[0, 8000], '2':[50, 5000], '4':[0, 8000]}    
    
    fieldstr = str(df)
    
    minr = red_dic[fieldstr][0]
    maxr = red_dic[fieldstr][1]
    
    ming = green_dic[fieldstr][0]
    maxg = green_dic[fieldstr][1]
    
    minb = blue_dic[fieldstr][0]
    maxb = blue_dic[fieldstr][1]
    
    falsecolor = np.zeros((blue.shape[0], blue.shape[1], 3), dtype=float)
    falsecolor2 = np.zeros((red.shape[0], red.shape[1], 3), dtype=float)
    
    falsecolor[:,:,0] = linear(red, minr, maxr)
    falsecolor[:,:,1] = linear(green, ming, maxg)
    falsecolor[:,:,2] = linear(blue, minb, maxb)

    falsecolor2[:,:,0] = linear(red, 0, limits[2]*0.7)
    #falsecolor2 = linear(red, 0, limits[2])
    
    #plt.savefig('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/Tal12/1116/color_plots/deep/testing/6_r'+str(minr)+'_'+str(maxr)+'_g'+str(ming)+'_'+str(maxg)+'_b'+str(minb)+'_b'+str(maxb)+'.png')
    return (falsecolor, falsecolor2, minr, maxr, ming, maxg, minb, maxb)
    
def linear(inputArray, scale_min=None, scale_max=None):
	"""Performs linear scaling of the input numpy array.

	@type inputArray: numpy array
	@param inputArray: image data array
	@type scale_min: float
	@param scale_min: minimum data value
	@type scale_max: float
	@param scale_max: maximum data value
	@rtype: numpy array
	@return: image data array
	
	"""		
	#print "img_scale : linear"
	imageData=np.array(inputArray, copy=True)
	
	if scale_min == None:
		scale_min = imageData.min()
	if scale_max == None:
		scale_max = imageData.max()

	imageData = imageData.clip(min=scale_min, max=scale_max)
	imageData = 100*(imageData -scale_min) / (scale_max - scale_min)
	indices = np.where(imageData < 0)
	imageData[indices] = 0.0
	indices = np.where(imageData > 1)
	imageData[indices] = 1.0
	
	return imageData

def fullcircle():
    
    centpx = 0
    centpy = 0
    
    patchesbi.Circle((centpx, centpy), 20)
    
    #cir = patchesbi.Circle(center=(centpx, centpy), radius=20)

def conversion():
     
     test = "10 00 56.96 +02 20 09.32"
     #test = "00 05 08.83239 +67 50 24.0135"
     #test = "02 25 59.665 -4 29 30.98"
     print('Original', test)
     
     ra, dec = pyasl.coordsSexaToDeg(test)
     print("Coordinates in [deg]: %010.6f  %+09.6f" % (ra, dec))
     
     # Convert back into sexagesimal representation
     sexa = pyasl.coordsDegToSexa(ra, dec)
     print("Coordinates back in [sexa]: ", sexa)
     

    
#color_plots_deep(field, i, centpx, centpy, scale, bbPath, ggroup, array, xall, yall, msall, xmass, ymass, msmass, xsfs, ysfs, xsfsb, ysfsb, marsf, marsfb, 1)
#test()
#fullcircle()
#conversion()
#color_plots_wide(fieldw, tpatch, tempng, tempx, tempy, scale, tbbpath, tggroup, sw)
