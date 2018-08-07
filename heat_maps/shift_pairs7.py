# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 12:01:16 2017

@author: osejo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.patches as patchesbi
import os
import timeit
import matplotlib.pyplot as plt

start=timeit.time.time()


'''
The difference between this script and shift_pairs.py is that here
everything is done in pixels, for consistency with real groups
'''

def allGroupsForElements(elements):
    
    allGroups = []
    
    for el in elements:
        if(el in dictElements):
            allGroups = np.r_[allGroups, dictElements[el]]
            del dictElements[el]
            
    return(allGroups)

###-----------------------------------------------------------------------------------------------------------------
    
def allElementsInGroups(group):
    
    allElements = []
    
    for item in group:
        allElements = np.r_[allElements, dictGroups[item]]
        
    return(allElements)
    
    
###-----------------------------------------------------------------------------------------------------------------


def mergeGroups(group, newGroup):
    
    if(len(group)==0):
        return
   

    fgroups = allGroupsForElements(allElementsInGroups(group))
    
    if(len(fgroups)==0):
        return
    
    mergeGroups(fgroups,newGroup)
    
    for  g in fgroups:
        if(g not in newGroup):
            newGroup.append(g)
    
###-------------------------------------------------------------------------------------------------------------------

def is_this_duplicate(final, obj):
    
    flag = True
    
    for el in final:
        flag = False
        
        if(el[1]==obj[1])and(el[2]==obj[2]):
            print('Duplicate')
            print(el[1], obj[1], el[2], obj[2])
            flag = True
            
    return(flag)
    
    
###--------------------------------------------------------------------------------------------------------------------
def show_groupings(final, cat):
    
    if(graphing==True):
        plt.figure(num=3, figsize= sizes)
        plt.scatter(final[:,12], final[:,13], marker = 'o', color='blue', alpha=0.5)
        
    for gal in cat:
        
        if(graphing==True):
            cir1 = patchesbi.Circle((gal[11], gal[12]), rsearch_pix, facecolor='None', edgecolor='cyan')
            plt.gca().add_patch(cir1)
    
    groups = np.unique(final[:,0])
    
    
    for el in groups:
        gg = final[final[:,0]==el] #select that group
        
        if(len(gg)>=1):        
            centx = np.mean(gg[:,12])
            centy = np.mean(gg[:,13])
            
            if(graphing==True):
                cir = patchesbi.Circle((centx, centy), rsearch_pix, facecolor='None', edgecolor='red')
                #plt.gca().add_patch(cir)
                #plt.text(centx, centy, str(el), color='black', fontsize=10, fontweight='bold')
            

###--------------------------------------------------------------------------------------------------------------------

def make_unique_ID(fullcatInitial):
    
    '''
    Assign a value from 1 to len(cat) to each galaxy that will be 
    its unique ID. This is due to the fact that each patch had an 
    individual ID so two different objects on different patches could
    end up with the same ID
    '''
    
    fullcat = np.c_[fullcatInitial, np.arange(1, len(fullcatInitial)+1,1)]
            
    return(fullcat)
###--------------------------------------------------------------------------------------------------------------------
def Add_Groups(gal, realnbs, i, final, cat_with_xyfull):
    
    '''
    From the list of neighbours (including the the galaxy itself with a distance = 0)
    Create an array that simply finds that index in the master catalog (cat_with_xyfull)
    and appends it to an array called final.
    
    i is a value that is given to each indiviual galaxy in the same group, it should start
    in 1 and it only increases once it passes from one group to another
    
    Returns cat_with_xyfull with an extra column at the beginning that assigns the value i
    to each group
    '''
        
    for nb in realnbs:
        
        if(len(final)==0):
            final = np.r_[i, cat_with_xyfull[int(nb[1])]]
        
        else:  
            final = np.vstack((final, np.r_[i,cat_with_xyfull[int(nb[1])]]))
        
    i = i + 1
    
    return(i, final)
###--------------------------------------------------------------------------------------------------------------------

def merge_groups(merged, cat, final):
    
    
    tmerged = merged#[135:137]#[135:137]#easy to see test #[137:144] #funny subsample ####select a subsample for testing;check indices on global variables
    farray = []
    idd_already = []
    
    
    for equal in tmerged: 
        
        if(len(equal)>1):
            #print('I have mroe than 1', equal)
            
            for item in equal:
                vec = []
                for gal in final:
                    if(gal[0]==item):
                        
                        vec = np.r_[equal[0], gal]
                        if(len(farray)==0):
                            farray = vec
                            idd_already.append(gal[11])
                            
                        if(len(farray)>0)and(gal[11] not in idd_already):
                            farray = np.vstack((farray, vec))
                            idd_already.append(gal[11])
                                        
                       
        else:   ### If there is only 1 object in group
            
            for gal in final:
                vec = []
                
                if(gal[0]==equal[0]):  ## Find that galaxy in the main catalog
                    vec = np.r_[equal[0], gal]         
                    
                    if(len(farray)==0):
                        farray = vec  
                        
                    else:
                        farray = np.vstack((farray, vec))
                    
                
    #output.close()
    #np.savetxt('test_merged2.cat', farray, fmt='%s')
    return(farray)

###--------------------------------------------------------------------------------------------------------------------

def show_merged_structures(merged, cat, final):
    
    cat_merged = merge_groups(merged, cat, final)
    unique_groups = np.unique(cat_merged[:,0])
    
    
    for g2 in unique_groups:
        
        iAmgroup = cat_merged[cat_merged[:,0]==g2]
        
        centx = np.mean(iAmgroup[:,13])
        centy = np.mean(iAmgroup[:,14])
        
        if(graphing==True):
            plt.figure(num=3, figsize= sizes)
            cir = patchesbi.Circle((centx, centy), rsearch_pix, facecolor='None', edgecolor='purple')
            #plt.gca().add_patch(cir)
            #plt.text(centx, centy, str(g2), color='purple', fontsize=10, fontweight='bold')
        
        if(graphing==True)and(len(iAmgroup)>2):
            plt.figure(num=3, figsize= sizes)
            cir = patchesbi.Circle((centx, centy), rsearch_pix, facecolor='None', edgecolor='magenta', lw=3)
            plt.gca().add_patch(cir)
            cir = patchesbi.Circle((centx, centy), rsearch_pix+500, facecolor='None', edgecolor='magenta')
            plt.gca().add_patch(cir)
    
    return(cat_merged)

###--------------------------------------------------------------------------------------------------------------------

def test_positions_mask(new_iAmgroupx, new_iAmgroupy):
    
    '''
    Based on the new random positions of objects in group
    Make sure that ALL elements in group, once shifted to
    new position (based on center of group) will be on good pixs
    '''
    
    allGood = False
    array = np.c_[new_iAmgroupx, new_iAmgroupy]
    #print('Checking flags in', array)
    
    tempflag = True
    
    for el in array:
        
        tempflag = tempflag*dwmask[int(round(el[0],0))][int(round(el[1],0))]
            
    allGood = tempflag  
    
    #print('Flag all 3', allGood)
    return(allGood)

    
###--------------------------------------------------------------------------------------------------------------------

def save_to_file(new_iAmgroup_x, new_iAmgroup_y, iAmgroup, final_catalog):
    '''
    Apend to catalog, as an extra 2 cols at the end the new randomized
    positions of small-scale structure
    '''

    array = np.c_[iAmgroup, new_iAmgroup_x, new_iAmgroup_y]
    #print(array)
    
    
    if(len(final_catalog)==0):
        final_catalog = array
        
    else:
        final_catalog = np.vstack((final_catalog, array))
    

    return(final_catalog)
###--------------------------------------------------------------------------------------------------------------------
  
def shuffle_positions(cat_merged, final_catalog):
    
    '''
    
    After IDing small scale structure, shuffle it to random positions in the field
    
    '''
    
    unique_groups = np.unique(cat_merged[:,0])
    
    for ugroup in unique_groups:
        #print('Working with group', ugroup)
        flag = False
        iAmgroup = cat_merged[cat_merged[:,0]==ugroup]
        
        centx = np.mean(iAmgroup[:,13])
        centy = np.mean(iAmgroup[:,14])
        
        #print('Real center group', centx, centy)
        
        if(graphing==True):
            plt.figure(num=3, figsize= sizes)
            plt.scatter(centx, centy, s = 80, marker='+', facecolors='purple', linewidths=4, alpha=0.7)
        
        ### Re-shuffle
        while(flag==False):
            
            new_random_groupCentx = np.random.randint(xmin, xmax)
            new_random_groupCenty = np.random.randint(ymin, ymax)
            #print('Random pos', new_random_groupCentx, new_random_groupCenty)
            
            
            transform_x = centx - new_random_groupCentx
            transform_y = centy - new_random_groupCenty
            #print('Transformation', transform_x, transform_y)
            
            new_iAmgroup_x = iAmgroup[:,13] - transform_x
            new_iAmgroup_y = iAmgroup[:,14] - transform_y
            
            #print('New position', new_iAmgroup_x, new_iAmgroup_y)
            
            flag = test_positions_mask(new_iAmgroup_x, new_iAmgroup_y)
        
        final_catalog = save_to_file(new_iAmgroup_x, new_iAmgroup_y, iAmgroup, final_catalog)
        
        if(graphing==True):
            plt.figure(num=2, figsize= sizes)
            plt.scatter(new_random_groupCentx, new_random_groupCenty, s = 80, marker='x', facecolors='purple', linewidths=4, alpha=0.7)
            plt.text(new_random_groupCentx, new_random_groupCenty, str(ugroup), color='green', fontsize=10, fontweight='bold')
            cir = patchesbi.Circle((new_random_groupCentx, new_random_groupCenty), 0.03, facecolor='None', edgecolor='purple')
            plt.gca().add_patch(cir)
            
            
            plt.figure(num=2, figsize= sizes)
            plt.scatter(new_iAmgroup_x, new_iAmgroup_y, s = 30, marker='s', facecolors='purple', linewidths=0.5, alpha=0.5)
            
            plt.figure(num=1, figsize= sizes)
            plt.scatter(new_iAmgroup_x, new_iAmgroup_y, s = 30, marker='s', facecolors='purple', linewidths=0.5, alpha=0.5)
               
   
    return(final_catalog)
###--------------------------------------------------------------------------------------------------------------------

def xyradec(fieldw):
    
    '''
    Find a suitable conversion x,y ->ra, dec
    
    '''
    
    if(os.path.isfile('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_'+str(fieldw)+'.dat')):
        data=np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/radec_xy/xyradec_'+str(fieldw)+'.dat')
    else:
        data = np.genfromtxt('/home/larcilao/nqs/fake_groups/small_scale_structure/xyradec_'+str(fieldw)+'.dat')
    ### This conversion, pixs are scaled by 100
    
    x = data[:,2]
    y = data[:,3]
    
    ra = data[:,4]
    dec = data[:,5]
    
    mx, bx = np.polyfit(x, ra, 1)
    my, by = np.polyfit(y, dec, 1)
    
    return(mx, bx, my, by)
###--------------------------------------------------------------------------------------------------------------------

def assign_pixvalues_full_field(catalog, field):
    
    '''
    This one is not being used for consistency with original data, see assign_pixvalues_full_field2()
    '''
    
    xfull = []
    yfull = []
    finc = []
    
    mx, bx, my, by = xyradec(field)
    xfull = (catalog[:,4]-bx)/mx
    yfull = (catalog[:,5]-by)/my
  
    
    plt.figure(num=1, figsize=sizes)
    plt.scatter(xfull, yfull, marker='o', color='blue', alpha=0.5)
    
     
 
###--------------------------------------------------------------------------------------------------------------------
   

def shift(fieldw, field):
    
    if(fieldw=='w1'):
        xpatch = {'20241041200':8, '20241050800':8,'20241060400':8, '20631041200':7, '20631050800':7, '20631060400':7, '20631070000':7, '21021041200':6, '21021050800':6, '21021060400':6, '21021070000':6, '21410041200':5, '21410050800':5, '21410060400':5, '21800041200':4, '21800050800':4, '21800060400':4, '22150041200':3, '22150050800':3, '22150060400':3, '22539050800':2, '22539060400':2, '22929041200':1, '22929050800':1, '22929060400':1, '23319041200':0, '23319050800':0, '23319060400':0}
        ypatch = {'20241041200':3, '20241050800':2,'20241060400':1, '20631041200':3, '20631050800':2, '20631060400':1, '20631070000':0, '21021041200':3, '21021050800':2, '21021060400':1, '21021070000':0, '21410041200':3, '21410050800':2, '21410060400':1, '21800041200':3, '21800050800':2, '21800060400':1, '22150041200':3, '22150050800':2, '22150060400':1, '22539050800':2, '22539060400':1, '22929041200':3, '22929050800':2, '22929060400':1, '23319041200':3, '23319050800':2, '23319060400':1}

        shiftx = xpatch[field]
        shifty = ypatch[field]
    
    elif(fieldw=='w4'):
        xpatch = {'220154011900':5, '220154021500':5,'220154031100':5, '220542011900':4, '220542021500':4, '220542031100':4, '220930003100':3, '220930002300':3, '220930011900':3, '220930021500':3, '221318003100':2, '221318002300':2, '221318011900':2, '221318021500':2, '221706002300':1, '221706011900':1, '221706021500':1, '222054002300':0, '222054011900':0}
        ypatch = {'220154011900':2, '220154021500':3,'220154031100':4, '220542011900':2, '220542021500':3, '220542031100':4, '220930003100':0, '220930002300':1, '220930011900':2, '220930021500':3, '221318003100':0, '221318002300':1, '221318011900':2, '221318021500':3, '221706002300':1, '221706011900':2, '221706021500':3, '222054002300':1, '222054011900':2}
    
        shiftx = xpatch[field]
        shifty = ypatch[field]
        
        
    return(shiftx, shifty)

 
###--------------------------------------------------------------------------------------------------------------------

def assign_pixvalues_full_field2(catalog, field, sizes):
    
    cat_with_xyfull = []
    
    for gal in catalog: #All passive galaxies with Ks < 20.5
    
        vec = []
        overlapx = 968
        overlapy = 1290
    
        patch = str(int(gal[0]))
        
        if(field=='w4')and((patch==220930003100)or(patch==221318002300)):
            overlapy=1741
        
        shiftx, shifty = shift(field, patch)
    
        # shift each position according to the position of their patch, each square degree is 19354 pixels minus the overlap between adjacent patches: 4 arcmin in DEC and 3 in RA
        xf = gal[2]+shiftx*(19354-overlapx)
        yf = gal[3]+shifty*(19354-overlapy)
        
        vec = np.r_[gal, xf, yf]
        
        if(len(cat_with_xyfull)==0):
            cat_with_xyfull = vec
        else:
            cat_with_xyfull = np.vstack((cat_with_xyfull, vec))
         
        ### For testing purposes
        #plt.figure(num=3, figsize=sizes)
        #plt.scatter(round(xf,0), round(yf,0), marker='^', color='red', alpha=0.5)
        
    return(cat_with_xyfull)


###--------------------------------------------------------------------------------------------------------------------
   
field = 'w1'
test = False #square patch small, not many flags
test_flags = False #test needs to be False for this one to work
graphing = True
its = 1
rsearch = 0.025 #in degrees, radius of search of nearest neighbours. Average size for groups is 0.025 deg
rsearch_pix =  rsearch*60*60/0.187 # Radius of search in pixels
nnbs = 10 # Number of nearest neighbours to return

dwmask = np.load('/Volumes/Liz/fake_env/pr/'+field+'_flip2.npy', mmap_mode='r')

m = 1
while(m<=its):
    
    print('Iteration', m)
    
    ####-------- Load files
    if(test==False):
        fullcatInitial = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/catalogs/pe_dw_cd_nostr_maglim_v_mu_2s01z05k_clean.dat') 
        sizes = (14,5)
    
    else:
        fullcatInitial = np.genfromtxt('/Users/osejo/Desktop/Tesis/massivend/k_catalog/automatize/0215/environment/sfandpe/density_plots/separate_populations/gaussian/dense_fake/pr/shift_pairs/w1_subset2.dat')
        sizes = (10,8)
        
    final = []
    final_catalog = []
    
    
    ####----- Define a dictionary wih both W1 and W4 catalogs
    dcat = {}
    w1cat = fullcatInitial[(fullcatInitial[:,0]>4) & (fullcatInitial[:,4]<40)]
    w4cat = fullcatInitial[(fullcatInitial[:,0]>4) & (fullcatInitial[:,4]>40)]    
    
    dcat['w1'] = w1cat
    dcat['w4'] = w4cat
    
    if(test_flags==True):
        dcat[field] = dcat[field][[(dcat[field][:,4]>34.5) & (dcat[field][:,4]<37.4)]]
        
    catalogIN = dcat[field]
    catalog = make_unique_ID(catalogIN) #Assign to each galaxy a number btw 1 and len(cat) which will be their unique ID
    
    # To each galaxy in this field, give it an x, y positon to the full combined field
    cat_with_xyfull =  assign_pixvalues_full_field2(catalog, field, sizes)    
    
    tree = cKDTree(zip(cat_with_xyfull[:,11], cat_with_xyfull[:,12]))
    

    i = 1 #Flag that assign groups together
    
    xmin = min(cat_with_xyfull[:,11])
    xmax = max(cat_with_xyfull[:,11])
    
    ymin = min(cat_with_xyfull[:,12])
    ymax = max(cat_with_xyfull[:,12])
    
    
        
    plt.figure(num=3, figsize=sizes)
    plt.xlim(xmax+100, xmin-100)
    plt.ylim(ymin-100, ymax+100)
        
        
        
    for gal in cat_with_xyfull:
        
        #print('Working with', gal[0], gal[1], gal[10], 'findig nbs')
        
        realnbs = []
        
        dist, indx = tree.query(np.hstack((gal[11], gal[12])), k = nnbs, distance_upper_bound=rsearch_pix) #Return xnbs + himself, upper_bound was found as an average radius for the real groups
    
        nbs = np.c_[dist, indx]
        realnbs = nbs[(nbs[:,0]!=np.inf)] # Gets rid of  objects outside of radius of search
        
        if(len(realnbs)>1)and(nnbs==1):
            print('Something weird', gal)
        #print('Found', len(realnbs), 'real nbs', realnbs[0])
        
        i, final = Add_Groups(gal, realnbs, i, final, cat_with_xyfull)
        
       
    
    show_groupings(final, cat_with_xyfull)
    
    ### Take care of the fact that one item can belong to two groups: Merge all groups with intersecting elements
    
    dictGroups = {}
    dictElements = {}
    merged = []
    
    
    ### Create 2 dictionaries: One has as the key the group and the other the elements.
    ### Complementary to each key is either the elements in that group or the groups that 
    ### element belongs to
    
    for gal in final:
        if(gal[0] not in dictGroups): #Dictionary of groups, Key: Groups
            dictGroups[gal[0]] = []
        dictGroups[gal[0]].append(gal[11])
        
        if(gal[11] not in dictElements): #Dictionary where the key are the element ID
            dictElements[gal[11]] = []
        dictElements[gal[11]].append(gal[0])
        
    #### Using both dictionaries, create a FOF or similar algorithm
    ### This algorithm returns a list of lists called merged
    ### Each list in merged tells you which groups should be merged
        
    for group in dictGroups:
        newGroup = []
        mergeGroups([group], newGroup)
        if(len(newGroup)!=0):
            merged.append(newGroup)
        
    ### With the information of which groups should be merged, 
    ### make a merged catalog, that has as its first column the real 
    ### group ID
    cat_merged = show_merged_structures(merged, cat_with_xyfull, final)
    
    '''
    final_catalog = shuffle_positions(cat_merged, final_catalog)
    
    if(graphing==True):
        plt.figure(num=1, figsize=sizes)
        plt.scatter(final_catalog[:,15], final_catalog[:,16], s = 30, marker='^', facecolors='black', linewidths=0.5, alpha=0.5)
        
    #np.savetxt('test_small_scstructure_pix_k1_it'+str(m)+'.dat', final_catalog, fmt='%s')
    '''
    m = m + 1
    
print('t_total=',timeit.time.time()-start, 's')
os.system('say "Liz python is done"')