# This folder contains most of the code used to identify over-dense regions of massive passive galaxies:

# density_[deep/wide]7.py = Creates heatmaps (modeling each galaxy as a gaussian) to identify over-dense regions of passive # galaxies.
# Over-dense regions are selected as a connected net of pixels above a certain established (by physical parameters) level.
# Once the code finds a peak above the given level, it finds its nearest contour from a contour map and derines a path for this contour (To identify the boundaries of the group).

# I calculate number densities and volumes in units of deg and physical Mpc. I correct for foreground/background contamination and output heatmaps, color images and combined catalogs of these over-densities
# Due to the difference in geometries betweeen the Deep and Wide surveys, the wide version includes my version of hopscotch to account for our survey being divided into overlapping tiles.

# density_dw_fake_rrr.py and shift_pairs7.py: To account for the effect of random projections and small scale clustering in our results
# we try to reproduce over-dense environments from heat maps using two kinds of fake catalogs:
# RRR = All galaxies are assigned a random position in the sky, calculating the probability of finding my groups by simple random projections. This was repeated 100 times.
# PR =  We assume that galaxies are more likely to be found in pairs or triplets and that and extra random projection could transform these small groups into detected over-dense regions.
# To define small-scale structure, we defined a linking length used to assume a galaxy to be gravitationally bound to another (based on research). After this linking length detected
# most likely pairs or triplets we performed a FoF analysis to link any remaining galaxies to this group and then randomized our sample maintaining established links. This procedure was 
# repeated 100 times.

# extras: Extra scripts create region files of the groups, color images, calculates number densities, etc.