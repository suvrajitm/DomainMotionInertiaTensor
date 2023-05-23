# DomainMotionInertiaTensor
Maji, S. et al. 2017. Quantitative Characterization of Domain Motions in Molecular Machines. 
https://doi.org/10.1021/acs.jpcb.6b10732

Dependencies:
1. To Perform Volume segmentation: DIP image processing library (for segmentation of the 3D density maps, executables for segmentation code written in Matlab/DIP image are provided). Download DIP toolbox and put it under the Utilities folder
2. Orient package in Tcl by Paul Grayson, See https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/ (included in the Utilities folder and a modified script)
4. LA library in Tcl (Linear algebra library already included in the Utilities folder) 


DomainMotionInertiaTensor.tcl is the starting script to run.

AbsoluteOrientation.tcl implements the procedures/functions and is the core of the motion characterization method along with the inertia tensor and principal axes calculation

You can set the domain definitions obtained either through segmentation or from a reference, in the script SetDomainSegments.tcl.

SegmentDensityMap.tcl is used for segmenting the map. There is a variable "mapSeglist_idx" in this script which is used to keep track of the segment indices. It can be set automatically but here it has been set manually for the data at hand with 4 different segments (0 - 3) to match the indices of the segmented domains. 


Notes: 
Pending issue:  Due to the random sign flip of eigenvectors , sometimes the axes sign flips in direction from one pdb to another pdb (small changes in the residue position or mass distribution). A fix for the sign flips is available in the literature, and needs to be implemented in the Orient script.
Some simple fix for sign flip is there but does not cover for all cases, so if there is flipped axes coordinate from the output files and then a custom change is required via the provided script for reversal of the sign of the coordinate axes. The adjusted axes should then be provide as input for calculating the motion between the sets of coordinate axes.
