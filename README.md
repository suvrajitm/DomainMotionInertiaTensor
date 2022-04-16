# DomainMotionInertiaTensor
Maji, S. et al. 2017. Quantitative Characterization of Domain Motions in Molecular Machines. 
https://doi.org/10.1021/acs.jpcb.6b10732

Dependencies:
1. DIP image library (for segmentation of the 3D density maps, executables for segmentation code written in Matlab/DIP image are provided). Download and put it under the Utilities folder
2. LA library in Tcl (Linear algebra library already included in the Utilities folder) 

Note that although this code was uploaded in April 2022, it has not been updated since 2017. 

Pending issue:  Due to the random sign flip of eigenvectors , sometimes the axes sign switches direction from one pdb to another pdb (small changes in the residue position or mass distribution). Although a robust fix for the flips is available in the literature, the fix was not implemented back in 2016-2017, as I moved on to a different project.
Some simple criteria is there in the code to check for the sign flip, but it does not cover for all cases, so if there is flipped axes coordinate from the output files and then a custom change via the provided script for reversal of the sign of the coordinate axes is required. The adjusted axes should then be provide as input for calculating the motion between the sets of coordinate axes.
