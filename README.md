# DomainMotionInertiaTensor
Maji, S. et al. 2017. Quantitative Characterization of Domain Motions in Molecular Machines. 
https://doi.org/10.1021/acs.jpcb.6b10732

Dependencies:
1. DIP image library (for segmentation, executables for segmentation needs the libraries). Download and put it under the Utilities folder
2. LA library in Tcl (already included in the Utilities folder) 

Pending issue:  Due to the random eigenvector sign flipping, sometimes the axes sign switches direction from one pdb to another pdb. A robust fix for this is available in the literature, but did not have the time to implement that, since I moved on to a different project.

Simple criteria is there to check for the sign flip, but it does not cover for all cases, so if there is flipped axes coordinate from the output files and then a custom change via the provided script for reversal of the sign of the coordinates is required and provide as input for calculating the motion between the sets of coordinates.
