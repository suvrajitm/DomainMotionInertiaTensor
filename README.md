# DomainMotionInertiaTensor
1. **Suvrajit Maji**, Rezvan Shahoei, Klaus Schulten, Joachim Frank. JPCB. 2017. Quantitative Characterization of Domain Motions in Molecular Machines. 
https://doi.org/10.1021/acs.jpcb.6b10732

2. Sayan Bhattacharjee, Xiangsong Feng, **Suvrajit Maji**, Prikshat Dadhwal, Zhening Zhang, Zuben P Brown, Joachim Frank. Cell. 2024. Time resolution in cryo-EM using a PDMS-based microfluidic chip assembly and its application to the study of HflX-mediated ribosome recycling. https://doi.org/10.1016/j.cell.2023.12.027


Dependencies:
1. To Perform Volume segmentation: DIP image processing library (for segmentation of the 3D density maps, executables for segmentation code written in Matlab/DIP image are provided). Download DIP toolbox and put it under the Utilities folder
2. Orient package in Tcl by Paul Grayson, See https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/ (included in the Utilities folder and a modified script)
4. LA library in Tcl (Linear algebra library already included in the Utilities folder) 

DomainMotionInertiaTensor.tcl is the starting script to run.

AbsoluteOrientation.tcl implements the procedures/functions and is the core of the motion characterization method along with the inertia tensor and principal axes calculation

You can set the domain definitions obtained either through segmentation or from a reference, in the script SetDomainSegments.tcl.

SegmentDensityMap.tcl is used for segmenting the map. There is a variable "mapSeglist_idx" in this script which is used to keep track of the segment indices. It can be set automatically but here it has been set manually for the data at hand with 4 different segments (0 - 3) to match the indices of the segmented domains. 

Few selected updates 2022-2023:
1. Determine the Hinge point through which the rotation axis and angle passes.
   Derivations and implementation in the Cell paper 2024:
   "Bhattacharjee, Sayan; Feng, Xiangsong; Maji, Suvrajit; Dadhwal, Prikshat; Zhang, Zhening; Brown, Zuben P.; Frank, Joachim, Time
   Resolution in Cryo-EM Using a Novel PDMS-Based Microfluidic Chip Assembly and Its Application to the Study of HflX-Mediated Ribosome 
   Recycling."
3. Option to align the reference domains between two pdbs using rmsd fit in vmd, since sometimes the alignment with just the three    
   principal axes is not good enough. The variable "align_force_tarmol_to_refmoldom" in ComputeIntertiaTensor.tcl can be set to 1 (rmsd 
   fit of target to reference domain) or 0 (fit of target to reference domain using the three principal axes) 
5. Manually flip the principal axes directions by setting the variable "man_flip1_a1a2a3" with the triplet values of 1, or -1 for each axis 
   direction, so that the     flips are also reflected in the visualized axes for the pdb structures on vmd, earlier we could only 
   calculate the rotation correctly after fixing the flip manually but the generated images would still be incorrect.  
  

Notes: 
Pending issue:  Due to the random sign flip of eigenvectors , sometimes the axes sign flips in direction from one pdb to another pdb (small changes in the residue position or mass distribution). A fix for the sign flips is available in the literature, and needs to be implemented in the Orient script.
Some simple manual fix for sign flip is there (Utilities/orient.tcl script and ComputeIntertiaTensor.tcl script) but does not cover for all cases, so if there is flipped axes coordinate from the output files and then a custom change is required via the provided script for reversal of the sign of the coordinate axes. The adjusted axes should then be provide as input for calculating the motion between the sets of coordinate axes.
