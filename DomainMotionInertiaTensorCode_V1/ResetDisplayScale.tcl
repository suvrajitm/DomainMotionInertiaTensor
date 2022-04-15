##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
##
## This part restores the VMD screen views as defined by the viewpoints
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Aug 26, 2016. Modified:Aug 26, 2016
##############################################################################################################

proc ResetDisplayScale {imscale vmdwindow} {
    display resetview
    display resize $vmdwindow $vmdwindow
    global viewpoints
    display update
    scale to $imscale
    molinfo top set rotate_matrix $viewpoints(0)
    molinfo top set center_matrix $viewpoints(1)
    molinfo top set scale_matrix $viewpoints(2)
    molinfo top set global_matrix $viewpoints(3)
}