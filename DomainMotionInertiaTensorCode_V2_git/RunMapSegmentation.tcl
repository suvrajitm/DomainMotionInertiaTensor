##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
##
## Tcl wrapper script for executing the segmentation program
## 
## Execute_RunWatershedSegmentMap.sh is the shell script 
## which executes the matlab standalone function for segmentation
## 
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Aug 26, 2016. Modified:Aug 26, 2016
##############################################################################################################

proc RunMapSegmentation {mcr_root mapFile Gauss_FilterSize walgo pixSize merge_level merge_size extra_merging} {

  exec ./Execute_RunWatershedSegmentMap.sh $mcr_root $mapFile $Gauss_FilterSize $walgo $pixSize $merge_level $merge_size $extra_merging
  #exec ./Execute_RunWatershedSegmentMap.sh "/guam.raid.cluster.software/matlab2013/" $mapFile $Gauss_FilterSize $walgo $pixSize $merge_depth $extra_merging

}