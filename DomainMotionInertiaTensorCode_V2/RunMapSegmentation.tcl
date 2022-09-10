## Tcl wrapper script for executing the segmentation program
## Execute_RunWatershedSegmentMap.sh is the shell script which executes the matlab standalone function for segmentation
## 
proc RunMapSegmentation {mcr_root mapFile Gauss_FilterSize walgo pixSize merge_level merge_size extra_merging} {

  exec ./Execute_RunWatershedSegmentMap.sh $mcr_root $mapFile $Gauss_FilterSize $walgo $pixSize $merge_level $merge_size $extra_merging
  #exec ./Execute_RunWatershedSegmentMap.sh "/guam.raid.cluster.software/matlab2013/" $mapFile $Gauss_FilterSize $walgo $pixSize $merge_depth $extra_merging

}