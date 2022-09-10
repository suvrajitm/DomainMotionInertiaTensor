##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
##
## This part uses EMAN 1/EMAN 2 function pdb2mrc/e2pdbmrc to create volumetric map from pdb
## 
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Aug 2015. Modified:Oct 30,2015
##############################################################################################################

proc MapSegment3d {mcr_root pdbFile mapFile mapBoxSize mapRes Gauss_FilterSize walgo pixSize merge_level merge_size extra_merging} {
  
  # it takes huge amount of time if a low resolution map >= 10.0 A is generated
  puts "\nGenerating map ($mapBoxSize x $mapBoxSize x $mapBoxSize) with a resolution = $mapRes A and pixel size = $pixSize A from the input .pdb file ..."

  exec pdb2mrc $pdbFile $mapFile box=$mapBoxSize apix=$pixSize res=$mapRes center het 
  #exec e2pdb2mrc.py $pdbFile $mapFile --box 240 --apix $pixSize --res $mapRes --het --verbose 2


  # This is a Watershed segmentation using Matlab/DipImage Library(c++ mex file) executable.To Do : platform independent python/c/c++ code for the watershed segmentation 
  puts "\n3D Watershed Segmentation of the map is in progress ...\n"
  puts "\nGaussian smoothing filter size: ($Gauss_FilterSize,$Gauss_FilterSize,$Gauss_FilterSize) ...\n"

  set Res_filter [expr $Gauss_FilterSize*$pixSize/(2.0*0.53)]
  set mapRes_effective [expr sqrt($mapRes*$mapRes + $Res_filter*$Res_filter)]
  puts "Effective resolution of the generated map after application of the Gaussian smoothing filter: [format "%.2f" $mapRes_effective A]\n"

  puts "Region-merging with max-merging-size = $merge_size, merging-level = $merge_level ...\n"
  puts "\nPlease wait for several seconds to verify the segmentation (scroll through the segmented slices). To proceed, the 'Slicer' window needs to be closed.\n"
  source RunMapSegmentation.tcl
  set v [RunMapSegmentation $mcr_root $mapFile $Gauss_FilterSize $walgo $pixSize $merge_level $merge_size $extra_merging]
  
  puts "\nSegmentation of the map is done ...\n\n"
}
