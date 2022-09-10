  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part caputures snapshots from the VMD screen for different views
  ##
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Jan 2016. Modified:Oct 13,2016
  ##############################################################################################################
	    
  ###########################################################################################
  # Generate views from different perspective
  # You can define the specific views as necessary 
  ###########################################################################################

  if {$showLabelFig == 1} { 
      LabelDisplay $pdbf $fixedmol_cen "y" $datasetlabeldisplayoffset $txtlblcol
  }
  puts "\nGet display viewpoint center :\n"

  if {$nfl==0} {
	
      global viewpoints

      set viewpoints(0) [molinfo top get rotate_matrix] 
      #set viewpoints(1) [molinfo top get center_matrix]
      set viewpoints(1) {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}
      set viewpoints(2) [molinfo top get scale_matrix]
      set viewpoints(3) [molinfo top get global_matrix]	
  }
  set cp $viewpoints(1)
  puts $cp


  set filename_def snap_def.[format "%04d" $nfl].$resfileaffix
  set filename_y180 snap_y180.[format "%04d" $nfl].$resfileaffix
  set filename_y90 snap_y90.[format "%04d" $nfl].$resfileaffix
  set filename_ym90 snap_ym90.[format "%04d" $nfl].$resfileaffix
  set filename_z90 snap_z90.[format "%04d" $nfl].$resfileaffix
  set filename_xm90 snap_xm90.[format "%04d" $nfl].$resfileaffix
  set filename_x90 snap_x90.[format "%04d" $nfl].$resfileaffix
  set filename_xm90_ym90 snap_xm90_ym90.[format "%04d" $nfl].$resfileaffix

  if {$nfl==$initState_fileno} {
      ## change the snapshot filename for the initial model pdb so that it appears at the first for 
      ## the animation frames
      set filename_def snap_def.00.[format "%04d" $nfl].$resfileaffix
      set filename_y180 snap_y180.00.[format "%04d" $nfl].$resfileaffix
      set filename_y90 snap_y90.00.[format "%04d" $nfl].$resfileaffix
      set filename_ym90 snap_ym90.00.[format "%04d" $nfl].$resfileaffix
      set filename_z90 snap_z90.00.[format "%04d" $nfl].$resfileaffix
      set filename_xm90 snap_xm90.00.[format "%04d" $nfl].$resfileaffix
      set filename_x90 snap_x90.00.[format "%04d" $nfl].$resfileaffix
      set filename_xm90_ym90 snap_xm90_ym90.00.[format "%04d" $nfl].$resfileaffix

  }

  if {$use_reference_segmentation==0} {
      set num_rep $num_alldomains
  } else {
      set num_rep $num_alldomains
  }

  #### define the views you want to capture

  puts "Save display perspective: default"
  ResetDisplayScale $imscale $vmdwindow
  ShowMolRep $num_rep "off"	    
  render $rendercmd $resultsrootdir/tmp/$filename_def-1.$renderfigext
  ShowMolRep $num_rep "on"
  render $rendercmd $resultsrootdir/tmp/$filename_def-0.$renderfigext



  #set show_more_view 0

  if {$show_more_view==1} {
    puts "Save display perspective: rotate y +180"
    ResetDisplayScale $imscale $vmdwindow
    rotate y by 180	
    render $rendercmd $resultsrootdir/tmp/$filename_y180-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_y180-1.$renderfigext
    ShowMolRep $num_rep "on"



    puts "save display perspective: rotate y +90"
    ResetDisplayScale $imscale $vmdwindow
    rotate y by 90	
    render $rendercmd $resultsrootdir/tmp/$filename_y90-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_y90-1.$renderfigext
    ShowMolRep $num_rep "on"


    puts "Save display perspective: rotate y -90"
    ResetDisplayScale $imscale $vmdwindow	
    rotate y by -90
    render $rendercmd $resultsrootdir/tmp/$filename_ym90-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_ym90-1.$renderfigext
    ShowMolRep $num_rep "on"


    puts "Save display perspective: rotate z +90"
    ResetDisplayScale $imscale $vmdwindow
    rotate z by 90
    render $rendercmd $resultsrootdir/tmp/$filename_z90-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_z90-1.$renderfigext
    ShowMolRep $num_rep "on"


    puts "Save display perspective: rotate x +90"
    ResetDisplayScale $imscale $vmdwindow
    rotate x by 90
    render $rendercmd $resultsrootdir/tmp/$filename_x90-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_x90-1.$renderfigext
    ShowMolRep $num_rep "on"


    puts "Save display perspective: rotate x -90"
    ResetDisplayScale $imscale $vmdwindow
    rotate x by -90
    render $rendercmd $resultsrootdir/tmp/$filename_xm90-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_xm90-1.$renderfigext
    ShowMolRep $num_rep "on"


    puts "Save display perspective: rotate x -90, y -90"
    ResetDisplayScale $imscale $vmdwindow
    rotate x by -90
    rotate y by -90
    render $rendercmd $resultsrootdir/tmp/$filename_xm90_ym90-0.$renderfigext
    ShowMolRep $num_rep "off"	
    render $rendercmd $resultsrootdir/tmp/$filename_xm90_ym90-1.$renderfigext
    ShowMolRep $num_rep "on"
  }

  ResetDisplayScale $imscale $vmdwindow
