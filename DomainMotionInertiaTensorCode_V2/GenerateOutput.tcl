  
  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ## 
  ## This part generates the ouput files
  ##
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Nov 2015. Modified:April 20,2017 
  ## 
  ##
  ##############################################################################################################

 ###########################################################################################

  if {$map_segment_only==0} {
      ######## save vmd state
      ######## don't use VMD's default save_state, it does not save the info to restore the 
      ######## visualization state in the correct coordiante frame
      #### save_state $resultsdir/resvisstate.vmd 


      #### use this extension instead 
      source Utilities/VMDextensions.tcl
      saveFullState $resultsdir/resVisState.vmd 



  }
  ###########################################################################################

  if {$map_segment_only==0 && $render=="ray"} {

    puts "Converting the rendered images into common format."
    puts "Running tachyon with AO lighting"
    # convert tachyon rendered images to bmp
    set rayimages [glob $resultsrootdir/tmp/*.$renderfigext]
    set lenviews [llength $rayimages]

    for {set imv 0} {$imv < $lenviews-1 } {incr imv 2} {
	puts "Image view: [expr $imv/2]"
	exec RenderVMD.sh [lindex $rayimages $imv] 1
	exec RenderVMD.sh [lindex $rayimages [expr $imv+1]] 2
    }
  }
  


  ###########################################################################################
  # Make animation gif of the individual snapshots
  ########################################################################################### 

  # del raw image file after rendering and conversion is done
  # if delrawimgfile (defined above when creating results directory) is set to 1
 
  if {$map_segment_only==0 && ($nFiles >= 1 || $nFrame >= 1)} {

      #convert to a gif movie and optionally delete the .rgb/.bmp files
      puts "\nCreating the gif files from raw/converted image files...\n"
      if {$show_more_view==1} {
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_xm90_ym90.*-0*.$figext $resultsrootdir/tmp/xm90_ym90-0-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_xm90.*-0*.$figext $resultsrootdir/tmp/xm90-0-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_x90.*-0*.$figext $resultsrootdir/tmp/x90-0-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_y180.*-0*.$figext $resultsrootdir/tmp/y180-0-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_y90.*-0*.$figext $resultsrootdir/tmp/y90-0-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_ym90.*-0*.$figext $resultsrootdir/tmp/ym90-0-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_z90.*-0*.$figext $resultsrootdir/tmp/z90-0-$resfileaffix.gif

	# axes only
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_xm90_ym90.*-1*.$figext $resultsrootdir/tmp/xm90_ym90-1-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_xm90.*-1*.$figext $resultsrootdir/tmp/xm90-1-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_x90.*-1*.$figext $resultsrootdir/tmp/x90-1-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_y180.*-1*.$figext $resultsrootdir/tmp/y180-1-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_y90.*-1*.$figext $resultsrootdir/tmp/y90-1-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_ym90.*-1*.$figext $resultsrootdir/tmp/ym90-1-$resfileaffix.gif
	exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_z90.*-1*.$figext $resultsrootdir/tmp/z90-1-$resfileaffix.gif

      }
      exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_def.*-0*.$figext $resultsrootdir/tmp/def-0-$resfileaffix.gif
      exec convert -delay 100 -loop 4 -quality 100 $resultsrootdir/tmp/snap_def.*-1*.$figext $resultsrootdir/tmp/def-1-$resfileaffix.gif

      if {$delrawimgfile==1} {
      # optionally delete the generated image files 
          puts "\nDeleting raw and converted image files..."
	  set imgfiles [glob -nocomplain $resultsrootdir/tmp/*.$figext]
	  set rawimgfiles [glob -nocomplain $resultsrootdir/tmp/*.$renderfigext]

	  set nimgfiles [llength $imgfiles]
	  set nrawimgfiles [llength $rawimgfiles]

	  for {set f 0} {$f < $nimgfiles} {incr f} {
	    file delete [lindex $imgfiles $f]
	    file delete [lindex $rawimgfiles $f]
	  }

      }
  }

# write all vars to file
  foreach var [info vars] {
      catch {
        puts $settingsfile "$var:"
      	puts $settingsfile [subst $$var]
      	puts $settingsfile ""
      	
      }
  }
  ###########################################################################################
  # close all open (*.dat) files
  ###########################################################################################
  if {$map_segment_only==0} {

    puts "\nClosing open .dat files ..."
    close $TensorAnglesFile
    close $PrAxesFile 
    close $QuatEulerAnglesFile 
    close $QuatRotAxisFile
    close $settingsfile

  }
 
  ########################################################################################### 
  # Move all the results file to the specified results directory $resultsdir
  ########################################################################################### 
  if {$map_segment_only==0} {
      puts "\nMoving all generated files to the Results folders..."
      if {$delrawimgfile==0} {
	  file rename {*}[glob $resultsrootdir/tmp/*$resfileaffix*.$figext] $resultsdir/images
	  if {$render=="ray"} {
	      file rename {*}[glob $resultsrootdir/tmp/*$resfileaffix*.$renderfigext] $resultsdir/rawimages
	  }
	}
      file rename {*}[glob $resultsrootdir/tmp/*$resfileaffix*.gif] $resultsdir/gifmovies
      file rename {*}[glob $resultsrootdir/tmp/*$resfileaffix*.dat] $resultsdir/datfiles
      if {$segment_map==1 } {
	file rename Maps $resultsdir/Maps
      }
      puts "\nAll result files (*.dat *.$figext *.renderfigext *.gif) generated in the working directory were moved to the directory $resultsdir ...\n\n\n"
 
  }
