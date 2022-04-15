  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part plots the graphics for rotation axis of a domain going from state A to state B 
  ## corresponding to the computed rotation unit quaternion from the absolute orientation method
  ## 
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Feb, 2015. Modified:Feb 28,2016 
  ##############################################################################################################

  proc ShowDomainAxisAngleQuat {nFiles AllModel_absorientQAxisAngle domname mol dom_cen rotaxiscol rotaxislblcol QuatRotAxisFile QuatEulerAnglesFile} {
      
      set showrotaxislbl 1

      for {set nfl 0} {$nfl < [expr $nFiles]} {incr nfl} {

	  if {$showrotaxislbl} {
	      #set axislabel "e($domname)$nfl"
	      set axislabel "e$nfl"

	  } else {
	      set axislabel ""
	  }

	  set quat_indx [list $nfl 0]
	  set axis_indx [list $nfl 1]
	  set angle_indx [list $nfl 2]
	  set eulerangles_indx [list $nfl 3] 

	  set qdomain [lindex $AllModel_absorientQAxisAngle $quat_indx]
	  set axisdomain [lindex $AllModel_absorientQAxisAngle $axis_indx]
	  set angledomain [lindex $AllModel_absorientQAxisAngle $angle_indx]
	  set euleranglesdomain [lindex $AllModel_absorientQAxisAngle $eulerangles_indx]

	  puts "Theta = $angledomain,\tAxis = $axisdomain,\tQ = $qdomain"
          # set the scale factors for drawing the axis vector as we see fit, it's kind of arbitrary right now 
	  # the scale I have used is proportional to the sqroot/cuberoot of the angle values 
          ### angles are in degrees

	  if {$angledomain < 0} {

	     set scale1 1.0
	     set scale2 1.01
  
	  } else {  

	     set scale1 [expr 20*(($angledomain) ** (1.0 / 2))]
	     set scale2 [expr 1.01*$scale1]
	  }


	  set rotOrigin $dom_cen
	  puts $scale1
	  puts $scale2
	  puts "Axis domain: $axisdomain $angledomain"
	  #graphics $mol color "green"
          graphics $mol color $rotaxiscol
	  vmd_draw_vector $mol $rotOrigin [vecscale $scale1 $axisdomain]

	  graphics $mol color $rotaxislblcol
	  graphics $mol text [vecadd $rotOrigin [vecscale $scale2 $axisdomain]] $axislabel size 1.5 thickness 2
      }
      
      # 
      puts $QuatRotAxisFile $domname
      puts $QuatEulerAnglesFile $domname
      puts "\n\nRotation Angle, Rotation Axis:\n"
      for {set nfl 0} {$nfl < [expr $nFiles]} {incr nfl} {

	  set quat_indx [list $nfl 0]
	  set axis_indx [list $nfl 1]
	  set angle_indx [list $nfl 2]
	  set eulerangles_indx [list $nfl 3] 

	  set qdomain [lindex $AllModel_absorientQAxisAngle $quat_indx]
	  set axisdomain [lindex $AllModel_absorientQAxisAngle $axis_indx]
	  set angledomain [lindex $AllModel_absorientQAxisAngle $angle_indx]
	  set euleranglesdomain [lindex $AllModel_absorientQAxisAngle $eulerangles_indx]

	  puts "Theta = $angledomain \tAxis = $axisdomain \tQ = $qdomain"
	  puts $QuatRotAxisFile "Theta = $angledomain \tAxis = $axisdomain \tQ = $qdomain"
	  puts $QuatEulerAnglesFile "EulerAngles = $euleranglesdomain \tQ = $qdomain"
      }
      puts "\n"
      puts $QuatRotAxisFile "\n\n"
      puts $QuatEulerAnglesFile "\n\n"
  }