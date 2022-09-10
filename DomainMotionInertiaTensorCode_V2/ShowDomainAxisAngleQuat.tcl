  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part plots the graphics for rotation axis of a domain going from state A to state B 
  ## corresponding to the computed rotation unit quaternion from the absolute orientation method
  ## Added in May/June 2022: Determine the hinge axis also known as screw axis in the literature for 
  ## 3D rigid body motion 
  ## 
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Feb, 2015. Modified: August,2022 
  ###############################################################################################################

  proc ShowDomainAxisAngleQuat {nFiles AllModel_absorientQAxisAngle domname molid allmodel_domains_cen initmodel_dom_cen rotaxiscol rotaxislblcol QuatRotAxisFile QuatEulerAnglesFile fixDomNum dn} {
      
       
       proc lmap {lambda m rlist} {
    	  set result {}
    	  foreach item $rlist {
        	lappend result [apply $lambda $m $item]
    	  }
    	  return $result
       }
      
       proc plotRotAxisPoint {molid rotaxiscol axislabel scale1 scale2 rotaxisOrigin rotaxisDomain contactDist} {
  	  puts "\n\nHinge axis: $rotaxisDomain" 
	  puts "Hinge point: $rotaxisOrigin"
	  graphics $molid color $rotaxiscol
  	  vmd_draw_vector $molid $rotaxisOrigin [vecscale $scale1 $rotaxisDomain]
  	  #graphics $molid text [vecadd $rotaxisOrigin [vecscale $scale2 $rotaxisDomain]] $axislabel size 1.5 thickness 2	
	  set minmax_coord [measure minmax [atomselect top all]]
	  lassign $minmax_coord minc maxc 
          set bs [vecdist $minc $maxc]
          set bs2 [expr 0.8*$bs/2.0]     
 	  set st_pt [vecadd $rotaxisOrigin [vecscale [expr -1.0*$bs2] $rotaxisDomain]]
	  set end_pt [vecadd $rotaxisOrigin [vecscale $bs2 $rotaxisDomain]]
	  graphics $molid line $st_pt $end_pt style dashed width 7
		
          set hingeaxisline [vecsub $end_pt $st_pt]

	  set n [expr round(2.0*$bs2/3.0)]
	  set h_selres {}
	  puts "\nFinding atoms close to the point on the axis ..."
         
          for {set j 0} {$j<$n} {incr j 1} {
              #puts "j=$j, n=$n"
	      set b [expr (-1.0 + (2.0*$j)/($n-1))*$bs2]
	      set hp [vecadd $rotaxisOrigin [vecscale $b $rotaxisDomain]]
	      #puts "hp=$hp"
 	      lassign $hp hx hy hz
	      set mol_selres [atomselect $molid "sqr(x-$hx)+sqr(y-$hy)+sqr(z-$hz)< sqr($contactDist)"]
	      lappend h_selind [$mol_selres get index]
	      graphics $molid sphere $hp radius 1.0
	  }
	
	  set hinge_selinds [lsort -unique -integer [concat {*}$h_selind]]  
	  set tmp_residuenums [atomselect $molid "same residue as index $hinge_selinds"]
	  set hinge_residuenums [lsort -unique [$tmp_residuenums get residue]]
	  
	  set lambda_resname {{x y} {lsort -unique [[atomselect $x "residue $y"] get resname]}}
	  set hinge_resnames [lmap $lambda_resname $molid $hinge_residuenums]
	  
	  set lambda_resid {{x y} {lsort -unique [[atomselect $x "residue $y"] get resid]}}
	  set hinge_resids [lmap $lambda_resid $molid $hinge_residuenums]
	  

	  puts "\nHinge residues:$hinge_residuenums"
	  mol addrep $molid 
          mol modselect 2 $molid "residue $hinge_residuenums"  
          mol modcolor 2 $molid "ColorID" 1
	  mol modstyle 2 $molid "VDW"
	  return [list $hinge_residuenums $hinge_resids $hinge_resnames]

      }
      
      set showrotaxislbl 1
      puts $QuatRotAxisFile $domname
      puts $QuatEulerAnglesFile $domname

      for {set nfl 0} {$nfl < [expr $nFiles]} {incr nfl} {
	
	
    	
          set dom_cen [lindex [lindex $allmodel_domains_cen $nfl] $dn]
            
	        
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

	  set qDomain [lindex $AllModel_absorientQAxisAngle $quat_indx]
	  set rotaxisDomain [lindex $AllModel_absorientQAxisAngle $axis_indx]
	  set rotangleDomain [lindex $AllModel_absorientQAxisAngle $angle_indx]
	  set euleranglesDomain [lindex $AllModel_absorientQAxisAngle $eulerangles_indx]

  	  set qDomain [lindex $AllModel_absorientQAxisAngle $quat_indx]
	  

          # added May 2022
	  global Pi eps contactDist
	  
	  # theta angle should be converted to radian, since rotangleDomain was converted to degrees
	  set theta [expr $rotangleDomain*$Pi/180]  
	  puts "l-theta, $theta"		
	  set lhat $rotaxisDomain
	  
	  puts "Theta = $rotangleDomain,\tAxis = $rotaxisDomain,\tQ = $qDomain"
          # set the scale factors for drawing the axis vector as we see fit, it's kind of arbitrary right now 
	  # the scale I have used is proportional to the sqroot/cuberoot of the angle values 
          ### angles are in degrees

	  if {$theta < 0.0} {
	     set scale1 1.0
	     set scale2 1.01
  
	  } else {  
	     set scale1 [expr 30.0*(($rotangleDomain) ** (1.0 / 2.0))]
	     set scale2 [expr 1.01*$scale1]
	  }
	
	  #set rotaxisOrigin $dom_cen
	  #puts "q-rotaxis plot through comB"
	  #plotRotAxisPoint $molid $rotaxiscol $axislabel $scale1 $scale2 $rotaxisOrigin $rotaxisDomain $contactDist
	  
	  
	  #### 
	  #center of mass of the domain in state A to state B
          set comA $initmodel_dom_cen
	  set comB $dom_cen
	  
	  # compute the COM of the moving domain in the two states A --> B, so B - A
          set t [vecsub $comB $comA] 
	  #set d [vecdist $comB $comA]
	  set d [veclength $t]

	  puts "com A:$comA,  comB:$comB"
	  puts "com-movement of domain $domname: $d."
			
	  # both implementations 1 and 2 have equivalent/exact expressions for the position of the screw axis 
	  # but expressed in two forms   	  
 	  set findScrewAxisOrigPoint1 1
 	  if {$findScrewAxisOrigPoint1 > 0 && $rotangleDomain > 0.0} {
		puts "\nfindScrewAxisOrigPoint1"
		if {$dn==$fixDomNum} {
	      	     break
	        }
		
		set cotthetaby2 [expr 1.0/(tan($theta/2.0))]
		set tc [vecscale 0.5 [vecscale $cotthetaby2 [veccross $lhat $t]]]
		set midAB [vecscale 0.5 [vecadd $comA $comB]]
	        set d1 [vecdot $lhat $midAB]
  	        set m1 [vecscale $d1 $lhat]
		
                set Cperp [vecadd [vecsub $midAB $m1] $tc]

		set rotaxisOrigin $Cperp
		#re-assign rotaxisDomain
		set rotaxisDomain $lhat 
			
	  	set rotaxiscol "green"
	  	set hinge_residues [plotRotAxisPoint $molid $rotaxiscol $axislabel $scale1 $scale2 $rotaxisOrigin $rotaxisDomain $contactDist]
	 	set hinge_residuenums [lindex $hinge_residues 0]
	 	set hinge_resids [lindex $hinge_residues 1]
	  	set hinge_resnames [lindex $hinge_residues 2]
	  	puts $QuatRotAxisFile "Hinge axis point = $rotaxisOrigin"
	  	puts $QuatRotAxisFile "Hinge residue numbers: $hinge_residuenums"
	  	puts $QuatRotAxisFile "Hinge resids: $hinge_resids"
		puts $QuatRotAxisFile "Hinge residue names: $hinge_resnames"
	  }
 	  
	  set findScrewAxisOrigPoint2 0
 	  if {$findScrewAxisOrigPoint2 > 0 && $rotangleDomain > 0.0} {
		puts "\nfindScrewAxisOrigPoint2"
		if {$dn==$fixDomNum} {
	      	     break
	        }
		
		set tpar [vecscale $lhat [vecdot $lhat $t]] 
		set tperp [vecsub $t $tpar]
		set c [veccross $lhat $tperp]
		set ccot [vecscale [expr 1.0/(tan($theta/2.0))] $c]
		puts "ccot:$ccot"
		set m0 [veccross $lhat [veccross $comA $lhat]]
		set C [vecadd $m0 [vecscale 0.5 [vecadd $tperp $ccot]]]	
	
		set rotaxisOrigin $C
		#re-assign rotaxisDomain
		set rotaxisDomain $lhat 
			
	  	set rotaxiscol "blue"
	  	set hinge_residues [plotRotAxisPoint $molid $rotaxiscol $axislabel $scale1 $scale2 $rotaxisOrigin $rotaxisDomain $contactDist]
		
		set hinge_residuenums [lindex $hinge_residues 0]
	 	set hinge_resids [lindex $hinge_residues 1]
	  	set hinge_resnames [lindex $hinge_residues 2]
	  	puts $QuatRotAxisFile "\nHinge axis point = $rotaxisOrigin"
	  	puts $QuatRotAxisFile "Hinge residue numbers: $hinge_residuenums"
	  	puts $QuatRotAxisFile "Hinge resids: $hinge_resids"
		puts $QuatRotAxisFile "Hinge residue names: $hinge_resnames"
		
		
	   }


	  
 			 
      }
      
      # 
      puts "\n\nRotation Angle, Rotation Axis:\n"
      for {set nfl 0} {$nfl < [expr $nFiles]} {incr nfl} {

	  set quat_indx [list $nfl 0]
	  set axis_indx [list $nfl 1]
	  set angle_indx [list $nfl 2]
	  set eulerangles_indx [list $nfl 3] 

	  set qDomain [lindex $AllModel_absorientQAxisAngle $quat_indx]
	  set rotaxisDomain [lindex $AllModel_absorientQAxisAngle $axis_indx]
	  set rotangleDomain [lindex $AllModel_absorientQAxisAngle $angle_indx]
	  set euleranglesdomain [lindex $AllModel_absorientQAxisAngle $eulerangles_indx]

	  puts "Theta = $rotangleDomain \tAxis = $rotaxisDomain \tQ = $qDomain"
	  puts $QuatRotAxisFile "Theta = $rotangleDomain \tAxis = $rotaxisDomain \tQ = $qDomain"
	  puts $QuatEulerAnglesFile "EulerAngles = $euleranglesdomain \tQ = $qDomain"

      }
      puts "\n"
      puts $QuatRotAxisFile "\n\n"
      puts $QuatEulerAnglesFile "\n\n"

    
  }
