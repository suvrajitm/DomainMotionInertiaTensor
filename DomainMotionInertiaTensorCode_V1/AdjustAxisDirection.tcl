  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ## 
  ## This part is used to adjust the direction of the principal axis, which might get flipped from one 
  ## to another model  
  ## 
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: May 2015. Modified:April 22,2017 
  ##############################################################################################################

proc AdjustAxisDirection {P_ref P} {

    set a1_ref [lindex $P_Ref 0]
    set a2_ref [lindex $P_Ref 1]
    set a3_ref [lindex $P_ref 2]

    set a1 [lindex $P 0]
    set a2 [lindex $P 1]
    set a3 [lindex $P 2]


    ### Ideally we do not need to do this extra steps for alignment, but in this case we have a problem of 
    ### flipping axis direction for the same domain due to slight change of mass distribution and 
    ### we have to enforce directionality constraint

    if { [vecdot $a1_ref $a1] < 0 } {

	 set a1 [vecinvert $a1]

    }

    if { [vecdot $a2_ref $a2] < 0 } {

         set a2 [vecinvert $a2]
    }

    if { [vecdot $a3_ref $a3] < 0 } {

	 set a3 [vecinvert $a3]
    }

    set P [list $a1 $a2 $a3]
  
	   
    return $P
	
}