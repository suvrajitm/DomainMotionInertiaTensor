  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part computes the unit quaternion and rotation axis and rotation angle for each domain going from 
  ## state A to state B using the abosolute orientation method implemented in the script "AbsoluteOrientation.tcl"
  ## 
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Dec 07, 2015. Modified:April 20,2017 
  ##############################################################################################################
  ## Customize this script to compute the coordination axes transformation between the two principal axes
  ## of a particular domain in two states, by calling the function "AbsoluteOrientation" 
  ## from within ComputeQuaternionAbsOrient.tcl
  ##############################################################################################################


  ###########################################################################################
  ## For any general "Domain", here is a template to call the ComputeQuaternionAbsOrient.tcl
  ###########################################################################################
  #for {set nfl 0} {$nfl < $nFiles} {incr nfl} {
  #    set absorientQAxisAngle_Domain [ComputeQuaternionAbsOrient $Domain_axes_nvh_all_model $nfl $initState_fileno]
  #    lappend AllModel_absorientQAxisAngle_Domain $absorientQAxisAngle_Domain   	
  #}



  # Calculating the SSU Head and Body movements and also the domains of the ligand eEF2 domain III , domain IV w.r.t domain I
      
  puts "\n\n***Computing the independent domain rotations with a closed form solution of Axes transformation using unit quaternions...\n"   

  ###########################################################################################
  
  # Extract the angles using Quaternion method
  # A different implementation on VMD 1.9.3 with TCL8.6 for the tclOO package
  source AbsoluteOrientation.tcl
  source  ComputeQuaternionAbsOrient.tcl 

  # puts "\nAxes optimal rotation using unit-quaternion based on Horn's quaternion method ...\n"
  # match the corresponding axes of reference and target which should be paired together for transformation    
  
  # reset the variables that stores the axes for each domain when starting a new run

  for {set dn 0} {$dn < $num_domains} {incr dn} {
  
    set domain_name [lindex $domain_sels_name $dn]

    set allmodel "AllModel_absorientQAxisAngle_$domain_name"
    set $allmodel {}
  }
  
  for {set dn 0} {$dn < $num_domains} {incr dn} {

       set domain_name [lindex $domain_sels_name $dn]

	# we can calculate the orientation axis-angle for individual domains here 
	if {$quaternion_calc > 0} {

	    puts "\n\nRotation axis using unit quaternions\n"

	    for {set nfl 0} {$nfl < $nFiles} {incr nfl} {
	        

		puts "\n\nRotation between Reference(or unrotated)-model: $initState_fileno - $nfl"

		set tvar "Tar_axes_nvh_all_model_$domain_name"
		set allmodeldomain [subst $$tvar]

		puts "\n\n*** Fixed Reference, Computing rotation motion for $domain_name...\n"
		set absorientQAxisAngle_Domain [ComputeQuaternionAbsOrient $allmodeldomain $nfl $initState_fileno]
		set allmodel "AllModel_absorientQAxisAngle_$domain_name"
		lappend $allmodel $absorientQAxisAngle_Domain    
			
	    }
	}
  }

 
