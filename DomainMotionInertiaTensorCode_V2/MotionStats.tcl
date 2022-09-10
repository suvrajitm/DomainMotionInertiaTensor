 
  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ## 
  ## This part summarizes the motion characterization results
  ##
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Nov 2015. Modified:April 20,2017 
  ## 
  ##
  ## ##############################################################################################################

  source ShowDomainAxisAngleQuat.tcl

  puts "\n\n\n\n*** Domain motion characterization ***"
  # two pairs of rotation motion for 3 input pdb models and so on ...

  puts "\nFixed Domain Molecule: $fixedDomain_name" 
 
  
  for {set dn 0} {$dn < $num_domains} {incr dn} {

      set domain_name [lindex $domain_sels_name $dn]
      set allmodelvar "AllModel_absorientQAxisAngle_$domain_name"
      set allmodel [subst $$allmodelvar]

      set initmodel_dom_cen [lindex $initmodel_domains_cen $dn]
      
      puts "\n\n*** Rotation quantification for $domain_name:\n"

      set rotaxislblcol $axislblcol
      set rotaxiscol "green"
     
      puts "allmodel, $allmodel"      

      ShowDomainAxisAngleQuat $nFiles $allmodel $domain_name $molid_model_full $allmodel_domains_cen $initmodel_dom_cen $rotaxiscol $rotaxislblcol $QuatRotAxisFile $QuatEulerAnglesFile $fixDomNum $dn
  }
 
  
  #mol showrep top $molid_model_full off

  #### TEMP
  #mol addrep top 
  #mol modselect 3 $molid_model_full "resid 1947 1948 1960 1961"  
  #mol modcolor 3 $molid_model_full "ColorID" 26 

  #mol addrep top 
  #mol modselect 4 $molid_model_full "resid 1848 1849 1896"  
  #mol modcolor 4 $molid_model_full "ColorID" 1 

  #mol addrep top 
  #mol modselect 5 $molid_model_full "resid 1913 1914 1918"  
  #mol modcolor 5 $molid_model_full "ColorID" 25 

  #mol addrep top 
  #mol modselect 6 $molid_model_full "resid 1836 1837"  
  #mol modcolor 6 $molid_model_full "ColorID" 3 
  



  ###########################################################################################
  ## Even after attempting to adjust the orientation of the axes for a domain when it is
  ## first calculated using the Orient script, we may not be able to fix the direction 
  ## for the same domain in the next model, due to slight variation.
  ## If one or more of the Principal axes gets flipped in direction from one model to the 
  ## next, then we can manually revert the sign of the flipped axes to calculate the rotation
  ## angle, although in the diplay, the axes will stay flipped
  ## Then we can call the absoluteOrientation for the re-oriented coordinate axes to 
  ## calculate the rotation
  ###########################################################################################

  ### Example 
  ### We can include this direction adjustment of the Principal axes of the two states
  ### in the script ComputeQuaternionAbsOrient.tcl 
  
  puts "\n\n***Example of manual adjustment of axes direction\n"
  
  #set domain_name "sel_ef2D1"

  set tm Tar_axes_nvh_all_model_$domain_name
  set Tarallmodel [subst $$tm]

  set stateA [lindex $Tarallmodel 0]
  set stateB [lindex $Tarallmodel 1]
             
  set PA1 [lindex $stateA 0]
  set PA2 [lindex $stateA 1]
  set PA3 [lindex $stateA 2]

             
  set PB1 [lindex $stateB 0]
  set PB2 [lindex $stateB 1]
  set PB3 [lindex $stateB 2]
   

  set dotPA1PB1 [vecdot $PA1 $PB1]
  set dotPA2PB2 [vecdot $PA2 $PB2]
  set dotPA3PB3 [vecdot $PA3 $PB3]



  if { $dotPA1PB1 < 0 } {
      set PB1 [vecscale -1 $PB1]
  }

  if { $dotPA2PB2 < 0 } {
      set PB2 [vecscale -1 $PB2]
  }
 
  if { $dotPA1PB1 < 0 } {
      set PB3 [vecscale -1 $PB3]
  }

  
  set PA1 [vecscale 1 $PA1]
  set PB1 [vecscale 1 $PB1]

  puts "pa1 = $PA1"
  set PA2 [vecscale 1 $PA2]
  set PB2 [vecscale 1 $PB2]

  set PA3 [vecscale 1 $PA3]
  set PB3 [vecscale 1 $PB3]


  set StateA_Axes [concat $PA1 $PA2 $PA3]
  set StateB_Axes [concat $PB1 $PB2 $PB3]

  puts "\n\nAdjusted Principal axes for state A--> stateB."
  set absorientQAxisAngle_dom  [AbsoluteOrientation $StateA_Axes $StateB_Axes]

  puts "\nCorrected Domain motion calculations for the specified domain $domain_name\n"
  

  ###########################################################################################
  ### Component of Rotation around a Specified axis: Quaternion Decomposition
  ###########################################################################################
  ### Example: Find the component of the quaternion 
  puts "\n\n*** Quaternion Decomposition ***"
  ### find the rotation component around a non-zero vector V 
  set V {1 0 0}


    
  set domain_name [lindex $domain_sels_name 1]
  set allmodeldom "AllModel_absorientQAxisAngle_$domain_name"
  puts "\nDomain : $domain_name"

  set Qrot [lindex $allmodeldom [list 1 0]]



  set domain_name [lindex $domain_sels_name 2]
  set allmodeldom "AllModel_absorientQAxisAngle_$domain_name"
  puts "\nDomain : $domain_name"

  set Qrot [lindex $allmodeldom [list 1 0]]
 

  ### test case 
  set Qrot {0.99 -0.09 -0.07 0.04}
  set Qrot [unitQuaternion $Qrot]

  puts "\nQ = $Qrot\n"
  puts "Specified Axis V = $V"

  #set qrotdecomp [QuatDecomposition $Qrot $V]
  set qrotdecomp [QuatDecompositionClifford $Qrot $V]

  set Qt [lindex $qrotdecomp 0]
  set qtrotaxisangle [lindex $qrotdecomp 1]

  set Qtaxis [lindex $qtrotaxisangle 0]
  set Qtangle [lindex $qtrotaxisangle 1]

  set Qs [lindex $qrotdecomp 2]
  set qsrotaxisangle [lindex $qrotdecomp 3]

  set Qsaxis [lindex $qsrotaxisangle 0]
  set Qsangle [lindex $qsrotaxisangle 1]


  puts "\n\nComponent of Q around V is Qt : $Qt"
  puts "Corresponding axis of rotation : $Qtaxis"
  puts "Corresponding angle of rotation : $Qtangle\n"
