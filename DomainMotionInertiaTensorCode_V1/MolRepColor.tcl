      
  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part shows the representation of the segments with different colors, styles etc
  ##  
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Aug 2015. Modified:April 20,2017 

  ##############################################################################################################

 
      set repno 1
      set dn 0
      foreach domain_name $all_domain_sels_name {
	  set str "_cartoon"
	  set domain_sel_cartoonrep [subst $$domain_name$str]

	  mol addrep $molid_model_full

	  set dom_colid [lindex $colorid_res $dn]

	  mol modselect $repno $molid_model_full $domain_sel_cartoonrep  
	  mol modcolor $repno $molid_model_full "ColorID" $dom_colid
	  mol modstyle $repno $molid_model_full $molrepstyle 0.5    
	  mol modmaterial $repno $molid_model_full $matrl
	  material change opacity $matrl $opqlevel 

	  set repno [expr $repno+1]
	  set dn [expr $dn+1]
      }
   
      ### instead of deleting the initial representation, 
      ### we can also turn it off
      mol delrep 0 $molid_model_full
      #mol showrep top 0 off

      display update on