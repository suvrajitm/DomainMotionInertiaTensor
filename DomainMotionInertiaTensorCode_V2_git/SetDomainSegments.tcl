##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
## 
## This is where you can set the domain definitions obtained either through segmentation or from a reference
## segmentation 
## 
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: May 2015. Modified:Oct 23,2017 
##############################################################################################################



# set the segments for target and ligand molecules properly, 
# We can either set the residues for each domain manually or through segmentation approach
#   # set sel_lsu "(segname 25S 5S 5D8S \"L.*\") and (not segname P1 P2 P5 L45 SO1 GDP) and (name P)" 
#   # set sel_lsu_cartoon "(segname 25S 5S 5D8S "\L.*\") and (not segname P1 P2 P5 L45 SO1 GDP)"
#   # set sel_ssu "(segname 18S "\S.*\" RACK) and (not segname P1 P2 P5 L45 SO1 GDP) and (name P)"
#   # set sel_ssu_cartoon "(segname 18S  "\S.*\" RACK) and (not segname P1 P2 P5 L45 SO1 GDP)"
#   # set sel_ef2 "(segname P1 P2 P5 L45 SO1 GDP)"
#   #	set sel_ef2_cartoon "(segname P1 P2 P5 L45 SO1 GDP)"

puts "\nSetting domain segments with the corresponding residue indices ...\n" 


set ef2data 0
#set ef2_full_ligand 0
set mData0 0
set mData 1

set all_atom_tensorcalc 0


if {$all_atom_tensorcalc==0} {
set atoms_to_use "P OP1 O1P OP2 O2P OP3 O3P O3' O5' C3' C4' C5'" 
}

# otherwise use only the P atoms for faster inertia tensor calculations


if {$segment_map==1 && $use_reference_segmentation==0} {

  puts "\nUsing volume segmented definitions for the domains.\n"
  
  set all_domain_sels_name [list "sel_lsu" "sel_ssu_body" "sel_ssu_head" "sel_ef2" "sel_ef2D1" "sel_ef2D2" "sel_ef2D3" "sel_ef2D4" "sel_ef2D5" "sel_ef2D1D2" "sel_ef2D3D5"]

  set domain_sels_from_volSegment [list "sel_lsu" "sel_ssu_body" "sel_ssu_head" "sel_ef2"]

  ### find the domains which are not determined by the volume segments but has to be defined from reference
  set domain_sels_refdef {}
  foreach all_dom_sel $all_domain_sels_name {
      if {[lsearch -exact $domain_sels_from_volSegment $all_dom_sel]==-1} {
          lappend domain_sels_refdef $all_dom_sel
      }
  }

  # select particular ligand domain to fix and also define the residue selection for each subdomain properly   
  # set sel_ef2 "(segname P1 P2 P5 L45 SO1 GDP)"
  # set sel_ef2 "(segname P1 P2 P3 P5 L45 SO1 GDP)"
  #set sel_ef2G "(segname P2 and resid 219 to 328)"
  set sel_ef2D1 "((segname P1) or (segname P2 and ((resid 67 to 218) or (resid 329 to 345))))"
  set sel_ef2D2 "(segname P2 and (resid 346 to 481))"
  set sel_ef2D3 "(segname P2 and (resid 482 to 558))"
  set sel_ef2D4 "(segname P2 and (resid 559 to 724) or (segname L45 and resid 801 to 842))"
  set sel_ef2D5 "((segname P5) or (segname L45) and (resid 763 to 800))"
  set sel_ef2D1D2 "((segname P1) or (segname P2 and ((resid 67 to 218) or (resid 329 to 345)))) or (segname P2 and (resid 346 to 481))"
  set sel_ef2D3D5 "(segname P2 and (resid 482 to 558)) or ((segname P5) or (segname L45) and (resid 763 to 800))"

  if {$num_segments != [llength $domain_sels_from_volSegment]} {
    #puts "Warning: The number of map segmented domains and number of defined domains do not match !"
    error "The number of map segmented domains and number of defined domains do not match !"
  }

  set domain_sels_from_volSegment_reshuffle {}

  for {set k 0} {$k < $num_segments} {incr k} {
    lappend domain_sels_from_volSegment_reshuffle [lindex $domain_sels_from_volSegment [lindex $mapSegSort_idx $k]]
  }


  if { $num_segments > 0 } {
    set segmented_domains_defined 0
  }

  for {set seg_no 0} {$seg_no < [llength $domain_sels_from_volSegment_reshuffle]} {incr seg_no} {

      set domain_seg_name [lindex $domain_sels_from_volSegment_reshuffle $seg_no]
      #set domain_sel [subst $$domain_name]
      if {([info exists $domain_seg_name]==1)} {
          unset $domain_seg_name
      }
      ### set only the variables that have not been set before 
      if {[lsearch -exact $domain_sels_refdef $domain_seg_name]==-1} {
          puts "\nSetting segment number $seg_no, $domain_seg_name definition.\n"
          #### assign the segmented volumes to the corresponding segment names
          set domResvar "residue_idx_vol$seg_no"
          set domRes_ids [subst $$domResvar] 
          set $domain_seg_name "(residue $domRes_ids)"
          set segmented_domains_defined 1
      }
  }
  
  if { $segmented_domains_defined == 0} {
    error "Segmented domains were not defined/assigned with the corresponding residue indices!"
  }
  
  

}



if {$use_reference_segmentation==1} { 

  puts "\nUsing reference definitions for the domains.\n"
  if {$ef2data==1} {

      #set all_domain_sels_name [list "sel_lsu" "sel_ssu_body" "sel_ssu_head" "sel_ef2" "sel_ef2G" "sel_ef2D1" "sel_ef2D2" "sel_ef2D3" "sel_ef2D4" "sel_ef2D5" "sel_ef2D1D2" "sel_ef2D3D5"]
      set all_domain_sels_name [list "sel_lsu" "sel_ssu_body" "sel_ssu_head" "sel_ef2" "sel_ef2D1" "sel_ef2D2" "sel_ef2D3" "sel_ef2D4" "sel_ef2D5" "sel_ef2D1D2" "sel_ef2D3D5"]

      #set all_domain_sels_name [list "sel_lsu" "sel_ssu_body" "sel_ssu_head" "sel_ef2"]

      #set sel_lsu "(segname 25S 5S 5D8S \"L.*\") and (not segname P1 P2 P5 L45 SO1 GDP)" 
      #set sel_ssu "(segname 18S \"S.*\" RACK) and (not segname P1 P2 P5 L45 SO1 GDP)"
      #set sel_ef2 "(segname P1 P2 P5 L45 SO1 GDP)"
      set sel_lsu  "(segname 25S 5S 5D8S \"L.*\") and (not segname L45)" 
      set sel_ssu "(segname 18S \"S.*\" RACK) and (not segname SO1)"
      set sel_ssu_body "((segname 18S and not (resid 1150 to 1635)) or (segname S11 S12 S15 S17 S1E S2 S21E S24E S26E S27E S30E S4 S4E S5 S6E S7E S8 S8E))"
      set sel_ssu_head "((segname 18S and (resid 1150 to 1635)) or (segname RACK S10 S10E S12E S13 S14 S17E S19 S19E S25E S28E S3 S7 S9))"
      #set sel_ef2 "(segname P1 P2 P3 P5 L45)"

      ## no head or body segmentation for ssu for reference segmentation
      #set sel_ssu_body 
      #set sel_ssu_head 
      set sel_ef2 "(segname P1 P2 P3 P5 L45 SO1 GDP)"
      set sel_ef2G "(segname P2 and resid 219 to 328)"
      set sel_ef2D1 "((segname P1) or (segname P2 and ((resid 67 to 218) or (resid 329 to 345))))"
      set sel_ef2D2 "(segname P2 and (resid 346 to 481))"
      set sel_ef2D3 "(segname P2 and (resid 482 to 558))"
      set sel_ef2D4 "(segname P2 and (resid 559 to 724) or (segname L45 and resid 801 to 842))"
      set sel_ef2D5 "((segname P5) or (segname L45) and (resid 763 to 800))"
      set sel_ef2D1D2 "((segname P1) or (segname P2 and ((resid 67 to 218) or (resid 329 to 345)))) or (segname P2 and (resid 346 to 481))"
      set sel_ef2D3D5 "(segname P2 and (resid 482 to 558)) or ((segname P5) or (segname L45) and (resid 763 to 800))"

  }


  if {$mData0==1} {

set all_domain_sels_name [list "sel_23S" "sel_5S" "sel_16S" "sel_ssu_dom1" "sel_ssu_dom2" "sel_ssu_dom3_major" "sel_ssu_dom3_minor" "sel_eftu_trna" "sel_psite" "sel_esite"]

#set sel_lsu "(chain A to Z) or (chain 0 to 6)"
set sel_23S "chain A"
set sel_5S "chain B"
set sel_16S "chain a"
set sel_ssu_dom1 "(chain a) and (resid 1 to 566)"
set sel_ssu_dom2 "(chain a) and (resid 567 to 912)"
set sel_ssu_dom3_major "(chain a) and (resid 913 to 1396)"
set sel_ssu_dom3_minor "(chain a) and (resid 1397 to 1539)"
set sel_eftu_trna "chain 'z' 'y'"  
set sel_psite "chain v"
      set sel_esite "chain w"
  }
  
  
  
    if {$mData==1} {

set all_domain_sels_name [list "sel_seg1" "sel_seg2"]
set sel_seg1 "chain A to Z 0 to 6"
set sel_seg2 "chain a to v"
  

  }

}

puts "\nDomains are now defined ...\n"

if {$map_segment_only == 0} {  

  if {[info exists PrAxesFile]} {
    puts $PrAxesFile "Initial model (left, State A): $initState_fileno"
  }
  if {[info exists settingsfile]} {
    puts $settingsfile "Initial (Reference) model: $initState_fileno"
  }
}



####### set the color id of the segments and axes

if {$ef2data==1} {
  set colorid_res {23 4 31 1 17 19 20 14 25 17 21 22 28 29} 

  if {$segment_map==1} {
    ##  set colorid_map {23 4 31 1 17}
    ##"colorid_map_reshuffle is set in the script SegmentDensityMap"
    set colorid_res_reshuffle $colorid_map_reshuffle
  }

  #set colorid_res $colorid_res_reshuffle
  
}

if {$mData0==1} {
  #mdata data
  #set colorid_res {23 22 4 18 31 14 1 25} 
  set colorid_res {4 22 4 0 31 14 27 1 12 19} 
}

if {$mData==1} {
  #mdata data
  set colorid_res {22 4} 
}

set axiscolors $colorid_res



#### Select the domains (from the list of all domains defined) for which the inertia tensors are actually
#### going to be computed  

if {$ef2data==1 && $segment_map==1} {
  ### ef2 data
  set which_dom_to_compute {0 1 2 3}
  set compute_allatomonly_dom {3}
}



if {$ef2data==1 && $use_reference_segmentation==1} {

  set which_dom_to_compute {0 1 2 3}
  set compute_allatomonly_dom {3}

  #set which_dom_to_compute {0 1 2 3}
  #set compute_allatomonly_dom {2 3 4 5 6 7 8 9 10}
}

if {$mData0==1} {
  ### mData
  set which_dom_to_compute {0 3 4 5 6 7 8 9}
  set compute_allatomonly_dom {6}

}

if {$mData==1} {
  ### mData
  set which_dom_to_compute {0 1}
  set compute_allatomonly_dom {0 1}

}


set fixDomNum 0
# align LSU to XYZ

set domain_sels_name {}

set num_alldomains [llength $all_domain_sels_name]


puts "\nSaving domains for cartoon representation.\n"
##### save the domain definition for cartoon representation
##### in case the all atom / not all atom definition is used 
for {set dn 0} {$dn < $num_alldomains} {incr dn} {  
  set domain_name [lindex $all_domain_sels_name $dn]
  set str "_cartoon"
  set $domain_name$str [subst $$domain_name]
}


puts "\nSelecting domains for all-atom only inertia tensor calculations, irrespective of other domains.\n"
#### when not using all atom for computing inertia tensors, just use the subset of atoms for 
#### each domain
#### Be careful of the choice, for example, if we choose only P atoms for calculations
#### eEF2 does not contain P atoms

if {$all_atom_tensorcalc==0} {
    
  foreach dnc $which_dom_to_compute {

      ### do this if the domain is not in the list of all atom only list 
      if {[lsearch -exact $compute_allatomonly_dom $dnc]==-1} {
          puts "Domain not in all atom only list:$dnc"
          set domain_n [lindex $all_domain_sels_name $dnc]
          set $domain_name "[subst $$domain_n] and (name $atoms_to_use)"
      }
  }
    
}


puts "\nSelecting domains for inertia tensor calculations.\n"
foreach dn $which_dom_to_compute {
  set domain_name [lindex $all_domain_sels_name $dn]
  lappend domain_sels_name $domain_name
}

set num_domains [llength $domain_sels_name]

set fixDomSel [lindex $all_domain_sels_name $fixDomNum]
# set fixDomSel "sel_23S"

puts "\nSave domains info and settings.\n"

if {$map_segment_only==0} {

  if {[info exists settingsfile]} {
    puts $settingsfile "All domain selection = $all_domain_sels_name"
    puts $settingsfile "Which domains to compute: $which_dom_to_compute"
    puts $settingsfile "Compute all atom only for domain number: $compute_allatomonly_dom"
    puts $settingsfile "Fixed domain number: $fixDomNum"
    if {$segment_map==1} {
      puts $settingsfile "colorid_map = $colorid_map"
    }
    puts $settingsfile "colorid_res = $colorid_res"

  }
}
