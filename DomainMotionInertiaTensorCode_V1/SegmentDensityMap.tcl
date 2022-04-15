  ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part performs segmentation of the Macromolecular structure into the relevant domains
  ## e.g. Ribosome structure segmented into LSU, SSU Head and SSU Body using the volumetric density map 
  ## 
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: May 2015. Modified:Oct 23,2017 
  ##############################################################################################################


  if {$segmentationtype == "map"} {

	puts "\n***Segmentation of the Full PDB structure using the volume map  ...\n"

	set Full_pdbFile $pdbf
	set Full_mapFile "FullMap.mrc"

	### use this representation for showing the segmented domains as extracted from the volume masks later  
        ### no residues are removed/curated for this representation
	set fullsegmolid [mol new $PDBdir/$Full_pdbFile waitfor all]

	if {$PSFfile==1} {
	    mol addfile $PDBdir/$psff 
	}

	########################################################################################### 
        ## 1. Segmentation of the 3D density map generated from the input pdb 
	###########################################################################################
	##### Segmentation cases
	set pixSize_list {1.00 1.25}
	set mapRes_list {1.4 1.8} 
	set Gauss_FilterSize_list {1.2}
	set mapBoxSize_list {380 240}

	set merge_level_list {11 10}
	set merge_size_list {30 50}
      
	
	#### worked: pix 1.00, res 1.4, merge 11, msz 30, gauss filter 1.2,  works fine with 4 segments for 0_ref.pdb 

	set seg_case 0

	# set the parameters for map generation using the executable pdb2mrc
	#set mapBoxSize 240
	#set mapBoxSize 380
	set mapBoxSize [lindex $mapBoxSize_list $seg_case]


	#set pixSize 1.00
	set pixSize [lindex $pixSize_list $seg_case]
	
	
	#set mapRes 1.4
	set mapRes [lindex $mapRes_list $seg_case]

	# set the parameters for the segmentation code 
	# to what depth should the watershed algorithm merge the oversegmented regions in a hierarchical way
	#set merge_level 10
	set merge_level [lindex $merge_level_list $seg_case]
	#set merge_level 11

	# any segment smaller than merge_size will be merged
	#set msz 50
	#set msz 35
	#set msz 30
	set msz [lindex $merge_size_list $seg_case]

	set merge_size [expr $msz*$msz*$msz]


	#walgo = 1 means Watershed algorithm with merging, otherwise basic watershed algorithm  
	set walgo 1

	# extra_merging = 1 means we do an additional step for merging
	# at present this criteria is not very efficient, so better to not use it
	set extra_merging 0

	# smoothing the map using a gaussian filter before running the segmentation 
	#set Gauss_FilterSize 1.2
	set Gauss_FilterSize [lindex $Gauss_FilterSize_list $seg_case]


	set Res_filter [expr $Gauss_FilterSize*$pixSize/(2.0*0.53)]
	set mapRes_effective [expr sqrt($mapRes*$mapRes + $Res_filter*$Res_filter)]

	
	### isovalue used for threholding the map to extract the atoms within that limit
	### higher threshold will remove more residues
	set isoval_thresh 0.5 

	### In case we want to use the previously segmented volumes and 
	### just proceed with extraction with new criteria. 
	### If we do map_segment_only==0, then the Maps folder is moved inside the resultsdir, 
	### so the correct Maps path needs to be used 
	### otherwise ExtractWithoutSegmentation = 0 to do the segmentation again 
	
	set ExtractWithoutSegmentation 0

	if {$ExtractWithoutSegmentation==0} {

	      if {$segment_map == 0} {
		  # choose to use the previously generated segmented map files    
		  puts "\nSegmentation flag is off, looking for already existing density maps from previous 3D volume segmentation ..."

		  set mapSeg_files [glob -nocomplain Maps/SegmentMap_*.mrc]
		  set num_segments [llength $mapSeg_files]
		  # now make the full model at top molecule for analysis
		  mol top $molid_model_full

		  if {$num_segments ==0} {

		      puts "No density files in the map directory, running the segmentation program with segmentation flag on."
		      set segment_map 1

		  } else {
		      # ideally the segmentation program should output the number of the segments,
		      # but since we are using a executable which only provides output through I/O
		      # we calculate the segmented maps through written mrc files
		      puts "The input map was segmented into $num_segments separate volumes ...\n"
		  }
	      }  

	      # the segmented volume slices are viewed with a matlab gui called "Slicer" see Mathworks FEX
	      # wait for several seconds and see if you can move the slider on the gui of the Slicer
	      # then close the Slicer window, 

	      if { $segment_map>0 } {

		  puts "\nNow performing 3D volume segmentation ..."
	  
		  # "MapSegment3d" is a tcl function which calls the script RunMapSegmentation.tcl
		  # RunMapSegmentation.tcl calls a shell script Execute_RunWatershedSegmentMap.sh
		  # Execute_RunWatershedSegmentMap.sh calls teh standalone matlab executable "RunWatershedSegmentMap" for segmentation
		  set v [MapSegment3d $mcr_root $PDBdir/$Full_pdbFile $Full_mapFile $mapBoxSize $mapRes $Gauss_FilterSize $walgo $pixSize $merge_level $merge_size $extra_merging]

		  # check if the map directory 'Maps' exists
		  if {[file isdirectory Maps]==0} {
		      file mkdir Maps
		  }

		  if {[file isdirectory Maps]==1} {	   
		      puts "\nMoving the generated map into the directory 'Maps' ..."
		      # move the segmented maps into the directory called 'Maps'
		      file rename -force {*}[glob *.mrc] Maps 	

		      set mapsegfiles [glob Maps/SegmentMap_*.mrc]
		      set mapSeg_files [lsort $mapsegfiles]
		      set num_segments [llength $mapSeg_files]

		      if {$num_segments ==0} {

		      error "No density files in the map directory. Rerun the program."

		      } else {
			  # ideally the segmentation program should output the number of the segments,
			  # but since we are using a executable which only provides output through I/O
			  # we calculate the segmented maps through written mrc files
			  puts "The input map was segmented into $num_segments separate volumes ...\n"
		      }

		  }

	      }
         }

	########################################################################################### 
        ## 2. Extraction of the residues from within the volume segments for each domain
	###########################################################################################
	# automatically assign correct labels to the segmented domains segments (e.g. based on 
	# the size of the segments?) TO DO 
	# This should only matter if we have more than two segmented domains and we would like to 
	# keep track of the individual domain index so that other properties are assigned to the 
	# correct segments

	set lsusegid  0 
	set bodysegid 3
	set headsegid 1
	set ligsegid  2

	###### set the color id of the segments and axes
	if {$segment_map==1} {
	  ##set colorid_map {23 4 31 18 22} 
	  set colorid_map {23 4 31 1 17}
	}
	
	### for now assign the segments manually
	set mapSeglist_idx [list $lsusegid $bodysegid $headsegid $ligsegid]
	set mapSeglistSort_idx [lsort -real -increasing -indices $mapSeglist_idx]

	set colorid_map_reshuffle {}

	for {set k 0} {$k < $num_segments} {incr k} {
	    lappend colorid_map_reshuffle [lindex $colorid_map [lindex $mapSeglistSort_idx $k]] 
	}
	
	set molid_segments {}
	set residue_seg {}
	set residue_idx {}
	set residue_idx_all {}
	
	###########################################################################################
	# the segmented maps are loaded as individual molecules for visualization
        for {set j 0} {$j < $num_segments} {incr j} {
	    set mapSegfile [lindex $mapSeg_files $j]
	    set molid_seg [mol new $mapSegfile waitfor all]      
	    set repno 0
	    mol modstyle $repno $molid_seg isosurface 0.2 0.0 0.0 0.0 1 1
	    mol modmaterial $repno $molid_seg Transparent
	    mol modcolor $repno $molid_seg colorid [lindex $colorid_map_reshuffle $j]
	    lappend molid_segments $molid_seg
	    puts "Segmented volume map file : [lindex $mapSeg_files $j] loaded as molecule $molid_seg ..." 
	    mol off $molid_seg
	}  
	puts "\n"

	# make the full model as 'top' molecule for analysis
        mol top $molid_model_full

	# load the segmented volumes to the existing pdb molecule molid_model_full, one by one for extracting the pdb coordinates
        puts "\nExtracting the residues from the individual volume segments with a Isoval thresh = $isoval_thresh\n"
	for {set j 0} {$j < $num_segments} {incr j} {
	    puts "Extracting from vol segment: $j"
	    set mapSegfile [lindex $mapSeg_files $j]
	    # above, we set the $molid_model_full as 'top' molecule, so even if we do not specify the molid at the end, 
	    # we still add the volume to the that molecule 
	    set molid_pdbseg [mol addfile $mapSegfile waitfor all $molid_model_full]  
	    
	    # extract the region of the pdb structure within the segmented volume
	    # we can modify the selection criteria for the residues here
	    set selpdbseg_mol_vol [atomselect $molid_model_full "vol$j > $isoval_thresh" frame $f]
	    set res_sel [$selpdbseg_mol_vol get residue] 

	    if { [llength $res_sel]  == 0} {
		puts "\n\nWarning:There were no residues found inside the mask. Possibly the pdb structure and volume map are misaligned or the density map has no model"
		puts "Using the previously derived residue list for the domain, in case it is an issue of the model being in a different coordinate system\n"
	    } else {

		set selpdbseg_mol_res [atomselect $molid_model_full "within 0.01 of residue $res_sel" frame $f]
	      
		lappend residue_seg [$selpdbseg_mol_res get residue] 
		#lappend residue_idx_all [lsort -unique -integer [lindex $residue_seg $j]]

		set residue_seg_j [$selpdbseg_mol_res get residue] 
		set residue_seg_j_sort [lsort -unique -integer $residue_seg_j]
		lappend residue_idx_all $residue_seg_j_sort

	    }

	    #### Use the residue indices from the first input pdb
	    if {$nfl > 0 } {
		set residue_idx_all $residue_idx_all_0 
	    }
	    
	    set residue_idx [lindex $residue_idx_all $j]

	    ###### loading the segment pdbs right after segmentation
	    mol addrep $fullsegmolid
	    set repno [expr $j+1]
	    #set repno [expr $j]
	    mol modselect $repno  $fullsegmolid "residue $residue_idx"  
	    mol modcolor $repno  $fullsegmolid colorid [lindex $colorid_map_reshuffle $j]
	    mol modstyle $repno $fullsegmolid $molrepstyle 0.5  
	    mol modmaterial $repno $fullsegmolid $matrl
	    material change opacity $matrl $opqlevel 
	    mol off $fullsegmolid
    
	    puts "Segmented volume map file : [lindex $mapSeg_files $j] loaded into pdb molecule $molid_model_full ...\n" 
	    if { [llength $res_sel]  != 0} {
		$selpdbseg_mol_vol delete
		$selpdbseg_mol_res delete 
	    }
	}

	#### Save the residue indices for the first input pdb only
	if {$nfl==0} {
	    set residue_idx_all_0 $residue_idx_all 
	    puts "\nSegmented map files saved ...\n"
	}
	
	##### extract the individual volume segments from the segmented full map
        ##### if you have a macromolecular structure with m segments, then do this in a loop over m segments
        ##### The residue indices for each domain is set to residue_idx_vol_j, j = 1 ... m
	for {set j 0} {$j < $num_segments} {incr j} {	
	    # get residues for the segmented domains idx_vol#  
	    set domresvar "residue_idx_vol$j"
	    set $domresvar [lindex $residue_idx_all $j]

        }


	if { $nfl==0 && $map_segment_only==0} {
	    puts $settingsfile "\nSegmentation Parameters::\nMatlabRootDir:$mcr_root\nMap-Box-Size:$mapBoxSize\nPixel-Size:$pixSize\nMap-Resolution:$mapRes\nGauss-Filter-Size:$Gauss_FilterSize\nmapRes-effective:$mapRes_effective\nSegmentation-Algo:$walgo\nMerge-Level:$merge_level\nMerge-Size:$merge_size\nExtraMerging:$extra_merging\nIsoval-thresh=$isoval_thresh\n"
	}

  } else {

	puts "Segmentation using other methods  ... TO DO"
  }