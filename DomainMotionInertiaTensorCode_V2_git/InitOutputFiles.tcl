
##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
## 
## This part sets up the files to write the output files
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Nov 2015. Modified:April 20,2017 
## 
##
##############################################################################################################


# Setup the results directory
set delrawimgfile 0


set resultsdir "$resultsrootdir/$dataset$timeStamp"

#set resfileaffix [file rootname [file tail $resultsdir]]
set resfileaffix $dataset$timeStamp
	
if {$map_segment_only==0} {
	
	file mkdir $resultsrootdir

	file mkdir $resultsrootdir/tmp
	
	file mkdir $resultsdir
	
	if {$delrawimgfile==0} {
		file mkdir $resultsdir/images
		file mkdir $resultsdir/rawimages
	}

	file mkdir $resultsdir/gifmovies
	file mkdir $resultsdir/datfiles

	########################################################################################### 
	# Setup the file names for the results to be written
	########################################################################################### 
	set settingsfilename "$resultsrootdir/tmp/Settings_$resfileaffix.dat"
	set settingsfile [open $settingsfilename {RDWR CREAT}]
	#puts $settingsfile "\nDataset : $dataset\nNumber of Files: $nFiles\nFileNames : $pdb_files\nFixed Molecule : $fixmol\nLigand present: $ligand_exist\nFull Ligand:$full_ligand"
	puts $settingsfile "\nDataset : $dataset\nNumber of Files: $nFiles\nFileNames : $pdb_files\n"
	#puts $settingsfile "Full SSU (0=Segmented):$full_SSU\nAll atom tensor: $all_atom_tensorcalc\nQuaternion calculations done : $quaternion_calc\nRotation Sequence for Quaternion to Euler: $rotSequence\n"
	


	set TensorAnglesFilename "$resultsrootdir/tmp/ThetaPsi_$resfileaffix.dat"
	set TensorAnglesFile [open $TensorAnglesFilename {RDWR CREAT APPEND}]
	puts $TensorAnglesFile "\t\ttheta_axis3 psi_axis3\ttheta_axis2 psi_axis2\ttheta_axis1 psi_axis1\t\n"

	set PrAxesFilename "$resultsrootdir/tmp/PrAxesinfo_$resfileaffix.dat"
	set PrAxesFile [open $PrAxesFilename {RDWR CREAT APPEND}]


	set QuatEulerAnglesFilename "$resultsrootdir/tmp/QPsiThetaPhi_$resfileaffix.dat"
	set QuatEulerAnglesFile [open $QuatEulerAnglesFilename {RDWR CREAT APPEND}]

	set QuatRotAxisFilename "$resultsrootdir/tmp/QRotAxisAngle_$resfileaffix.dat"
	set QuatRotAxisFile [open $QuatRotAxisFilename {RDWR CREAT APPEND}]

}
