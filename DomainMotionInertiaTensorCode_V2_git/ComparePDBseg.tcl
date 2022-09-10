##############################################################################################################
## Script for Ribosome domain conformation change with EF2/EF-G study using Inertia Tensor Analysis
## This part performs the error calculation of a volume segmented pdb structure compared 
## to a reference segmented pdb. It calculates the precision and recall values for the segmented atoms
## 
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Sep 21, 2016. Modified:Sept 22,2016
##############################################################################################################

########################################################################################### 
### Reference model selections
########################################################################################### 
set sel_lsu_Ref "((segname 25S 5S 5D8S \"L.*\") and (not segname P1 P2 P5 L45 SO1 GDP))" 
set sel_ssu_Ref "((segname 18S \"S.*\" RACK) and (not segname P1 P2 P5 L45 SO1 GDP))"
set sel_ef2_Ref "(segname P1 P2 P5 L45 SO1 GDP)"
set sel_lsu_ssu_Ref "((segname 25S 5S 5D8S \"L.*\" 18S \"S.*\" RACK) and (not segname GDP SO1 P1 P2 P5 L45))"

########################################################################################### 
### Classified/Segmented model selections
########################################################################################### 
set sel_lsu "((residue $residue_idx_lsu) and (not segname P1 P2 P5 L45 SO1 GDP))"
set sel_ssu_body "((residue $residue_idx_body) and (not segname P1 P2 P5 L45 SO1 GDP))"
set sel_ssu_head "((residue $residue_idx_head) and (not segname P1 P2 P5 L45 SO1 GDP))"
set sel_ssu "((residue [concat $residue_idx_body $residue_idx_head]) and (not segname P1 P2 P5 L45 SO1 GDP))"



###########################################################################################
if {$compare_reference==1} {
	    set Full_pdbFile $pdbf
    
	    set fullcorrectedmolid [mol new $PDBdir/$Full_pdbFile waitfor all]
	    mol addfile $PDBdir/$psff 
	    # set the representation type 
	    mol modcolor 0 $fullcorrectedmolid colorid 10 
	    mol modstyle 0 $fullcorrectedmolid Bonds 0.5  
	    mol modmaterial 0 $fullcorrectedmolid Opaque
	    mol off $fullcorrectedmolid
 }

#mol delrep 0 $molid_model_full
mol top $molid_model_full

if {$compare_reference==1} {
  package require clonerep
  ::CloneRep::clone_reps $molid_model_full $fullcorrectedmolid
}

	
 if {$map_segment_only == 1 && $compare_reference==1} {
	### compare with the reference class2-40S,60S,EF2 
	cd ../Class2-Reference
	set reflsumolid [mol new class2-60S.pdb waitfor all]
	set refssumolid [mol new class2-40S.pdb waitfor all]
	set refef2molid [mol new class2-EF2.pdb waitfor all]
	

	set sel_reflsu "all"
	set sel_refssu "all"
	set sel_refef2 "all"

	set repno_reflsu 1
	set repno_refssu 1
	set repno_refef2 1 

	set reflsucolid 20
	set refssucolid 20
	set refef2colid 20

	mol addrep $reflsumolid
	mol addrep $refssumolid
	mol addrep $refef2molid

	mol modselect $repno_reflsu $reflsumolid $sel_reflsu  
	mol modcolor $repno_reflsu $reflsumolid colorid $reflsucolid
	mol modstyle $repno_reflsu $reflsumolid $molrepstyle 0.5    
	mol modmaterial $repno_reflsu $reflsumolid $matrl
	material change opacity $matrl $opqlevel 

	mol modselect $repno_refssu $refssumolid $sel_refssu  
	mol modcolor $repno_refssu $refssumolid colorid $refssucolid
	mol modstyle $repno_refssu $refssumolid $molrepstyle 0.5    
	mol modmaterial $repno_refssu $refssumolid $matrl
	material change opacity $matrl $opqlevel 

	mol modselect $repno_refef2 $refef2molid $sel_refef2  
	mol modcolor $repno_refef2 $refef2molid colorid $refef2colid
	mol modstyle $repno_refef2 $refef2molid $molrepstyle 0.5    
	mol modmaterial $repno_refef2 $refef2molid $matrl
	material change opacity $matrl $opqlevel 


	mol delrep 0 $reflsumolid
	mol delrep 0 $refssumolid	
        mol delrep 0 $refef2molid

	cd ../CompareSegPDBMaps
 }
  


########################################################################################### 
### We can only compare the atoms for lsu and ssu 
### get the atom indices from the atom selection for the lsu and ssu 
###
########################################################################################### 

set atomIndexList_lsu_Ref_mol [atomselect $molid_model_full "$sel_lsu_Ref"]
set atomIndexList_ssu_Ref_mol [atomselect $molid_model_full "$sel_ssu_Ref"] 

set atomIndexList_lsu_mol [atomselect $molid_model_full "$sel_lsu"]
set atomIndexList_ssu_mol [atomselect $molid_model_full "$sel_ssu"]


set atomIndexList_lsu_Ref [$atomIndexList_lsu_Ref_mol list] 
set atomIndexList_ssu_Ref [$atomIndexList_ssu_Ref_mol list] 

set atomIndexList_lsu [$atomIndexList_lsu_mol list] 
set atomIndexList_ssu [$atomIndexList_ssu_mol list] 


set num_atom_lsu_Ref [$atomIndexList_lsu_Ref_mol num]
set num_atom_ssu_Ref [$atomIndexList_ssu_Ref_mol num]

set num_atom_lsu [$atomIndexList_lsu_mol num]
set num_atom_ssu [$atomIndexList_ssu_mol num]

puts "\nnum_atom_lsu_Ref = $num_atom_lsu_Ref"
puts "\nnum_atom_ssu_Ref = $num_atom_ssu_Ref"
set num_atom_all_Ref [expr ($num_atom_lsu_Ref+$num_atom_ssu_Ref)]
puts "\nnum_atom_all_Ref = $num_atom_all_Ref"

puts "\n\nnum_atom_lsu = $num_atom_lsu"
puts "\nnum_atom_ssu = $num_atom_ssu"
set num_atom_all [expr ($num_atom_lsu+$num_atom_ssu)]
puts "\nnum_atom_all = $num_atom_all"



### Difference in reference set and segmented set
set atomIndexList_lsuRef_not_lsu_mol [atomselect $molid_model_full "(($sel_lsu_Ref) and (not $sel_lsu))"] 
set atomIndexList_ssuRef_not_ssu_mol [atomselect $molid_model_full "(($sel_ssu_Ref) and (not $sel_ssu))"]
 
set atomIndexList_not_lsuRef_lsu_mol [atomselect $molid_model_full "(not ($sel_lsu_Ref) and ($sel_lsu))"] 
set atomIndexList_not_ssuRef_ssu_mol [atomselect $molid_model_full "(not ($sel_ssu_Ref) and ($sel_ssu))"]

set FN_lsu_atomIdx [$atomIndexList_lsuRef_not_lsu_mol list]
set FN_ssu_atomIdx [$atomIndexList_ssuRef_not_ssu_mol list]

set FP_lsu_atomIdx [$atomIndexList_not_lsuRef_lsu_mol list]
set FP_ssu_atomIdx [$atomIndexList_not_ssuRef_ssu_mol list]

set sel_FN_lsu "index $FN_lsu_atomIdx" 
set sel_FP_lsu "index $FP_lsu_atomIdx"


set sel_FN_ssu "index $FN_ssu_atomIdx" 
set sel_FP_ssu "index $FP_ssu_atomIdx"


mol addrep $molid_model_full
mol addrep $molid_model_full
mol addrep $molid_model_full
mol addrep $molid_model_full

set repnoFN_lsu 3
set repnoFP_lsu 4

set FN_lsucolid 20
set FP_lsucolid 13



set repnoFN_ssu 5
set repnoFP_ssu 6

set FN_ssucolid 20
set FP_ssucolid 13


mol modselect $repnoFN_lsu $molid_model_full $sel_FN_lsu  
mol modcolor $repnoFN_lsu $molid_model_full colorid $FN_lsucolid
mol modstyle $repnoFN_lsu $molid_model_full $molrepstyle 0.5    
mol modmaterial $repnoFN_lsu $molid_model_full $matrl
material change opacity $matrl $opqlevel 


mol modselect $repnoFP_lsu $molid_model_full $sel_FP_lsu  
mol modcolor $repnoFP_lsu $molid_model_full colorid $FP_lsucolid
mol modstyle $repnoFP_lsu $molid_model_full $molrepstyle 0.5    
mol modmaterial $repnoFP_lsu $molid_model_full $matrl
material change opacity $matrl $opqlevel 




mol modselect $repnoFN_ssu $molid_model_full $sel_FN_ssu  
mol modcolor $repnoFN_ssu $molid_model_full colorid $FN_ssucolid
mol modstyle $repnoFN_ssu $molid_model_full $molrepstyle 0.5    
mol modmaterial $repnoFN_ssu $molid_model_full $matrl
material change opacity $matrl $opqlevel 


mol modselect $repnoFP_ssu $molid_model_full $sel_FP_ssu  
mol modcolor $repnoFP_ssu $molid_model_full colorid $FP_ssucolid
mol modstyle $repnoFP_ssu $molid_model_full $molrepstyle 0.5    
mol modmaterial $repnoFP_ssu $molid_model_full $matrl
material change opacity $matrl $opqlevel 




### intersection of ssu and ssu_ref

set atomIndexListIntersect_lsuRef_lsu_mol [atomselect $molid_model_full "(($sel_lsu_Ref) and ($sel_lsu))"] 
set atomIndexListIntersect_ssuRef_ssu_mol [atomselect $molid_model_full "(($sel_ssu_Ref) and ($sel_ssu))"] 

set atomIndexListIntersect_lsuRef_lsu [$atomIndexListIntersect_lsuRef_lsu_mol list]
set atomIndexListIntersect_ssuRef_ssu [$atomIndexListIntersect_lsuRef_lsu_mol list]

set num_atomIntersect_lsu [$atomIndexListIntersect_lsuRef_lsu_mol num]
set num_atomIntersect_ssu [$atomIndexListIntersect_ssuRef_ssu_mol num]
set num_atom_all_intersect [expr ($num_atomIntersect_lsu+$num_atomIntersect_ssu)]

puts "\n\nnumIntersect_lsu = $num_atomIntersect_lsu"
puts "\nnumIntersect_ssu = $num_atomIntersect_ssu"
puts "\n\nnum_atom_all_intersect = $num_atom_all_intersect"

### write a function for precision/recall/F-measure .. TO DO 

set alpha 0.5
set precision_lsu [expr $num_atomIntersect_lsu/($num_atom_lsu + 0.0)]
set recall_lsu [expr $num_atomIntersect_lsu/($num_atom_lsu_Ref + 0.0)]
set Fmeasure_lsu [expr 1.0/($alpha*(1.0/$precision_lsu)+(1-$alpha)*(1.0/$recall_lsu))]


set precision_ssu [expr $num_atomIntersect_ssu/($num_atom_ssu + 0.0)]
set recall_ssu [expr $num_atomIntersect_ssu/($num_atom_ssu_Ref + 0.0)]
set Fmeasure_ssu [expr 1.0/($alpha*(1.0/$precision_ssu)+(1-$alpha)*(1.0/$recall_ssu))]

set precision_allSeg [expr ($num_atomIntersect_lsu + $num_atomIntersect_ssu)/($num_atom_lsu + $num_atom_ssu + 0.0)]
set recall_allSeg [expr ($num_atomIntersect_lsu + $num_atomIntersect_ssu)/($num_atom_lsu_Ref + $num_atom_ssu_Ref + 0.0)]
set Fmeasure_allSeg [expr 1.0/($alpha*(1.0/$precision_allSeg)+(1-$alpha)*(1.0/$recall_allSeg))]

puts "\n\nprecision_lsu = $precision_lsu"
puts "\nrecall_lsu = $recall_lsu"
puts "\nFmeasure_lsu = $Fmeasure_lsu"

puts "\n\nprecision_ssu = $precision_ssu"
puts "\nrecall_ssu = $recall_ssu"
puts "\nFmeasure_ssu = $Fmeasure_ssu"

puts "\n\nprecision_allSeg = $precision_allSeg"
puts "\nrecall_allSeg = $recall_allSeg"
puts "\nFmeasure_allSeg = $Fmeasure_allSeg"

set segaccuracyfilename "SegAccuracy-$dataset$timeStamp.dat"
set segaccuracyfile [open $segaccuracyfilename {RDWR CREAT}]
puts $segaccuracyfile "Segmentation Parameters:\nmapBoxSize=$mapBoxSize\nmapRes=$mapRes\npixSize=$pixSize\nmerge_level=$merge_level\nmerge_size=$merge_size\nGauss_FilterSize=$Gauss_FilterSize\nisoval_thresh=$isoval_thresh\n"
puts $segaccuracyfile "\nSegmentation Statistics:\nnum_atom_lsu_Ref = $num_atom_lsu_Ref\nnum_atom_ssu_Ref = $num_atom_ssu_Ref\nnum_atom_lsu = $num_atom_lsu\nnum_atom_ssu = $num_atom_ssu"
puts $segaccuracyfile "\nnumIntersect_lsu = $num_atomIntersect_lsu\nnumIntersect_ssu = $num_atomIntersect_ssu"
puts $segaccuracyfile "\n\nprecision_lsu = $precision_lsu\nrecall_lsu = $recall_lsu\nFmeasure_lsu = $Fmeasure_lsu\n\nprecision_ssu = $precision_ssu\nrecall_ssu = $recall_ssu\nFmeasure_ssu = $Fmeasure_ssu\n\nprecision_allSeg = $precision_allSeg\nrecall_allSeg = $recall_allSeg\nFmeasure_allSeg = $Fmeasure_allSeg"

close $segaccuracyfile

########################################################################################### 
### Write the selections into pdb files
###########################################################################################
[atomselect $molid_model_full $sel_lsu_Ref] writepdb $pdbf-lsu-Ref.pdb
[atomselect $molid_model_full $sel_ssu_Ref] writepdb $pdbf-ssu-Ref.pdb

[atomselect $molid_model_full $sel_lsu] writepdb $pdbf-lsu.pdb
[atomselect $molid_model_full $sel_ssu] writepdb $pdbf-ssu.pdb

set pdbFile_lsu $pdbf-lsu.pdb
set pdbFile_ssu $pdbf-ssu.pdb

#set pdbFile_lsu class2-60S.pdb
#set pdbFile_ssu class2-40S.pdb

#set pdbFile_lsu class2.pdb-lsu.pdb
#set pdbFile_ssu class2.pdb-ssu.pdb

#mol new $pdbf-lsu.pdb
#mol new $pdbf-ssu.pdb


###########################################################################################
### Compare with maps 
###########################################################################################
set compare_maps 0

if {$compare_maps==1} {
    set mapFile_lsu $pdbFile_lsu.mrc
    set mapFile_ssu $pdbFile_ssu.mrc


    set mapBoxSize_m 240
    set pixSize_m 1.25
    set mapRes_m 2.0

    #exec pdb2mrc $pdbFile_lsu $mapFile_lsu box=$mapBoxSize_m apix=$pixSize_m res=$mapRes_m center het 
    #exec pdb2mrc $pdbFile_ssu $mapFile_ssu box=$mapBoxSize_m apix=$pixSize_m res=$mapRes_m center het 


    exec e2pdb2mrc.py $pdbFile_lsu $mapFile_lsu --box $mapBoxSize_m --apix $pixSize_m --res $mapRes_m --het --verbose 2
    exec e2pdb2mrc.py $pdbFile_ssu $mapFile_ssu --box $mapBoxSize_m --apix $pixSize_m --res $mapRes_m --het --verbose 2

    ###########################################################################################
    mol new $pdbFile_lsu.mrc
    mol new $pdbFile_ssu.mrc
}

set comparedir "class2-RefVs-Segmented-b$mapBoxSize-p$pixSize-r$mapRes-f$Gauss_FilterSize-ml$merge_level-msz$msz"
file mkdir $comparedir
file rename -force {*}[glob *.dat] $comparedir
file rename -force {*}[glob *.pdb] $comparedir
