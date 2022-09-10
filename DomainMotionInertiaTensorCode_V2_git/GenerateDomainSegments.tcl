##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
## 
## This script generates the domain segmentation  
## 
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: May 2015. Modified:April 20,2017 
##############################################################################################################


set molid_top [molinfo top] 
set molid_model_full [expr $molid_top]
puts "\nMolid of the full molecular machine pdb model: $molid_model_full" 

###########################################################################################
### First step is to get the pdb structure segmented out, using the volume 
### segmentation or otherwise
###########################################################################################  
set f $nfl

if {$segmentationtype=="map" && $use_reference_segmentation==0} {

  # delete the previously loaded segmented molecules
  if {$f>0} {

      mol delete $fullsegmolid
      for {set k 0} {$k < $num_segments } {incr k} {
          mol delete [lindex $molid_segments $k]
      }
    }
}
  

display update 
    
#go to the current frame
animate goto $f

#frame #
puts "\n\n***Analysing molecule in frame : $f"
# convert frame number to time (ns)
set time [expr $f*(.001)]
    

### First segment out the head and body of the small subunit using automatic map segmentation or 
### (check the segmented volumes using segementation only mode with map_segment_only = 1 
### and if they are fine then rerun again with map_segment_only = 0)
if {$segmentationtype=="map"} {

  #### do this only for first input model, for rest of the model, the residues from first model 
  #### are used for the corresponding domain 
  if {$nfl==0} {
      source SegmentDensityMap.tcl
  }
}
# set the full model molecule to "top" for further analysis
mol top $molid_model_full

###########################################################################################  
    
# use the segmented residue ids from first segmentation and use it for rest of the subsequent pdbs from the dataset
# since we expect the pdbs in the same dataset to have the same order. This would ensure we get the principal axes of the 
# segments to be consistent with the residues constituting the individual segments
# always verify the segment index visually

# As such the new residue list are created every time the file SegmentRibosomeCryoEMmap.tcl is called
# In general the updated list should be used for problems where new domain residue list is required
# For our case we want the residue list to be same as the first model we analyze. 
# So we ignore the updated list of residues for each domain even if they are calculated. 
# Just customize the script for your own purpose

# Also if for some reason, the segmentation doesn't work properly and the domain residues cannot be obtained automatically
# we can set the residues for each domain manually in this section, if they are available from other sources 
# We were not able to obtain the segmentation for eEF2, and used the know residues 
if {$nfl==0} {
    source SetDomainSegments.tcl 
}

###
### domain_sels_name is set in the script setDomainSegments.tcl



# color and draw the molecule segments appropriately with additional representation
# manually turn off/on the full representation
source MolRepColor.tcl
display update on


if {$map_segment_only>0} {
    return
}


