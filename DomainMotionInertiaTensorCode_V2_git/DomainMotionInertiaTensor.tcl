##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
##
## This part is the main execution script where all the necessary input models and parameters are set
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: May 2015. Modified:Oct 13,2016; September 08, 2022
##############################################################################################################

########################################################################################### 
# setup the paths to the different scripts and utility programs 
###########################################################################################
set workdir [pwd]

set mcr_root "./"

# add the 'orient' folder to the path 
lappend auto_path Utilities/
lappend auto_path Utilities/la1.0
lappend auto_path Utilities/orient

source Utilities/orient/orient.tcl
package require Orient
namespace import Orient::orient

###########################################################################################  
global Pi eps
set Pi 3.141592653897931
set eps 0.001
####
# contact atomic distance from a point
# set to 6-8 angstrom 
global contactDist
set contactDist 6.0 
###########################################################################################
####### Input dataset
###########################################################################################
#set dataset "splitting-12"
#set dataset "splitting-13"
#set dataset "splitting-23"
#set dataset "splitting-123"
#set dataset "splitting-123_control"
#set dataset "splitting-c1"
#set dataset "splitting-c2"
#set dataset "splitting-c3"

#set dataset "30Sstate12"
#set dataset "asso-3state"

#set dataset "splitting-12"
#set dataset "splitting-13"
#set dataset "70Snohflxcontrol"

#set dataset "new_splitting-c1"
#set dataset "new_splitting-12"
#set dataset "new_splitting-13"

set dataset "new-70S-split-C1-03Oct2023"
#set dataset "new-70S-split-C2-03Oct2023"
#set dataset "new-70S-split-C3-03Oct2023"
#set dataset "new-70S-split-12-03Oct2023"
#set dataset "new-70S-split-13-03Oct2023"

#set PDBdir "../../../DomainMotionProjectData/EF2Data/PDBData/$dataset"
set PDBdir "../pdb_files/$dataset" 
set pdbfile_ext "pdb"

#### for mdff generated structure, we have to use psf+pdb by default
set PSFfile 0

###########################################################################################
# time stamp of run
set timeStamp [clock format [clock seconds] -format -D%Y%m%d-T%H%M%S]

###########################################################################################
####### VMD Rendering
###########################################################################################
set render "snap"
#set render "ray"

if {$render=="snap"} {

    set rendercmd "snapshot"
    set renderfigext "rgb"
    set figext $renderfigext 

} else {

set rendercmd "Tachyon"
set renderfigext "tga"
set figext "bmp" 
}

###########################################################################################
## Clear files and reset variables
source ClearFilesReset.tcl


########################################################################################### 
# source the other required and helper scripts
###########################################################################################

source MapSegment3d.tcl
source ResetDisplayScale.tcl 
source LabelDisplay.tcl

########################################################################################### 
#setup some display settings on VMD
########################################################################################### 
display projection orthographic
display depthcue on
display rendermode GLSL
# for tachyon rendering
display ambientocclusion on

axes location lowerleft

#change the display background to backgroundcol
set backgroundcol "white"
color Display Background $backgroundcol

set vmdaxiscol "black"
color Axes Labels $vmdaxiscol

set vmdwindow 1000
display resize $vmdwindow $vmdwindow
set imscale 0.07

#set molDisplayStyle "NewCartoon"
#set molDisplayStyle "lines"
set molDisplayStyle "Bonds"

set ligDisplayStyle "Bonds"

if {$backgroundcol=="white"} {
    set axislblcol "black"
} else {
    set axislblcol "white"
}
set txtlblcol $axislblcol


set compare_reference 0

if {$compare_reference==1} {

    set repstyle $molDisplayStyle
    set molrepstyle $molDisplayStyle
    set ligrepstyle $ligDisplayStyle

    #set matrl "Glossy"
    set matrl "Transparent"
    #set matrl "Opaque"
    set opqlevel 1
    
} else {
    set repstyle $molDisplayStyle
    set molrepstyle $molDisplayStyle
    set ligrepstyle $ligDisplayStyle
    
    #set matrl "Glossy"
    set matrl "Transparent"
    #set matrl "Opaque"
    set opqlevel 1
}




########################################################################################### 
# Setup some settings for current run
########################################################################################### 

#set segmentationtype "map"
set segmentationtype "pdb"

# settings for volume map segmentation 
set segment_map 0
set map_segment_only 0

set use_reference_segmentation 1

set datasetlabeldisplayoffset -180
set showLabelFig 1


### In case we want to use the previously segmented volumes and 
### just proceed with extraction with new criteria.  
### If we do map_segment_only==0, then the Maps folder is moved inside the resultsdir, 
### so the correct Maps path needs to be used 
### otherwise ExtractWithoutSegmentation = 0 to do the segmentation again 
set ExtractWithoutSegmentation 0

# settings for quaternion calculations
set quaternion_calc 1
# rotation sequence for Euler Angle calculations
# although for now we are using a formula for conversion to Euler angles with the default sequence "ZYX" 
# not using the longer version of the conversion, so the variable rotSequence is not in use for now
set rotSequence "zyx"

# Choose Initial Model State A or Canonical (unrotated) model pdb file number (from sorted list of the files in the folder)
# Arrange the pdb file list such that State A model is the first/last file in the list
#set initState_fileno [expr $total_pdbfiles-1]
set initState_fileno 0

if {$segmentationtype == "map"|| $use_reference_segmentation==1} {

#set outputlog [open $dataset.log w]
cd $PDBdir
set pdb_files [glob "*.$pdbfile_ext"]
set pdb_files [lsort $pdb_files]
    
if {$PSFfile==1} {
    set psf_files [glob *.psf]  
    set psf_files [lsort $psf_files]
}

cd $workdir

set total_pdbfiles [llength $pdb_files]
set nFiles $total_pdbfiles


puts "\nInput Data : $dataset\n"  
    

##### initialize the output files
source InitOutputFiles.tcl



########################################################################################### 
# Processing the pdb structures
########################################################################################### 
for {set nfl 0} {$nfl < $nFiles} {incr nfl} {

    set pdbf [lindex $pdb_files $nfl]	

    if {$PSFfile==1} {
        set psff [lindex $psf_files $nfl]
        puts "\n\n\n\n*** Currently loading $pdbf $psff ***"

    } else {
        puts "\n\n\n\n*** Currently loading $pdbf***"
    }
    mol new $PDBdir/$pdbf  

    if {$PSFfile==1} {
        mol addfile $PDBdir/$psff 
    }

    ########################################################################################### 
    # Perform Segmentation/Generate the domain segments 
    ###########################################################################################
    source GenerateDomainSegments.tcl


    ###########################################################################################
    # Do segmentation and perform Tensor calculations
    ###########################################################################################
    if {$map_segment_only == 0 } {
    source ComputeInertiaTensor.tcl
    }

    
    ###########################################################################################
    # Generate views from different perspective
    ###########################################################################################
    if {$map_segment_only == 0} {  

        source ShowMolRep.tcl
        set show_more_view 1

        source GenerateVMDScreenSnaps.tcl

        puts "\nEnd of tensor analysis for $pdbf structural data ...\n\n\n"
    }
    
    if {$map_segment_only == 0} {
        #### Turn off all the other molecules except the last one in the list 
        if {$nfl != $nFiles-1} {
            mol off top 
        }
    }

    ### to bring the views of loaded molecules to the same scale
    display resetview

    if {$map_segment_only > 0} {
        ### When only segmentation is performed (e.g. initially to check)
        ### stop after the first pdb (used for generating the volume and segmentation)
        break
    }

}


if {$map_segment_only == 0 && $quaternion_calc >0} {
    # do post processing on the saved tensor principal axes for each domain
    source ComputeAbsOrientAllDomain.tcl

}


if {$map_segment_only > 0} {
    puts "\n\nFinished segmentation of map and domain assignment ...\n\n"

    if {$compare_reference==1} {
        # error analysis
        cd ../CompareSegPDBMaps
        source ComparePDBseg.tcl
        cd $workdir
    }

} else {
    puts "\n\nFinished calculations for all pdb structures ...\n\n"
}

}




########################################################################################### 
# Tabulating the rotation angles and rotation axis using the rotation unit quaternion
########################################################################################### 
if {$map_segment_only==0} {
    if {$quaternion_calc >0} {
        source MotionStats.tcl	
    }
}

### assemble the outputs and write to files 
source GenerateOutput.tcl
