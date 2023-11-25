##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
## 
## This part first runs the segmentation script and then computes the Inertia Tensor 
## and Principal Axes for each domain defined by segmentation or otherwise
## 
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: May 2015. Modified:October 23,2023 
##############################################################################################################

## Note:
## Customize this script to compute the inertia tensors and principal axes for any particular 
## domain of a structure. 
## Note: This section can most likely be made lot simpler, if simple criteria for adjusting the direction of 
## of the computed principal axes can be determined. (To be done later )  
## Also some parts of the code here can be further modularized for clarity
##############################################################################################################


###########################################################################################  
## Inertia Tensor and Principal Axes Computation
########################################################################################### 

### do this for all selected domains
### reference and all target domains
set domain_sels_mol {}

puts "frame: $f"
for {set dn 0} {$dn < $num_domains} {incr dn} {
    set domain_name [lindex $domain_sels_name $dn]
    set domain_sel [subst $$domain_name]
    puts "Domain selection name:$domain_name"
    set domain_sel_mol [atomselect $molid_model_full $domain_sel frame $f]
    lappend domain_sels_mol $domain_sel_mol
    unset domain_sel_mol

}
set all_mol {}
set all_mol [atomselect $molid_model_full "all" frame $f]
### there is another variable 'fixDomNum',which is for the entire list
### fixedDomain_ind would be for the list of domains to be analyzed only
set fixedDomain_ind [lsearch $domain_sels_name $fixDomSel]
set fixedDomain_name [lindex $domain_sels_name $fixedDomain_ind]

set fixedMol {}
set fixedMol [lindex $domain_sels_mol $fixedDomain_ind]

        
###########################################################################################
    
########################################################################################### 
##
## Align the principal axes of reference molecule to the XYZ-axis 
##  
########################################################################################### 

set calc_tensor 1

# October 23, 2023
# if align_force_tarmol_to_refmoldom=1, 
# we are forcing to alignment of the domain/molecule using [measure fit]
# otherwise the reference domains are aligned with just the Principal axis

# use rmsd fit for reference domain alignment
#set align_force_tarmol_to_refmoldom 1

# use the principal axes to align instead 
set align_force_tarmol_to_refmoldom 0 
    


if {$calc_tensor > 0} {
    # remove all the previous graphics objects
    graphics top delete all  

    if {$all_atom_tensorcalc == 0} {

        puts "\nConsidering only the $atoms_to_use atoms for fast computing the Inertia Tensors ...\n"

    } else {

        puts "\nConsidering all atoms for computing the Inertia Tensors (so a bit slower) ...\n"

    }
    
    ## Step 1. Calculate the initial Principal Moment of Inertia for the reference subunit unit.
    ## Align the reference principal axes to the XYZ cartesian coordinate system. This would help us to characterize the motion of the 
    ## othe subunit relative to a fixed axis system

    # translate the full molecule such that the centroid of reference domain is at origin (0,0,0)

    puts  "Measuring center of the fixed molecule"
    set fixedmol_c [measure center $fixedMol weight mass] 
    puts "\n\nIntial fixed molecule centroid : $fixedmol_c"

    set all_c [measure center $all_mol weight mass] 
    puts "Initial full mol centroid : $all_c"
    
    set cdist [vecscale -1.0 $fixedmol_c]
    puts "moveby coordinates : $cdist"
    $all_mol moveby $cdist

    # align few more times, for better accuracy of principal axes alignment to XYZ, to origin at (0,0,0)
    # idealy only one "moveby" displacement fitting should work, but not sure why it takes few round to achieve 
    # the alignment here, alternatively use a tolarence level to stop the iteration
    for {set i 0} {$i<2} {incr i} {
        set fixedmol_c [measure center $fixedMol weight mass] 	  
        set cdist [vecscale -1.0 $fixedmol_c] 
        puts "moveby coordinates : $cdist"    
        $all_mol moveby $cdist
    }


    unset all_mol
    set all_mol [atomselect $molid_model_full "all" frame $f]
    
    puts "\n***Computing the principal axes for the reference molecule $fixedDomain_name ...\n"
    #(a) First Align the molecule-to-be-fixed, e.g Principal Axis 3 to the Z-axis {0 0 1} 
    set adjust_direction 0

    # December 28, 2022. Manually flip axis/axes at this step so that output visualization of axes is consistent with other pdb models in case the axes are 
    # randomly flipped due to small changes 
    #set man_flip1_a1a2a3 {-1 -1 1}
    #set man_flip2_a1a2a3 {-1 -1 1}
    set man_flip1_a1a2a3 {1 1 1}
    set man_flip2_a1a2a3 {1 1 1}
    # no flip if 1 1 1
    set Iref {}
    set Jref {}

    set Iref [Orient::calc_principalaxes $fixedMol $adjust_direction $man_flip1_a1a2a3]
    puts "\nIref=$Iref"
    
    set A [orient $fixedMol [lindex $Iref 2] {0 0 1}]


    #### align entire molecule according to the alignment of the fixedMol 
    set all_mol [atomselect $molid_model_full "all" frame last]
    $all_mol move $A
        
    #(b) Next Align the molecule-to-be-fixed, Principal Axis 2 to the Y-axis {0 1 0} 

    set Jref [Orient::calc_principalaxes $fixedMol $adjust_direction $man_flip2_a1a2a3] 
    puts "\nJref=$Jref"



    set B [orient $fixedMol [lindex $Jref 1] {0 1 0}]

    #### align entire molecule accroding to the alignment of the fixedMol 
    $all_mol move $B
        
    # October 23, 2023
    # if align_force_tarmol_to_refmoldom=1, we are forcing to align the fixed domain of second pdb to the fixed domain of 
    # the first pdb , even when the principal axes are very slightly mis-aligned
    # Note that approximating the entire big domain by three principle axes may not always represent the 
    # domain atom distribution the way we expect.
    puts "fixed domain selection: [subst $$fixDomSel]"
    if {$align_force_tarmol_to_refmoldom > 0} {
        set fixedDomSel [subst $$fixDomSel]
        set fixedDomMol [atomselect $molid_model_full $fixedDomSel frame $f]
        if {$nfl==0} {
            set fixedDomMolref $fixedDomMol
        } else {
            set fixedDomMoltar $fixedDomMol
            set transMat [measure fit $fixedDomMoltar $fixedDomMolref]
            $all_mol move $transMat
        }
        
    }

    # do we need to re select the atoms afer the alignment?
    set fixdomain_cen [measure center $fixedMol weight mass]  
    puts "Fixed Domain centroid : $fixdomain_cen\n"



    ######################################################
    ### Alignment done and then we draw the principal axes for the different selected domains
    puts "\nReference segment $fixedDomain_name, Principal Axes are now aligned with X Y Z.\n"


    #########################################################
    ## Now get the atom selection for each domain of interest
    #########################################################
    # modified Sep 2022
    set All_mol [atomselect $molid_model_full "all" frame $f]

    set domain_sels_align_mol {}

    if {$nfl==0} {
        set allmodel_domains_cen {}
    }

    set domains_cen {}
    for {set dn 0} {$dn < $num_domains} {incr dn} {

        set domain_name [lindex $domain_sels_name $dn]
        set domain_sel [subst $$domain_name]

        puts "Domain selection name:$domain_name"
        set domain_sel_align_mol [atomselect $molid_model_full $domain_sel frame $f]
        lappend domain_sels_align_mol $domain_sel_align_mol

        set domain_cen [measure center $domain_sel_align_mol weight mass]  
        puts "New Domain centroid : $domain_cen\n"
        lappend domains_cen $domain_cen
        
    } 
    lappend allmodel_domains_cen $domains_cen

    # added Sep 2022
    if {$nfl==$initState_fileno} {
        set initmodel_domains_cen [lindex $allmodel_domains_cen $nfl]
    }	


    # To save time, we can skip the computation of the pricipal axes for reference domain

    ## Step 2. Re-calculate the principal axes of the aligned domains.

    ### Careful with the direction of principal axes, some of them flip for a different model even with a slight change in the mass distribution
    ### I have tried to enforce a simple rule for keeping the direction of the principal axes consistent for a particular domain in two different
    ### model, but it seems like it doesn't work all the time and it seems to be difficult to find a general criteria 


    puts "\n\n***Computing and drawing principal axes for XYZ aligned Referenece domain $fixedDomain_name ...\n"

    set Principal_axes_domain_names {}

    for {set dn 0} {$dn < $num_domains} {incr dn} {
        if {$dn==$fixedDomain_ind} {
            set adjust_direction 0

        } else {
            set adjust_direction 1
        }

        
        set domain_sel_align_mol [lindex $domain_sels_align_mol $dn]
        set domain_name [lindex $domain_sels_name $dn]
        puts "\n\n\nDomain selection name:$domain_name"

        puts "***Computing and drawing principal axes for $domain_name ...\n"
        set adjust_scale 0

        set dom_comp_num [lindex $which_dom_to_compute $dn]
        set axiscol [lindex $axiscolors $dom_comp_num]


        set ts "_Paxes"
        set dom_tensor $domain_name$ts
        ### note that $dom_tensor (instead of just dom_tensor) creates the variable for storing the principal axes for the domain
        ### so the syntax "set $dom_tensor [ ] ... " is correct
        
 
        set manual_flip_tensor 1
        set adjust_direction 0
        ### specify the models $nfl for which the axes needs to be flipped: $nfl==0, $nfl==1 etc
        if {$nfl==1} {
           set manual_flip_tensor 1
        } 

        if {$manual_flip_tensor > 0} { 

               puts "do manual flip of $dom_tensor"
                   # December 28, 2022. Manually flip axis/axes at this step so that output visualization of axes is consistent with other pdb models in case the axes are randomly flipped due to small changes 
               if {$dn==$fixedDomain_ind} {
                  set man_flip_a1a2a3 {-1 1 -1}
               } else {
                  puts "pdbMol: $nfl, domain: $dn"
                  set man_flip_a1a2a3 {-1 1 1}
               }	  
               
               
         } else {
               set man_flip_a1a2a3 {1 1 1}
         }

        
        set $dom_tensor [draw principalaxes $domain_sel_align_mol $axiscol $axislblcol $adjust_direction $man_flip_a1a2a3 $adjust_scale]
        
	
        lappend Principal_axes_domain_names $dom_tensor
    }
}


puts "\n\n************************************************************"
puts "*** Domain motion results ***"
puts "************************************************************\n\n"

###########################################################################################
# Polar angles Theta and Psi directly from the principal axes coordinates
###########################################################################################


### Check the theta and phi angle calculation formula from the axes, seems to be incorrect 
puts "*** Polar angles from Tensor Principal Axes ***\n\n"
source TensorPolarAngles.tcl

puts $TensorAnglesFile "$pdbf"	

for {set dn 0} {$dn < $num_domains} {incr dn} {

    set domain_name [lindex $domain_sels_name $dn]
    set dom_t [lindex $Principal_axes_domain_names $dn]
    set domT [subst $$dom_t]

    puts $TensorAnglesFile "\n$domain_name"	
    TensorPolarAngles $domT $domain_name $TensorAnglesFile
}

puts $TensorAnglesFile "\n\n"
puts $TensorAnglesFile "$nfl $pdbf"

###########################################################################################
# Euler Angle Calculations for Local coordinate system (LCS) expressed in terms of the 
# World Coordinate System(WCS) XYZ 
# When $L is the fixed reference domain it is the WCS or XYZ
# Just shown for a selected domain (LCS) rotation relative to a reference domain (WCS)
# Others can be calculated similarly
###########################################################################################
#### domT contains the principal axes for the listed domain in "$domain_sels_name"
source ComputeEuler.tcl
set EulerDomainT [ComputeEuler $domT]

###########################################################################################

###########################################################################################
# Calculating the selected domain motions
########################################################################################### 

if {$quaternion_calc > 0} {

    # just use the default ordering of the principal axes for each domain to compute the transformation
    # since we are assuming that the overall directions will not change much from one state to other and even if it does 
    # we can argue about the large rotations because of the change in mass ? but overall the principal 
    # axes order should remain same. Manually check the ordering, in case of inconsistencies  > TO DO 

    puts "*** Principal axes storage for Model $nfl ***"
    set axisorder "default"
    puts "\n\nPrincipal Axes Ordering for Computing the Rotational angles: $axisorder\n"

    puts $PrAxesFile "\n$pdbf"

    ### prepare the variables for storing the axes for models (state a-->b) for each domain
    if {$nfl==0} {

        for {set dn 0} {$dn < $num_domains} {incr dn} {

            set domain_name [lindex $domain_sels_name $dn]
            set rvar Ref_axes_nvh_all_model_$domain_name
            set tvar Tar_axes_nvh_all_model_$domain_name
                puts "rvar :$rvar"
            puts "tvar :$tvar"
            set $rvar {} 
            set $tvar {}
        }

    } else {

        puts "Storage variables already initialized for earlier model 0, this is pdb model $nfl."
    }



    set RefT [subst $[lindex $Principal_axes_domain_names $fixedDomain_ind]]

    for {set dn 0} {$dn < $num_domains} {incr dn} {


        set domain_name [lindex $domain_sels_name $dn]

        set TarT [subst $[lindex $Principal_axes_domain_names $dn]]

        ### three principal axes in n, v, h directions, relative to the reference domain
        if {$axisorder=="default"} {
            set RefT_n 0  
            set RefT_v 1
            set RefT_h 2

            set TarT_n 0
            set TarT_v 1
            set TarT_h 2 

        }

        set RefT_axis_n [lindex $RefT $RefT_n]
        set RefT_axis_v [lindex $RefT $RefT_v]
        set RefT_axis_h [lindex $RefT $RefT_h]

        set Ref "[concat $RefT_axis_n $RefT_axis_v $RefT_axis_h]"
        puts "\nRefT_n = $RefT_n, RefT_v = $RefT_v, RefT_h = $RefT_h"
        puts $PrAxesFile "\nRefT_n = $RefT_n, RefT_v = $RefT_v, RefT_h = $RefT_h\n"
        
        set Ref_axis_n $RefT_axis_n
        set Ref_axis_v $RefT_axis_v
        set Ref_axis_h $RefT_axis_h

        set TarT_axis_n [lindex $TarT $TarT_n]
        set TarT_axis_v [lindex $TarT $TarT_v]
        set TarT_axis_h [lindex $TarT $TarT_h]

        set Tar "[concat $TarT_axis_n $TarT_axis_v $TarT_axis_h]"
        puts "TarT_n = $TarT_n, TarT_v = $TarT_v, TarT_h = $TarT_h"
        puts $PrAxesFile "TarT_n = $TarT_n, TarT_v = $TarT_v, TarT_h = $TarT_h\n"

        set Tar_axis_n $TarT_axis_n
        set Tar_axis_v $TarT_axis_v
        set Tar_axis_h $TarT_axis_h

        puts "\nReference axes coordinates : $Ref"
        puts "\nTarget axes coordinates : $Tar"
        
        
        set Ref_axes_nvh_model [list $Ref_axis_n $Ref_axis_v $Ref_axis_h]
        puts $PrAxesFile "\nReference ($fixedDomain_name) Principal Axes\n"
        puts $PrAxesFile $Ref_axes_nvh_model
        puts $PrAxesFile "\n"

        set Tar_axes_nvh_model [list $Tar_axis_n $Tar_axis_v $Tar_axis_h]
        puts $PrAxesFile "\nTarget ($domain_name) Principal Axes\n"
        puts $PrAxesFile $Tar_axes_nvh_model
        puts $PrAxesFile "\n"


        puts "\ndomain name: $domain_name"
        set rvar Ref_axes_nvh_all_model_$domain_name
        set tvar Tar_axes_nvh_all_model_$domain_name

        # Store the Pairwise Axes for computing the domain rotations 
        lappend $rvar $Ref_axes_nvh_model
        lappend $tvar $Tar_axes_nvh_model
            

    }

}

#######################################################################################################################
### use the center position to label the display with current pdb file name
set fixedmol_cen [measure center $fixedMol weight mass] 
# draw a red sphere at the origin (0,0,0) , which will coincide with the molecule being fixed or aligned to the XYZ coordinate axes
draw color red
draw sphere $fixedmol_cen radius 1.0 


display update on


# optionally delete the variables to free up space
for {set dn 0} {$dn < $num_domains} {incr dn} {
    
    [lindex $domain_sels_mol $dn] delete

}    
