
##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
## 
## This part summarizes the motion characterization results
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Dec 2022. Modified:Dec 28,2022 
## 
##
## ##############################################################################################################

###########################################################################################
## If one or more of the Principal axes gets flipped in direction from one model to the 
## next, then we can manually revert the sign of the flipped axes
###########################################################################################

proc FlipAxes {Axes flipA1 flipA2 flipA3} {
     set PA1 [lindex $Axes 0]
     set PA2 [lindex $Axes 1]
     set PA3 [lindex $Axes 2]

     set PA1 [vecscale $flipA1 $PA1]
     set PA2 [vecscale $flipA2 $PA2]
     set PA3 [vecscale $flipA3 $PA3]

     set Axes [list $PA1 $PA2 $PA3]
     return $Axes
}

