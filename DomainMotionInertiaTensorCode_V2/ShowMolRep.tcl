##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
##
## This part toggles off/on for the set of representations in a model
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Aug 26, 2016. Modified:Aug 26, 2016
##############################################################################################################

proc ShowMolRep {totalrepno repstatus} {

  if {$repstatus=="off"} {

      for {set i 0} { $i < $totalrepno } {incr i} {
	  mol showrep top $i off
      }

  } else {

      for {set i 0} { $i < $totalrepno } {incr i} {
	  mol showrep top $i on
      }
  }
}