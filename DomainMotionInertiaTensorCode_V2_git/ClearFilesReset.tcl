##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
##
## Clear residual files from previous execution and reset created varibles on VMD
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: May 2015. Modified:Oct 13,2016
##############################################################################################################


###########################################################################################
# delete old files
###########################################################################################
# use the delete files option with caution
# find a way to delete all the files using wild card search 
# for some reason if I do : file delete [glob -nocomplain *.dat] 
# it says long file name as the all the file names are concatenated into one file name 
# Also , file delete {*}[glob -nocomplain *.dat] , does not work 

set resultsrootdir "../DMResults"

set imgfiles [glob -nocomplain $resultsrootdir/tmp /*.$figext]
set rawimgfiles [glob -nocomplain $resultsrootdir/tmp /*.$renderfigext]
set giffiles [glob -nocomplain $resultsrootdir/tmp /*.gif]
set datfiles [glob -nocomplain $resultsrootdir/tmp /*.dat]

set nimgfiles [llength $imgfiles]
set nrawimgfiles [llength $rawimgfiles]
set ngiffiles [llength $giffiles]
set ndatfiles [llength $datfiles]

# delete the tmp results folder
exec rm -vrf $resultsrootdir/tmp

  
mol delete all

# How to clear all variables from previous run ?? 
# should use it to make sure every variable is reset/unset
# this does not remove all user created variables
foreach var [info locals] {
    unset $var
}
#unset -nocomplain var
eval unset -nocomplain --[info locals]
# clear screen
clear  
