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

  set imgfiles [glob -nocomplain $resultsrootdir /*.$figext]
  set rawimgfiles [glob -nocomplain $resultsrootdir /*.$renderfigext]
  set giffiles [glob -nocomplain $resultsrootdir /*.gif]
  set datfiles [glob -nocomplain $resultsrootdir /*.dat]

  set nimgfiles [llength $imgfiles]
  set nrawimgfiles [llength $rawimgfiles]
  set ngiffiles [llength $giffiles]
  set ndatfiles [llength $datfiles]


  if {$nimgfiles > 1 || $ngiffiles > 1 || $ndatfiles > 1} {

      puts "There are incomplete set of files (.dat or .rgb/.bmp/.jpg/.tga or .gif) from previous run in the working directory.\nDeleting all of them before starting a new run"
      
      for {set f 0} {$f < $nimgfiles} {incr f} {
	  file delete [lindex $imgfiles $f]
	  file delete [lindex $rawimgfiles $f]
      }

      for {set f 0} {$f < $ngiffiles} {incr f} {  
	  file delete [lindex $giffiles $f]
      }

      for {set f 0} {$f < $ndatfiles} {incr f} { 
	  file delete [lindex $datfiles $f]
      }
  }
  
  mol delete all

  # How to clear all variabless from previous run ?? 
  # should use it to make sure every variable is reset/unset
  # this does not remove all user created variables
  foreach var [info local] {
      unset $var
  }
  foreach var [info local] {
      unset $var
  }
  #unset -nocomplain var
  eval unset -nocomplain --[info local]
  # clear screen
  clear  