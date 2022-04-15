###########################################################################################
# Polar angles Theta and Psi directly from the principal axes coordinates
###########################################################################################

 proc TensorPolarAngles {I domname TensorAnglesFile_num} {

      variable Pi
      set theta_axis3 [expr [expr acos([lindex $I {2 2}])]*180/$Pi]
      set psi_axis3 [expr [expr atan2([lindex $I {2 1}],[lindex $I {2 0}])]*180/$Pi]


      set theta_axis2 [expr [expr acos([lindex $I {1 2}])]*180/$Pi]
      set psi_axis2 [expr [expr atan2([lindex $I {1 1}],[lindex $I {1 0}])]*180/$Pi]


      set theta_axis1 [expr [expr acos([lindex $I {0 2}])]*180/$Pi]
      set psi_axis1 [expr [expr atan2([lindex $I {0 1}],[lindex $I {0 0}])]*180/$Pi]

      set polarAngles [list [list $theta_axis3 $psi_axis3] [list $theta_axis2 $psi_axis2] [list $theta_axis1 $psi_axis1]]

      puts "Theta for $domname axis3: $theta_axis3"
      puts "Psi for $domname axis3: $psi_axis3\n"

      puts "Theta for $domname axis2: $theta_axis2"
      puts "Psi for $domname axis2: $psi_axis2\n"

      puts "Theta for $domname axis1: $theta_axis1"
      puts "Psi for $domname axis1: $psi_axis1\n"
      puts "\n\n"
      puts $TensorAnglesFile_num $polarAngles 
      return $polarAngles
  }
