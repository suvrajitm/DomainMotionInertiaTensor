# Script to calculate Euler angles for a Local Coordinate System (LCS) when it is expressed in 
# terms of the World Coordinate System (WCS)

proc ComputeEuler {LCS} { 

  variable Pi 
  set tol 0.0001
  # Coordinate Axes for the local coordinate system
  set Xl [lindex $LCS 0]
  set Yl [lindex $LCS 1]
  set Zl [lindex $LCS 2]
  
  set Xlx [lindex $Xl 0]
  set Xly [lindex $Xl 1]
  set Xlz [lindex $Xl 2]

  set Ylx [lindex $Yl 0]
  set Yly [lindex $Yl 1]
  set Ylz [lindex $Yl 2]

  set Zlx [lindex $Zl 0]
  set Zly [lindex $Zl 1]
  set Zlz [lindex $Zl 2]

  set Zlxy [expr sqrt($Zlx*$Zlx + $Zly*$Zly)]
    
  if {$Zlxy > $tol} {

      set Alpha [expr atan2($Ylx*$Zly - $Yly*$Zlx, $Xlx*$Zly - $Xly*$Zlx)*180/$Pi]
      set Beta [expr atan2($Zlxy, $Zlz)*180/$Pi]
      set Gamma [expr -atan2(-$Zlx,$Zly)*180/$Pi]

  } else {

      set Alpha 0.0;
      
      if { $Zlz > 0} {

	  set Beta 0.0 

      } else {

	  set Beta 180
      }

      set Gamma [expr -atan2($Xly,$Xlx)]
  }

  set EulerAngles [list $Alpha $Beta $Gamma]

  return $EulerAngles
}  