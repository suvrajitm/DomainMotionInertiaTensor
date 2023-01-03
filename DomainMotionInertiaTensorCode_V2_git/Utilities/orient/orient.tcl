package provide Orient 1.0
package require La

namespace eval ::Orient:: {
    namespace export orient
}

# Find a proper solution to consistent axis direction, currently it works partially
# last modified Nov 03, 2015, Suvrajit Maji, sm4073@cumc.columbia.edu, suvrajit@gmail.com
# modified to keep the orientation of the principal axes in x,y,x consistent 
# I have tried to keep the convention that the principal axes always points towards the positive direction for each axis, if needed
# it still does not give correct result everytime due to the variation of mass
 
# Example
# package require Orient
# namespace import Orient::orient
# ... load your molecules and make a selection ...
#
# set I [draw principalaxes $sel]           <--- show/calc the principal axes
# set A [orient $sel [lindex $I 2] {0 0 1}] <--- rotate axis 2 to match Z
# $sel move $A
# set I [draw principalaxes $sel]           <--- recalc principal axes to check
# set A [orient $sel [lindex $I 1] {0 1 0}] <--- rotate axis 1 to match Y
# $sel move $A
# set I [draw principalaxes $sel]           <--- recalc principal axes to check#


proc Orient::sel_com { sel weights } {
    set x [ $sel get x ]
    set y [ $sel get y ]
    set z [ $sel get z ]
    set m $weights
    
    set comx 0
    set comy 0
    set comz 0
    set totalm 0
    foreach xx $x yy $y zz $z mm $m {
        # use the abs of the weights
        set mm [expr abs($mm)]
	set comx [ expr "$comx + $xx*$mm" ]
	set comy [ expr "$comy + $yy*$mm" ]
	set comz [ expr "$comz + $zz*$mm" ]
	set totalm [ expr "$totalm + $mm" ]
    }
    set comx [ expr "$comx / $totalm" ]
    set comy [ expr "$comy / $totalm" ]
    set comz [ expr "$comz / $totalm" ]
    puts "Total weight: $totalm"
    return [list $comx $comy $comz]
}

proc Orient::sel_it { sel COM weights} {
    set x [ $sel get x ]
    set y [ $sel get y ]
    set z [ $sel get z ]
    set m $weights

    # compute I
    set Ixx 0
    set Ixy 0
    set Ixz 0
    set Iyy 0
    set Iyz 0
    set Izz 0
    foreach xx $x yy $y zz $z mm $m {
        # use the abs of the weights
        set mm [expr abs($mm)]
        
        # subtract the COM
        set xx [expr $xx - [lindex $COM 0]]
        set yy [expr $yy - [lindex $COM 1]]
        set zz [expr $zz - [lindex $COM 2]]

        set rr [expr $xx + $yy + $zz]

        set Ixx [expr $Ixx + $mm*($yy*$yy+$zz*$zz)]
        set Ixy [expr $Ixy - $mm*($xx*$yy)]
        set Ixz [expr $Ixz - $mm*($xx*$zz)]
        set Iyy [expr $Iyy + $mm*($xx*$xx+$zz*$zz)]
        set Iyz [expr $Iyz - $mm*($yy*$zz)]
        set Izz [expr $Izz + $mm*($xx*$xx+$yy*$yy)]

    }
    
    return [list 2 3 3 $Ixx $Ixy $Ixz $Ixy $Iyy $Iyz $Ixz $Iyz $Izz]
}

proc vmd_draw_arrow {mol start end} {
    set scaling [expr [veclength [vecsub $end $start]]/100]
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius [expr 2*$scaling]
    puts [list cone $middle $end radius [expr 5*$scaling]]
    graphics $mol cone $middle $end radius [expr 5*$scaling]
}

proc vmd_draw_vector { mol pos val } {
    set end   [ vecadd $pos [ vecscale +1 $val ] ]
    vmd_draw_arrow $mol $pos $end
}    

# find the max of some numbers
proc Orient::max { args } {
    set maxval [lindex $args 0]
    foreach arg $args {
        if { $arg > $maxval } {
            set maxval $arg
        }
    }
    return $maxval
}

# draws the three principal axes
proc vmd_draw_principalaxes { mol sel axiscol axislblcol adjustdirection man_flip_a1a2a3 adjustscale {weights domass} } {
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    set I [Orient::calc_principalaxes $sel $adjustdirection $man_flip_a1a2a3 $weights ]
    set a1 [lindex $I 0]
    set a2 [lindex $I 1]
    set a3 [lindex $I 2]

    # find the size of the system
    set minmax [measure minmax $sel]
    set ranges [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    set scale [expr .7*[Orient::max [lindex $ranges 0] \
                             [lindex $ranges 1] \
                             [lindex $ranges 2]]]
    if {$adjustscale > 0} {
      set scale $adjustscale
    }
    set scale2 [expr 1.05 * $scale]

    puts "Drawing the principal components..."
    # draw some nice vectors
    #graphics $mol delete all

    graphics $mol color $axiscol
    set COM [Orient::sel_com $sel $weights]
    vmd_draw_vector $mol $COM [vecscale $scale $a1]
    vmd_draw_vector $mol $COM [vecscale $scale $a2]
    vmd_draw_vector $mol $COM [vecscale $scale $a3]

    graphics $mol color $axislblcol
    #draw material Transparent
    graphics $mol text [vecadd $COM [vecscale $scale2 $a1]] "1" size 1 thickness 2
    graphics $mol text [vecadd $COM [vecscale $scale2 $a2]] "2" size 1 thickness 2
    graphics $mol text [vecadd $COM [vecscale $scale2 $a3]] "3" size 1 thickness 2
    
    return [list $a1 $a2 $a3]
}

# returns the three principal axes
proc Orient::calc_principalaxes { sel adjustdirection man_flip_a1a2a3 {weights domass} } {
    puts "\nCalculating principal axes."
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    puts "Getting the center-of-mass..."
    # get the COM
    set COM [Orient::sel_com $sel $weights]
    puts "Computing the inertia tensor..."
    # get the I
    set I [Orient::sel_it $sel $COM $weights]
    La::mevsvd_br I evals
     
    set a1x [lindex $I 3]
    set a1y [lindex $I 6]
    set a1z [lindex $I 9]

    set a2x [lindex $I 4]
    set a2y [lindex $I 7]
    set a2z [lindex $I 10]

    set a3x [lindex $I 5]
    set a3y [lindex $I 8]
    set a3z [lindex $I 11]

    set s1 [lindex $evals 3]
    set s2 [lindex $evals 4]
    set s3 [lindex $evals 5]
    
    puts "Inertia Tensor Eigenvector, Eigenvalue: ($a1x $a1y $a1z) , $s1"
    puts "Inertia Tensor Eigenvector, Eigenvalue: ($a2x $a2y $a2z) , $s2"
    puts "Inertia Tensor Eigenvector, Eigenvalue: ($a3x $a3y $a3z) , $s3"

    
    
    

    #based on max value of all axis for each of x,y,z coordinate values
    # make some adjustments to the directions if necessary
    if {$adjustdirection > 0} {
	set max_coordinate_axis 0
	set max_axis_coordinate 1
	set axis_coordinate_sign 0
	set axisprojectionval 0.70

        
	if {$max_coordinate_axis == 1} {

	      #find the max absolute values of the x,y,z - coords

	      set maxx [expr max (abs($a1x), abs($a2x), abs($a3x))] 
	      puts "\na1x = $a1x"  
	      puts "a2x = $a2x"  
	      puts "a3x = $a3x\n"  

	      if {(([expr abs($a1x) == $maxx]) || ([expr abs($a1x)] > $axisprojectionval)) && $a1x < 0}  {
		    
		    set a1x [expr -$a1x]
		    set a1y [expr -$a1y]
		    set a1z [expr -$a1z]
		    puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive X direction.\n"

	      } 
	      if {(([expr abs($a2x) == $maxx]) || ([expr abs($a2x)] > $axisprojectionval))  && $a2x < 0} {

		    set a2x [expr -$a2x]
		    set a2y [expr -$a2y]
		    set a2z [expr -$a2z]
		    puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive X direction.\n"

	      }
	      if {(([expr abs($a3x) == $maxx]) || ([expr abs($a3x)] > $axisprojectionval)) && $a3x < 0} {

		    set a3x [expr -$a3x]
		    set a3y [expr -$a3y]
		    set a3z [expr -$a3z]
		    puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive X direction.\n"
	      }  


	      
	      set maxy [expr max (abs($a1y), abs($a2y), abs($a3y))] 
	      puts "\na1y = $a1y"  
	      puts "a2y = $a2y"  
	      puts "a3y = $a3y\n"  

	      if {(([expr abs($a1y) == $maxy]) || ([expr abs($a1y)] > $axisprojectionval)) && $a1y < 0}  {
		    
		    set a1x [expr -$a1x]
		    set a1y [expr -$a1y]
		    set a1z [expr -$a1z]
		    puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive Y direction.\n"

	      } 
	      if {(([expr abs($a2y) == $maxy]) || ([expr abs($a2y)] > $axisprojectionval))  && $a2y < 0} {

		    set a2x [expr -$a2x]
		    set a2y [expr -$a2y]
		    set a2z [expr -$a2z]
		    puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive Y direction.\n"

	      }
	      if {(([expr abs($a3y) == $maxy]) || ([expr abs($a3y)] > $axisprojectionval)) && $a3y < 0} {

		    set a3x [expr -$a3x]
		    set a3y [expr -$a3y]
		    set a3z [expr -$a3z]
		    puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive Y direction.\n"
	      }  



	      set maxz [expr max (abs($a1z), abs($a2z), abs($a3z))] 
	      puts "\na1z = $a1z"  
	      puts "a2z = $a2z"  
	      puts "a3z = $a3z\n"  

	      if {(([expr abs($a1z) == $maxz]) || ([expr abs($a1z)] > $axisprojectionval)) && $a1z < 0}  {
		    
		    set a1x [expr -$a1x]
		    set a1y [expr -$a1y]
		    set a1z [expr -$a1z]
		    puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive Z direction.\n"

	      } 
	      if {(([expr abs($a2z) == $maxz]) || ([expr abs($a2z)] > $axisprojectionval)) && $a2z < 0} {

		    set a2x [expr -$a2x]
		    set a2y [expr -$a2y]
		    set a2z [expr -$a2z]
		    puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive Z direction.\n"

	      }
	      if {(([expr abs($a3z) == $maxz]) || ([expr abs($a3z)] > $axisprojectionval)) && $a3z < 0} {

		    set a3x [expr -$a3x]
		    set a3y [expr -$a3y]
		    set a3z [expr -$a3z]
		    puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive Z direction.\n"
	      }  
	
	}




	#### based on max components of each axis
	if { $max_axis_coordinate == 1} {
	    set maxa1 [expr max (abs($a1x), abs($a1y), abs($a1z))] 
	    puts "\na1x = $a1x"  
	    puts "a1y = $a1y"  
	    puts "a1z = $a1z\n"  
		
	    if { (([expr abs($a1x) == $maxa1]) || ([expr abs($a1x)] > $axisprojectionval)) && $a1x < 0}  {
		  
		  set a1x [expr -$a1x]
		  set a1y [expr -$a1y]
		  set a1z [expr -$a1z]
		  puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive X direction.\n"

	    } 
	    if { (([expr abs($a1y) == $maxa1]) || ([expr abs($a1y)] > $axisprojectionval)) && $a1y < 0} {

		  set a1x [expr -$a1x]
		  set a1y [expr -$a1y]
		  set a1z [expr -$a1z]
		  puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive Y direction.\n"

	    }
	    if { (([expr abs($a1z) == $maxa1]) || ([expr abs($a1z)] > $axisprojectionval)) && $a1z < 0} {

		  set a1x [expr -$a1x]
		  set a1y [expr -$a1y]
		  set a1z [expr -$a1z]
		  puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive Z direction.\n"
	    }  


	    set maxa2 [expr max (abs($a2x), abs($a2y), abs($a2z))] 
	    puts "\na2x = $a2x"  
	    puts "a2y = $a2y"  
	    puts "a2z = $a2z\n" 

	    if { (([expr abs($a2x) == $maxa2]) || ([expr abs($a2x)] > $axisprojectionval)) && $a2x < 0}  {
		  
		  set a2x [expr -$a2x]
		  set a2y [expr -$a2y]
		  set a2z [expr -$a2z]
		  puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive X direction.\n"

	    } 
	    if { (([expr abs($a2y) == $maxa2]) || ([expr abs($a2y)] > $axisprojectionval)) && $a2y < 0} {

		  set a2x [expr -$a2x]
		  set a2y [expr -$a2y]
		  set a2z [expr -$a2z]
		  puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive Y direction.\n"

	    }
	    if { (([expr abs($a2z) == $maxa2]) || ([expr abs($a2z)] > $axisprojectionval)) && $a2z < 0} {

		  set a2x [expr -$a2x]
		  set a2y [expr -$a2y]
		  set a2z [expr -$a2z]
		  puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive Z direction.\n"
	    }  


	    set maxa3 [expr max (abs($a3x), abs($a3y), abs($a3z))] 
	    puts "\na3x = $a3x"  
	    puts "a3y = $a3y"  
	    puts "a3z = $a3z\n" 

	    if { (([expr abs($a3x) == $maxa3]) || ([expr abs($a3x)] > $axisprojectionval)) && $a3x < 0}  {
		  
		  set a3x [expr -$a3x]
		  set a3y [expr -$a3y]
		  set a3z [expr -$a3z]
		  puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive X direction.\n"

	    } 
	    if { (([expr abs($a3y) == $maxa3]) || ([expr abs($a3y)] > $axisprojectionval)) && $a3y < 0} {

		  set a3x [expr -$a3x]
		  set a3y [expr -$a3y]
		  set a3z [expr -$a3z]
		  puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive Y direction.\n"

	    }
	    if { (([expr abs($a3z) == $maxa3]) || ([expr abs($a3z)] > $axisprojectionval)) && $a3z < 0} {

		  set a3x [expr -$a3x]
		  set a3y [expr -$a3y]
		  set a3z [expr -$a3z]
		  puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive Z direction.\n"
	    }  
	}
	

	proc flipSignv1 {v1 vx vy vz} {
	   if {$v1 < 0.0} {
		set vsgn -1
	   } else {
	        set vsgn 1 
	   }
	   set vx [vecscale $vx $vsgn]
	   set vy [vecscale $vy $vsgn]
	   set vz [vecscale $vz $vsgn]
			
	   return [list $vx $vy $vz $vsgn]
	}
	
	#### based on max components of each axis
	if { $axis_coordinate_sign == 1} {
	    #set maxa1 [expr max (abs($a1x), abs($a1y), abs($a1z))] 
	    puts "\na1x = $a1x"  
	    puts "a1y = $a1y"  
	    puts "a1z = $a1z\n"  
            puts "\nAttempting sign flip a1 ..."
            
	    lassign [flipSignv1 $a1x $a1x $a1y $a1z] a1x a1y a1z sgna1
	    if {$sgna1 < 0} {
		puts "\nNote: Principal axis 1 has been inverted to keep the orientation consistent with positive X direction.\n"
	    }
	  
	    #set maxa2 [expr max (abs($a2x), abs($a2y), abs($a2z))] 
	    puts "\na2x = $a2x"  
	    puts "a2y = $a2y"  
	    puts "a2z = $a2z\n" 
	    puts "\nAttempting sign flip a2 ..."
		    
	    lassign [flipSignv1 $a2x $a2x $a2y $a2z] a2x a2y a2z sgna2
	    if {$sgna2 < 0} {
		puts "\nNote: Principal axis 2 has been inverted to keep the orientation consistent with positive X direction.\n"
	    }	  

	    

	    #set maxa3 [expr max (abs($a3x), abs($a3y), abs($a3z))] 
	    puts "\na3x = $a3x"  
	    puts "a3y = $a3y"  
	    puts "a3z = $a3z\n" 
            puts "Attempting sign flips a3 ..."	
		  
	    lassign [flipSignv1 $a3x $a3x $a3y $a3z] a3x a3y a3z sgna3
	    if {$sgna3 < 0} {
   		puts "\nNote: Principal axis 3 has been inverted to keep the orientation consistent with positive X direction.\n"
   	    } 
       }
    
   }
    
    

    
    # now $I holds in its columns the principal axes
    #set a1 "[lindex $I 3] [lindex $I 6] [lindex $I 9]"
    #set a2 "[lindex $I 4] [lindex $I 7] [lindex $I 10]"
    #set a3 "[lindex $I 5] [lindex $I 8] [lindex $I 11]"

    set a1 "$a1x $a1y $a1z"    
    set a2 "$a2x $a2y $a2z"
    set a3 "$a3x $a3y $a3z"
    
    set Axes [list $a1 $a2 $a3]
      
    lassign $man_flip_a1a2a3 flipA1 flipA2 flipA3
  
    
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
         
    set Axes [FlipAxes $Axes $flipA1 $flipA2 $flipA3]
	
	
    return $Axes
    #return [list $a1 $a2 $a3]
}


# rotate a selection about its COM, taking <vector1> to <vector2>
# e.g.: orient $sel [lindex $I 2] {0 0 1}
# (this aligns the third principal axis with z)
proc Orient::orient { sel vector1 vector2 {weights domass}} {
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }
    set COM [Orient::sel_com $sel $weights]
    set vec1 [vecnorm $vector1]
    set vec2 [vecnorm $vector2]
    # compute the angle and axis of rotation
    set rotvec [veccross $vec1 $vec2]
    set sine   [veclength $rotvec]
    set cosine [vecdot $vec1 $vec2]
    set angle [expr atan2($sine,$cosine)]

    # return the rotation matrix
    return [trans center $COM axis $rotvec $angle rad]
    puts $trans 
    puts $center
    puts $COM
    puts $axis
    puts $rotvec
    puts $angle
    puts $rad 
}
