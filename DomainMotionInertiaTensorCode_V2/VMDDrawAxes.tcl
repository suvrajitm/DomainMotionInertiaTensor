   ####################### 
   #### VMD Draw Axes
   #######################  

proc VMDDrawAxes {COM} {
    
    ### X,Y,Z-axis 
    set a1 {1 0 0}
    set a2 {0 1 0}
    set a3 {0 0 1}
    #set COM {0 0 0}

    set axis1col "red"
    set axis2col "green"
    set axis3col "blue"
    set axislblcol "black"

    set scale1 50
    set scale2 [expr 1.02*$scale1]

    graphics top color $axis1col
    vmd_draw_vector top $COM [vecscale $scale1 $a1]
    graphics top color $axis2col
    vmd_draw_vector top $COM [vecscale $scale1 $a2]
    graphics top color $axis3col
    vmd_draw_vector top $COM [vecscale $scale1 $a3]

    graphics top color $axislblcol
    #draw material Transparent
    graphics top text [vecadd $COM [vecscale $scale2 $a1]] "X" size 2 thickness 3
    graphics top text [vecadd $COM [vecscale $scale2 $a2]] "Y" size 2 thickness 3
    graphics top text [vecadd $COM [vecscale $scale2 $a3]] "Z" size 2 thickness 3

}

    
proc vmd_draw_arrow {mol start end} {
      set scaling [expr [veclength [vecsub $end $start]]/100]
      # an arrow is made of a cylinder and a cone
      set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
      graphics $mol cylinder $start $middle radius [expr 1.5*2*$scaling]
      puts [list cone $middle $end radius [expr 5*$scaling]]
      graphics $mol cone $middle $end radius [expr 5*2.5*$scaling]
}

proc vmd_draw_vector { mol pos val } {
	set end   [ vecadd $pos [ vecscale +1 $val ] ]
	vmd_draw_arrow $mol $pos $end
}  