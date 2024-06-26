##############################################################################################################
# Script to print text on the vmd screen to label the structure on the frame
# Not sure how to print the text at the bottom of the screen with rotated views 
# right now the position of the text is rotated along with the whole view
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Nov 2015. Modified:April 20,2017 
##############################################################################################################

proc LabelDisplay {pdbf fixedmol_cen axis offset lblcol} { 

    # label the figure with current pdb name
    switch $axis {
		"x" {
			set xp [expr $offset + [lindex $fixedmol_cen 0]]
			set textpos [vecadd $fixedmol_cen [list $xp 0 0]]
		}

		"y" {
			set yp [expr $offset + [lindex $fixedmol_cen 1]]
			set textpos [vecadd $fixedmol_cen [list 0 $yp 0]]
		}
		
		"z" {
			set zp [expr $offset + [lindex $fixedmol_cen 2]]
			set textpos [vecadd $fixedmol_cen [list 0 0 $zp]]
		}

		default {
			puts "Wrong axis specification !"
		}
    }

	# for now, add text at the position defined by below, irrespective of the above textpos locations 
	set xp [expr $offset + [lindex $fixedmol_cen 0]]
	set yp [expr $offset + [lindex $fixedmol_cen 1]]
	set zp [expr $offset + [lindex $fixedmol_cen 2]]
	set textposxyz [vecadd $fixedmol_cen [list $xp $yp $zp]]
    draw color $lblcol
    graphics top text $textposxyz "[format "%s" $pdbf ]" size 1 thickness 2
   
}