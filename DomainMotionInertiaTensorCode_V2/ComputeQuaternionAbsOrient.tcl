 ##############################################################################################################
  ## Script for Quantitative Characterization of Domain Motions in Molecular Machines
  ##
  ## This part calls the absolute orientation procedure to calculate the transformation of the coordinate axes
  ## 
  ##
  ## Suvrajit Maji,sm4073@cumc.columbia.edu 
  ## Columbia University
  ## Created: Nov 2015. Modified:Oct 13,2016
  ##############################################################################################################


proc ComputeQuaternionAbsOrient {dom_axes_nvh_all_model nfl initState_fileno} {

    set StateA_fileno $initState_fileno 
    set StateB_fileno $nfl

    set dom_axes_nvh_model_StateA [lindex $dom_axes_nvh_all_model $StateA_fileno]
    set dom_axes_nvh_model_StateB [lindex $dom_axes_nvh_all_model $StateB_fileno]
	
    set StateA_Axes [concat [lindex $dom_axes_nvh_model_StateA 0] [lindex $dom_axes_nvh_model_StateA 1] [lindex $dom_axes_nvh_model_StateA 2]]
    set StateB_Axes [concat [lindex $dom_axes_nvh_model_StateB 0] [lindex $dom_axes_nvh_model_StateB 1] [lindex $dom_axes_nvh_model_StateB 2]]


    puts "State A axes coords : $StateA_Axes"
    puts "State B axes coords : $StateB_Axes"
    
    # StateA --> StateB : this will give us the rotation axis and direction in which the model has rotated 
    # from the state A to state B (transformed)
    set absorientQAxisAngle_dom  [AbsoluteOrientation $StateA_Axes $StateB_Axes]

    # StateB --> StateA should provide us the same rotation angle but opposite direction
    #set absorientQAxisAngle_dom  [AbsoluteOrientation $StateB_Axes $StateA_Axes]


    return $absorientQAxisAngle_dom


} 