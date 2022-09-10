##############################################################################################################
## Script for Quantitative Characterization of Domain Motions in Molecular Machines
## 
## This part is the tcl implementation of the paper (B.K.P.Horn 1986/87) which solves the least-squares 
## absolute orientation problem based on unit quaternions given a pair of coordinate axes system 
## (Principal axes in our case) A & B . Note that this is just the basic implementation of the actual method
## Future plan is to include other methods in here e.g. Umeyama's method with rotation matrix etc.
## 
## Also implemented are few operations for quaternion operatios, 
## including Quaternion decomposition with Clifford algebra and vector projection
## The implemented procedures in this script is the core of the motion characterization method along with the 
## inertia tensor and principal axes calculation
##
## Suvrajit Maji,sm4073@cumc.columbia.edu 
## Columbia University
## Created: Nov 2015. Modified:October 25,2016 
## 
##
##############################################################################################################

lappend auto_path Utilities/la1.0
lappend auto_path Utilities/orient

#source Utilities/la1.0/la.tcl

package require La
catch { namespace import La::*}
#namespace import ::math::complexnumbers::*

proc AbsoluteOrientation {A B} {
     
    puts "\n\nHorn's (1987) Closed Form Solution to the Least Square problem of Absolute Orientation Using Unit Quaternions"
    puts "\nMatrices of object and reference coordinates"
    
    # make the array compatible with the La package
    set A [concat 2 3 3 $A]
    set B [concat 2 3 3 $B]
    # The input vectors are usually row vectors in tcl 
    # Make the vectors as columns in the matrix, since the eigensolvers produces column vectors, 
    # so it will be easier to keep the same convention, to read the matrix as column vectors 
    
    set A [La::transpose $A]
    set B [La::transpose $B]
  
    ## can also use A & -A , B & -B to make the centroid at origin, remove the comments below 
    #set A [La::join_cols $A [La::mscale $A -1]]
    #set B [La::join_cols $B [La::mscale $B -1]]

    # 3 x N matrix (N points)
    puts "\nState A (axes domain model in frame k):"
    puts [La::show $A]

    # 3 x N matrix (N points)
    puts "\nState B (axes domain model in frame k+1):"
    puts [La::show $B]

    set dim 3
    # -3 because first 3 numbers {2 3 3} are the information needed 
    #  for constructing the array  in La package
    set numA [La::mcols $A]
    set numB [La::mcols $B]
    puts "Number of vectors (A) = $numA"
    puts "Number of vectors (B) = $numB"
    if {$numA!= $numB} {
	puts "Matrix dimensions doesn't match."
	return
    } else {
	set numVecs $numA
    }
    
    puts "\nCompute centroids as column means (on transpose)"
    set Atr [La::transpose $A]
    set Btr [La::transpose $B]
    La::mnorms_br Atr centroidL centroidL_sig
    La::mnorms_br Btr centroidR centroidR_sig

    #set centroidL [La::transpose $centroidL]
    #set centroidL [La::transpose $centroidR]

    puts "\nCentroid left = [show $centroidL]"
    puts "Centroid right = [show $centroidR]"

    #puts "\nCentroid std dev left = [show $centroidL_sig]"
    #puts "Centroid std dev right = [show $centroidR_sig]"
   
    # any easy way to do matlab's "repmat" like operation ? to compose the centroid matrix for subtraction ? 
    puts "\nCompose the centroid matrix for subtraction:"
    if { $numVecs == 3 } {
	set matcenL [La::join_cols $centroidL [La::join_cols $centroidL $centroidL]]
	set matcenR [La::join_cols $centroidR [La::join_cols $centroidR $centroidR]]

    } elseif {$numVecs == 6 } {
	set matcenL [La::join_cols $centroidL [La::join_cols $centroidL $centroidL]]
	set matcenL [La::join_cols $matcenL $matcenL]

	set matcenR [La::join_cols $centroidR [La::join_cols $centroidR $centroidR]]
	set matcenR [La::join_cols $matcenR $matcenR]
    }

    
    #puts "\nMatrix of centroid left:" 
    puts [La::show $matcenL]
    #puts "Matrix of centroid right:" 
    puts [La::show $matcenR]

    set LeftM [La::msub $A $matcenL]
    set RightM [La::msub $B $matcenR]
  

    # 3 x N matrix 
    puts "\nState A matrix LeftM with centroid subtracted:"
    puts [La::show $LeftM]

    # 3 x N matrix 
    puts "\nState B matrix RightM with centroid subtracted:"
    puts [La::show $RightM]


    puts "\n\nComputing the matrix M "
    # (3 x N matrix)  x  (3 x N)^transpose  matrix
    # = 3 x N * N x 3 matrix = 3 x 3 matrix 
    set M [La::mmult $LeftM [La::transpose $RightM]]
    puts [La::show $M]

    puts "\nGetting the matrix elements of M: Sxx, Sxy, Sxz, etc ...\n"
    set Sxx [lindex $M 3]; set Sxy [lindex $M 4]; set Sxz [lindex $M 5];
    set Syx [lindex $M 6]; set Syy [lindex $M 7]; set Syz [lindex $M 8];
    set Szx [lindex $M 9]; set Szy [lindex $M 10]; set Szz [lindex $M 11];

    puts "\nCreating the matrix elements of N\n"
    set N11 [expr ($Sxx + $Syy + $Szz)];  set N12 [expr ($Syz - $Szy)];          set N13 [expr ($Szx - $Sxz)]; 		set N14 [expr ($Sxy - $Syx)];
    set N21 [expr ($Syz - $Szy)];         set N22 [expr ($Sxx - $Syy - $Szz)]; 	 set N23 [expr ($Sxy + $Syx)]; 		set N24 [expr ($Szx + $Sxz)];	
    set N31 [expr ($Szx - $Sxz)];         set N32 [expr ($Sxy + $Syx)];          set N33 [expr (-$Sxx + $Syy - $Szz)]; 	set N34 [expr ($Syz + $Szy)];
    set N41 [expr ($Sxy - $Syx)];         set N42 [expr ($Szx + $Sxz)];          set N43 [expr ($Syz + $Szy)];          set N44 [expr (-$Sxx - $Syy + $Szz)]; 	

    puts "N = "
    set N [list 2 4 4 $N11 $N12 $N13 $N14 $N21 $N22 $N23 $N24 $N31 $N32 $N33 $N34 $N41 $N42 $N43 $N44]
    puts [La::show $N]

    set V $N


    # It appears that the eigensolver in the 'La package' is not very accurate (at least from what I have observed), 
    # Future work : try to use a different solver, preferably which can also handle complex numbers
    # In real applications, the eigenvalues/eigenvectors may have complex components due to 
    # numerical errors and hence should be handled accordingly
    # For now the eigenvectors corresponding to the max and sometimes the min (if unique) eigenvalues seems 
    # to be consistent with other solvers, check ?

    puts "\n\nPerforming the eigen decomposition of N"
    La::mevsvd_br V evals
    puts "\nEigenvectors V = "    
    puts [La::show $V]

    # Is there a easier way to access the rows and columns of a matrix in La ? 
    # For now get the rows or columns by direct indexing 
    set V1w [lindex $V 3]
    set V1x [lindex $V 7]
    set V1y [lindex $V 11]
    set V1z [lindex $V 15]
    set eV1 [list $V1w $V1x $V1y $V1z]

    set V2w [lindex $V 4]
    set V2x [lindex $V 8]
    set V2y [lindex $V 12]
    set V2z [lindex $V 16]
    set eV2 [list $V2w $V2x $V2y $V2z]

    set V3w [lindex $V 5]
    set V3x [lindex $V 9]
    set V3y [lindex $V 13]
    set V3z [lindex $V 17]
    set eV3 [list $V3w $V3x $V3y $V3z]

    set V4w [lindex $V 6]
    set V4x [lindex $V 10]
    set V4y [lindex $V 14]
    set V4z [lindex $V 18]
    set eV4 [list $V4w $V4x $V4y $V4z]
    

    set eigV [list $eV1 $eV2 $eV3 $eV4]
  

    puts "\nEigenvalues D ="
    puts [La::show $evals]
    
    set D [lreplace $evals 0 2]
    puts "D = $D"
    set evalmax [lindex [lsort -real -decreasing $D] 0] 
    set evalmaxid [lindex [lsearch -real $D $evalmax] 0]
    puts "Maximum eigenvalue id : $evalmaxid "
  
    set q [lindex $eigV $evalmaxid]
    set w [lindex  $q 0]
    set x [lindex  $q 1]
    set y [lindex  $q 2]
    set z [lindex  $q 3]


    # resolve the sign ambiguity for the eigenvector
    # by making the sign of largest (magnitude wise) component of the eigenvector as positive
    # This would give a consistent directions for the eigenvector : +q and -q represents the same rotation 
    set absq [list [expr abs($w)] [expr abs($x)] [expr abs($y)] [expr abs($z)]]
    set absqmax [lindex [lsort -real -decreasing $absq] 0] 
    set absqmaxid [lindex [lsearch -real $absq $absqmax] 0]
    set qabsqmax [lindex $q $absqmaxid] 

    if {$qabsqmax < 0} {
	set qmaxsgn -1
    } else {
        set qmaxsgn 1 
    }
    set q [vecscale $q $qmaxsgn]
    set Q [unitQuaternion $q]
    puts "\nRotation Unit Quaternion Q = $Q"
    

    puts "\nUnit Quaternion to Rotation Matrix..."
    set RotMat [quatToRotMat $Q]

    # can also extract the centroid subtracted vectors and then use VecSumSquare
    # to get the sum of norm squares 

    # sum of the squares of the vector norms is same as the 
    # sum of the square of matrix elements constituting the vectors
    set MprodL [La::mprod $LeftM $LeftM]
    set MprodR [La::mprod $RightM [La::mmult $RotMat $LeftM]]

    # ssL,R = sum of squares on the left/Right 
    set ssL [La::msum [La::msum $MprodL]]
    set ssR [La::msum [La::msum $MprodR]]
    puts "\nLeft(A) sum of squares, ssL = [La::show $ssL]"
    puts "Right(B) sum of squares, ssR = [La::show $ssR]"

    # optimal transformation scale 
    set transformScale [expr ($ssR/$ssL)]
    puts "\nOptimal transformation (A -> B) scale: $transformScale"

    # optimal translation offset
    # t = centroidR - scale *(R * centroidL)
    set translationOffset [La::msub $centroidR  [La::mscale [La::mmult $RotMat $centroidL] $transformScale]]
    puts "Optimal transformation (A -> B) translational offset: [La::show [La::transpose $translationOffset]]"

    # any easy way to make a matlab's "repmat" like operation ?
    set translationOffsetMat [La::join_cols [La::join_cols $translationOffset $translationOffset] $translationOffset]

    puts "\nMultiplying the left coordinates (A) with rotation matrix R"

    # rotated left coordinates should be close to the corresponding right coordinates
    set RotLeftM [La::mmult [La::mscale $RotMat $transformScale] $LeftM]
    set RotLeftM [La::msub $RotLeftM $translationOffsetMat]

    puts "\nRotated left coordinates R(A):" 
    puts [La::show $RotLeftM]

    puts "Right Coordinates B:"
    puts [La::show $RightM]
    
    # compute the euclidean residual error
    set DiffMat [La::msub $RightM $RotLeftM]
    set SumSquareDiffMat [La::msum [La::msum [La::mprod  $DiffMat $DiffMat]]]
    set TransformationError [expr sqrt($SumSquareDiffMat)]
    puts "\nTransformation Error: [La::show $TransformationError]"

    # compute the rotation axis and rotation angle corresponding to the Unit Quaternion Q obtained from 
    # solving the absolute orientation problem
    puts "\n\nCalculating Angle-Axis"
    set AxisAngle [ComputeAxisAngle $Q]
   
    set rotAxis [lindex $AxisAngle 0]
    set rotAngle [expr [lindex $AxisAngle 1]]
    puts "\nAxis-Angle: ($rotAxis) $rotAngle"
    
    # Obtain Euler anlges from the unit Quaternion
    set EulA [QuatToEuler $Q]
    puts "Quaternion to Euler Angles : $EulA" 
    
    set QuatAxisAngle [list $Q $rotAxis $rotAngle $EulA]
    return $QuatAxisAngle

}


#### We can also separate out the methods on just the quaternion to a separate .tcl file 
#### and source it for AbsoluteOrientation.tcl

#### Quaternion Swing-Twist decomposition in Clifford algebra
#### Dobrowolski 2015
#### unit quaternion is a spinor in 3-D, so we replaced the algorithm for spinor 
#### with a unit quaternion 
#### spinor s = a + b*e12 + c*e23 + d*e31
#### bivectors e23,e31,e12 are respectively equal to i,j,k for quaternions
  
proc QuatDecompositionClifford {Q V} {
    
   ### Qt:component of Q around a specified non-zero vector V
   ### Q = Qs*Qt
   
   lassign $Q qw qx qy qz
   lassign $V vx vy vz

   set a $qw
   set b $qz
   set c $qx
   set d $qy

   ### Twist component:Qt
   # set u [expr $vx*$qx + $vy*$qy + $vz*$qz]
   set u [expr $vx*$c + $vy*$d + $vz*$b]

   set n [expr $vx*$vx + $vy*$vy + $vz*$vz]
   #can also just use veclength2 $V

   set m [expr $a*$n]

   set l [expr sqrt($m*$m + $u*$u*$n)]

   set qt [list [expr $m/$l] [expr ($vx*$u/$l)]  [expr ($vy*$u/$l)] [expr ($vz*$u/$l)]] 

   set Qt [unitQuaternion $qt]
   set QtAxisAngle [ComputeAxisAngle $Qt]

   ### Swing component: Qs
   set Qtinv [QuatInverse $Qt]
   set Qs [QuatMultiplication $Q $Qtinv]
   set QsAxisAngle [ComputeAxisAngle $Qs]

   return [list $Qt $QtAxisAngle $Qs $QsAxisAngle]
}


#### Quaternion Swing-Twist decomposition
#### Derived from the Clifford algebra algorithm above
#### but represents a more geometric intuition
proc QuatDecomposition {Q V} {
    
   ### Qt:component of Q around a specified non-zero vector V
   ### Q = Qs*Qt
   
   lassign $Q qw qx qy qz

   ### rotation axis 
   set ra [list $qx $qy $qz]
    
   ### Twist componen:Qt
   set radotv [vecdot $ra $V]
   set v_radotv [vecscale $radotv $V]
   set vnorm2 [veclength2 $V]
   set Vp [vecscale $v_radotv [expr 1.0/$vnorm2]]

   set qtw $qw 
   set qtx [lindex $Vp 0]
   set qty [lindex $Vp 1]
   set qtz [lindex $Vp 2]

   set qt [list $qtw $qtx $qty $qtz]
   set Qt [unitQuaternion $qt]
   set QtAxisAngle [ComputeAxisAngle $Qt]

   ### Swing component: Qs
   set Qtinv [QuatInverse $Qt]
   set Qs [QuatMultiplication $Q $Qtinv]
   set QsAxisAngle [ComputeAxisAngle $Qs]

   return [list $Qt $QtAxisAngle $Qs $QsAxisAngle]
    
}




proc quatToRotMat {Q} {
  
   puts "\nConverting the unit quaternion into an equivalent orthonormal rotation matrix..."

   #lassign $Q qw qx qy qz

   set qw [lindex  $Q 0]
   set qx [lindex  $Q 1]
   set qy [lindex  $Q 2]
   set qz [lindex  $Q 3]
   
   #set R11 [expr ($qw*$qw + $qx*$qx - $qy*$qy - $qz*$qz)]
   #set R12 [expr 2*($qx*$qy - $qw*$qz)]
   #set R13 [expr 2*($qx*$qz + $qw*$qy)]
   
   #set R21 [expr 2*($qy*$qx + $qw*$qz)] 
   #set R22 [expr ($qw*$qw - $qx*$qx + $qy*$qy - $qz*$qz)]
   #set R23 [expr 2*($qy*$qz - $qw*$qx)]
   
   #set R31 [expr 2*($qz*$qx - $qw*$qy)]
   #set R32 [expr 2*($qz*$qy + $qw*$qz)]
   #set R33 [expr ($qw*$qw - $qx*$qx - $qy*$qy + $qz*$qz)]
   
   set R11 [expr (1 - 2*($qy*$qy + $qz*$qz))]
   set R12 [expr 2*($qx*$qy - $qw*$qz)]
   set R13 [expr 2*($qx*$qz + $qw*$qy)]
   
   set R21 [expr 2*($qy*$qx + $qw*$qz)] 
   set R22 [expr (1 - 2*($qx*$qx + $qz*$qz))]
   set R23 [expr 2*($qy*$qz - $qw*$qx)]
   
   set R31 [expr 2*($qz*$qx - $qw*$qy)]
   set R32 [expr 2*($qz*$qy + $qw*$qx)]
   set R33 [expr (1 - 2*($qx*$qx + $qy*$qy))]

   #rotation matrix
   set RotMat [list 2 3 3 $R11 $R12 $R13 $R21 $R22 $R23 $R31 $R32 $R33]
   puts "\nRotation matrix created..."
   puts [La::show $RotMat]
   #set RotMat [La::transpose $RotMat]
 
   return $RotMat
}



proc VecSumSquare {vecs} {
   set sumSq 0
   foreach vec $vecs {
      set vnormSquare [veclength2 $vec]
      set sum [expr {$sumSq + $vnormSquare}]
   }
   return $sumSq
}

proc QuatNormSquare {Q} {

   #assign Q values to qw qx qy qz
   lassign $Q qw qx qy qz

   set normSquareQ [expr ($qw*$qw + $qx*$qx + $qy*$qy + $qz*$qz)]

   return $normSquareQ

}


proc QuatNorm {Q} {

   set normsqQ [QuatNormSquare $Q]
  
   set normQ [expr sqrt(normsqQ)]
   return $normQ

}


proc unitQuaternion {Q} {
   set uQ  [vecnorm $Q]	
   return $uQ
}

proc QuatConjugate {Q} {
   
   lassign $Q qw qx qy qz
    
   set Qconj [list $qw [expr -1*$qx] [expr -1*$qy] [expr -1*$qz]] 
  
   return $Qconj

}

proc QuatInverse {Q} {

  set q_normSqr [QuatNormSquare $Q]
  set q_conj [QuatConjugate $Q]

  lassign $q_conj qwc qxc qyc qzc
  
  set Qinv [list [expr $qwc/$q_normSqr] [expr $qxc/$q_normSqr] [expr $qyc/$q_normSqr] [expr $qzc/$q_normSqr]]

  return $Qinv

}

proc QuatMultiplication {Q1 Q2} {
  
  lassign $Q1 qw1 qx1 qy1 qz1
  lassign $Q2 qw2 qx2 qy2 qz2 

  set qwp [expr {$qw1*$qw2 - $qx1*$qx2 - $qy1*$qy2 - $qz1*$qz2}]
  set qxp [expr {$qw1*$qx2 + $qx1*$qw2 + $qy1*$qz2 - $qz1*$qy2}] 
  set qyp [expr {$qw1*$qy2 + $qy1*$qw2 + $qz1*$qx2 - $qx1*$qz2}] 
  set qzp [expr {$qw1*$qz2 + $qz1*$qw2 + $qx1*$qy2 - $qy1*$qx2}]

  set Qp [list $qwp $qxp $qyp $qzp]

}

proc ComputeAxisAngle {Q} {

   variable Pi
   #just to be sure, normalize to get unit quaternion first  
   set Q  [unitQuaternion $Q]	    

   #assign Q values to qw qx qy qz
   #lassign $Q qw qx qy qz

   set qw [lindex  $Q 0]
   set qx [lindex  $Q 1]
   set qy [lindex  $Q 2]
   set qz [lindex  $Q 3]
   
   #calculate the angles from the quaternion component qw  
   #set Angle [expr (2 * acos($qw))*180/$Pi]

   #more numerically stable version for the angle
   set Angle [expr (2 * atan2([veclength [list $qx $qy $qz]],$qw))*180/$Pi];
   
   set s [expr sqrt(1 - $qw*$qw)]

   set tol 0.0001

   if { $s < $tol } {
       set x $qx
       set y $qy
       set z $qz

   } else {   
       set x [expr $qx/$s]
       set y [expr $qy/$s]
       set z [expr $qz/$s]
   }
   set Axis [list $x $y $z]

   # renormalize the axis if there is still numerical error
   set axislen [veclength $Axis] 
   if { $axislen != 0} {
      set Axis [vecnorm $Axis]
   } else {
     puts "Rotation Angle : $Angle, Rotation Axis : $Axis ... implies no rotation." 
   }
  
   set AxisAngle [list $Axis $Angle]
   ### angle is already converted from radian to degree
   return $AxisAngle
    
}


proc QuatToEuler {Q} {
    
    variable Pi 

    #lassign $Q qw qx qy qz

    set qw [lindex  $Q 0]
    set qx [lindex  $Q 1]
    set qy [lindex  $Q 2]
    set qz [lindex  $Q 3]

    # this corresponds to rotation sequence of "ZYX". We can also get the values for all sequences
    # There is a longer version for calculating the Euler angles corresponding to all possible 
    # rotation sequences "ZYX", "XYZ", ...

    # the definitions may vary
    #Roll (Phi), also called Heading , rotation about X-axis   
    set Phi [expr atan2(2*($qw * $qx + $qy * $qz), ($qw * $qw - $qx * $qx - $qy * $qy + $qz * $qz))*180/$Pi]

    #Pitch (Theta), also called Attitude , rotation about Y-axis 
    set Theta [expr asin(2*($qw * $qy - $qz * $qx))*180/$Pi]
	
    # Yaw (Psi), also called Bank , rotation about Z-axis
    set Psi [expr atan2(2*($qw * $qz + $qx * $qy), ($qw * $qw + $qx * $qx - $qy * $qy - $qz * $qz))*180/$Pi]
	
    set EulA [list $Phi $Theta $Psi]
    
    return $EulA
}