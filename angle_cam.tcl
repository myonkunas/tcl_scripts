# Set of scripts to:
#  find the angle made by that vector with another vector; and
#  do this for all timesteps in a trajectory.
# Justin Gullingsrud
# A. Saladino
# M. Yonkunas
#  Fit the points x to x = ai + b, i=0...N-1, and return the value of a 
# a = 12/( (N(N^2 - 1)) ) sum[ (i-(N-1)/2) * xi]
# reference: Bevington
proc lsq { x } {
  set N [llength $x]
  set xtot 0
  set d [expr {0.5*($N-1)}]
  
  set i 0.0
  foreach elem $x {
    set xtot [expr {$xtot + ($i - $d) * $elem}]
    set i [expr {$i + 1.0}]
  }
  
# no need to normalize if all we want is the direction
#  set xtot [expr $xtot * 12 / ($N * ($N * $N - 1))]
  return $xtot
}    
#-------------------------------------------------------------------
proc sel_sel_angle { sel1 sel2 } {
  set x1 [$sel1 get x]
  set y1 [$sel1 get y]
  set z1 [$sel1 get z]

  set x2 [$sel2 get x]
  set y2 [$sel2 get y]
  set z2 [$sel2 get z]

  set xa1 [lsq $x1]
  set ya1 [lsq $y1]
  set za1 [lsq $z1]

  set xa2 [lsq $x2]
  set ya2 [lsq $y2]
  set za2 [lsq $z2]

  set anorm1 [expr sqrt($xa1 * $xa1 + $ya1 * $ya1 + $za1 * $za1)]
  set anorm2 [expr sqrt($xa2 * $xa2 + $ya2 * $ya2 + $za2 * $za2)]
   
  set xa1 [expr $xa1 / $anorm1]
  set ya1 [expr $ya1 / $anorm1]
  set za1 [expr $za1 / $anorm1]
 
  set xa2 [expr $xa2 / $anorm2]
  set ya2 [expr $ya2 / $anorm2]
  set za2 [expr $za2 / $anorm2]

  set costheta [expr $xa1 * $xa2 + $ya1 * $ya2 + $za1 * $za2]
  set angle [expr 180 * acos($costheta) / 3.14159]
# puts $angle
  return $angle
}

#----------------------------------------------------------------------- 
# Compute the angle for all frames
#-----------------------------------------------------------------------

proc sel_sel_angle_frames { helix1_a helix1_b helix2_a helix2_b filename } {

set sel1 [atomselect top "protein and resid $helix1_a to $helix1_b"]
set sel2 [atomselect top "protein and resid $helix2_a to $helix2_b"]

set N [molinfo top get numframes]
set fileId [open $filename a]
  for { set i 0 } { $i < $N } { incr i} {
    $sel1 frame $i
    $sel2 frame $i
set crap [sel_sel_angle $sel1 $sel2]
puts $fileId "$i [format "%2.3f" $crap]" 
  }
close $fileId 
}

#Compute the cross product between two sets of vectors and take thier dot product for every frame
#--------------------------------------------------------------------------
proc sel_cross_sel_angle_frames { helix1_a helix1_b helix2_a helix2_b helix3_a helix3_b helix4_a helix4_b filename } {

set sel1 [atomselect top "protein and resid $helix1_a to $helix1_b"]
set sel2 [atomselect top "protein and resid $helix2_a to $helix2_b"]
set sel3 [atomselect top "protein and resid $helix3_a to $helix3_b"]
set sel4 [atomselect top "protein and resid $helix4_a to $helix4_b"]

set N [molinfo top get numframes]
set fileId [open $filename a]
  for { set i 0 } { $i < $N } { incr i} {
    $sel1 frame $i
    $sel2 frame $i
    $sel3 frame $i
    $sel4 frame $i
set crap [sel_cross_sel_angle $sel1 $sel2 $sel3 $sel4]
puts $fileId "$i [format "%2.3f" $crap]" 
  }
close $fileId 
}


#-------------------------------------------------------------------------
proc sel_cross_sel_angle { sel1 sel2 sel3 sel4 } {
  set x1 [$sel1 get x]
  set y1 [$sel1 get y]
  set z1 [$sel1 get z]

  set x2 [$sel2 get x]
  set y2 [$sel2 get y]
  set z2 [$sel2 get z]

  set x3 [$sel3 get x]
  set y3 [$sel3 get y]
  set z3 [$sel3 get z]

  set x4 [$sel4 get x]
  set y4 [$sel4 get y]
  set z4 [$sel4 get z]

  set xa1 [lsq $x1]
  set ya1 [lsq $y1]
  set za1 [lsq $z1]

  set xa2 [lsq $x2]
  set ya2 [lsq $y2]
  set za2 [lsq $z2]

  set xa3 [lsq $x3]
  set ya3 [lsq $y3]
  set za3 [lsq $z3]

  set xa4 [lsq $x4]
  set ya4 [lsq $y4]
  set za4 [lsq $z4]

  set anorm1 [expr sqrt($xa1 * $xa1 + $ya1 * $ya1 + $za1 * $za1)]
  set anorm2 [expr sqrt($xa2 * $xa2 + $ya2 * $ya2 + $za2 * $za2)]
  set anorm3 [expr sqrt($xa3 * $xa3 + $ya3 * $ya3 + $za3 * $za3)]
  set anorm4 [expr sqrt($xa4 * $xa4 + $ya4 * $ya4 + $za4 * $za4)]

  set xa1 [expr $xa1 / $anorm1]
  set ya1 [expr $ya1 / $anorm1]
  set za1 [expr $za1 / $anorm1]
 
  set xa2 [expr $xa2 / $anorm2]
  set ya2 [expr $ya2 / $anorm2]
  set za2 [expr $za2 / $anorm2]

  set xa3 [expr $xa3 / $anorm3]
  set ya3 [expr $ya3 / $anorm3]
  set za3 [expr $za3 / $anorm3]

  set xa4 [expr $xa4 / $anorm4]
  set ya4 [expr $ya4 / $anorm4]
  set za4 [expr $za4 / $anorm4]

  set a1_cross_a2_x [expr $ya1 * $za2 - $ya2 * $za1]
  set a1_cross_a2_y [expr $xa1 * $za2 - $xa2 * $za1]
  set a1_cross_a2_z [expr $xa1 * $ya2 - $xa2 * $ya1]

  set a3_cross_a4_x [expr $ya3 * $za4 - $ya4 * $za3]
  set a3_cross_a4_y [expr $xa3 * $za4 - $xa4 * $za3]
  set a3_cross_a4_z [expr $xa3 * $ya4 - $xa4 * $ya3]

  set E_dot_F [expr $a1_cross_a2_x * $a3_cross_a4_x + $a1_cross_a2_y * $a3_cross_a4_y + $a1_cross_a2_z * $a3_cross_a4_z]  
  set cross_angle [expr 180 * acos($E_dot_F) / 3.14159]
# puts $cross_angle
  return $cross_angle
}


