proc distance {trajname} {

mol new 1k4c_patch_TM1.pdb waitfor all
mol addfile $trajname.trj type crdbox first 1 last -1 waitfor all

set atom1 [atomselect top "protein and resid 45 and name CA"]
set atom2 [atomselect top "protein and resid 94 and name CA"]
set atom3 [atomselect top "protein and resid 143 and name CA"]
set atom4 [atomselect top "protein and resid 192 and name CA"]

set num [molinfo top get numframes]
	for {set frame 0} {$frame < $num} {incr frame} {
	$atom1 frame $frame
  	set a1 [lindex [$atom1 get { x y z }] 0] 
	$atom2 frame $frame
  	set a2 [lindex [$atom2 get { x y z }] 0] 
	$atom3 frame $frame
  	set a3 [lindex [$atom3 get { x y z }] 0] 
	$atom4 frame $frame
  	set a4 [lindex [$atom4 get { x y z }] 0] 

	set dist1 [vecdist $a1 $a2]
	set dist2 [vecdist $a3 $a4]

puts "$dist1 $dist2"
set fileId [open $trajname.dat a]
puts $fileId "$frame [format "%2.3f" $dist1] [format "%2.3f" $dist2]" 
close $fileId
}
}

