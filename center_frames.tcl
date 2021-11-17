# This script is used to move center of mass to origin in each frame 
# of multiple dcd files. 

	mol load psf gA_062602.psf dcd all_.dcd
	trace variable vmd_trajectory_read w move_box_to_center

# This script is used to move the whole box to center, to remove the 
# translation of center of mass during simulation
# Sep 2, 2002 by Zhanwu

proc move_box_to_center { } { 
	set num [molinfo top get numframes]
# TIP3 water molecules were excluded because wrapping.
# 
	set measure [atomselect top "not resname TIP3"]
	set all [atomselect top all]

	for { set i 0 } { $i < $num } { incr i } {
		$measure frame $i
		$all frame $i 
		set move_mat [measure center $measure weight mass]
		$all moveby [vecscale -1.0 $move_mat]
	}
animate write dcd centered_control_3.dcd
}

