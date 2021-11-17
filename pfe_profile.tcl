# 07-26-2002 Zhanwu Liu, To calculate the profile of PFE in membrane
# 08-13-2002 Zhanwu Liu, include two functions:
# pfe_profile is used to calculate profile in z-axis
# pfe_profile_xy is used to calculate profile in x-y plane
#VMD  --- start of VMD description block
#Name:
# Trajectory path
#Synopsis:
# Draws the path of the center of mass of a selection through an animation
#Version:
# 1.0
#Uses VMD version:
# 1.1
#Ease of use:
# 2
#Procedures:
# <li>trajectory_path selection {color blue} -- follows the center of
# mass of the given selection.  Color is a solid color, or "scale" for color scale.
#Description:
# For each step in the animation, the center of mass of the selection is
# calculated.  A new graphics molecule is created containing lines connected
# successive coordinates.  The color is a solid color (default is blue) or
# they are mapped to the color scale from lowest (= first trajectory frame)
# to highest (= last trajectory frame).
#Example:
# <pre>
# set water [atomselect top "retopsid 5243"]
# trajectory_path $water scale
#Files: 
# <a href="trajectory_path.vmd">trajectory_path</a>
#Author: 
# Andrew Dalke &lt;dalke@ks.uiuc.edu&gt;
#\VMD  --- end of block

# each PFE molecule has different resid 
#set sel [atomselect 0 "segname ANE"]
#trajectory_path $sel $i

proc pfe_profile_z { filename } {
 # save the current selection frame number
 # make the list of coordinates
 set num_frames [molinfo top get numframes]
 set fileId [open $filename a 0644]
#initialize the list for coordinate
 set framenum {}
 for { set i 0 } { $i < 10 } { incr i 1 } { 
	set coord$i {}
	}

 for {set i 0} {$i < $num_frames} {incr i 1 } {
	lappend framenum $i
	for { set j 0  } { $j < 10 } { incr j} {
		set sel [atomselect top  "segname PFE and resid $j"]
   		$sel frame $i
   # compute the center of mass and save it on the list
   	set center [measure center $sel weight mass]
# 2 is z -axis, if y, chose 1
	lappend coord$j [lindex $center 2]

 }
puts "frame $i finished"
}
puts "Data acquiring finished"
 ##### now make the graphics
 # make a new molecule for this selection

foreach framenumb $framenum p0 $coord0 p1 $coord1 p2 $coord2 p3 $coord3 p4 $coord4 p5 $coord5 p6 $coord6 p7 $coord7 p8 $coord8 p9 $coord9 {
puts $fileId "$framenumb $p0 $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 "
	} 
close $fileId
 # use the color scale?
}



# The function is used to plot the pfe_profile on x-y plane
proc pfe_profile_xy { filename } {
 # save the current selection frame number
 # make the list of coordinates
 set num_frames [molinfo top get numframes]
 set fileId [open $filename w]
#initialize the list for coordinate
 #set framenum {}
 #for { set i 0 } { $i < 10 } { incr i 1 } { 
#	set coords$i {}
#	}

 for {set i 0} {$i < $num_frames} {incr i 1 } {
	lappend framenum $i
	for { set j 0  } { $j < 10 } { incr j} {
		set sel [atomselect top  "segname PFE and resid $j"]
   		$sel frame $i
   # compute the center of mass and save it on the list
   	set center [measure center $sel weight mass]
# 2 is z -axis, if y, chose 1, and x chose 0
	set x [lindex $center 0]
	set y [lindex $center 1]
	set z [expr sqrt([expr $x*$x + $y*$y])]
	lappend coord$j $z
puts " Frame $i, resdue $j, coord coord$j"

 }
}
puts "Data acquiring finished"
 ##### now make the graphics
 # make a new molecule for this selection

foreach framenumb $framenum p0 $coord0 p1 $coord1 p2 $coord2 p3 $coord3 p4 $coord4 p5 $coord5 p6 $coord6 p7 $coord7 p8 $coord8 p9 $coord9 {
puts $fileId "$framenumb $p0 $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 "
	flush $fileId
	} 

 # use the color scale?
}



# This function is used to get the velocity profile of C2F6 movement. (distance/unit time, maybe 5 ps).
proc pfe_movement {filename } {

	set num_frames [molinfo top get numframes]
 	set fileId [open $filename w]

# This loop to get all the coordinates out
 	for { set i 0 } { $i < $num_frames } { incr i 1 } { 
		lappend framenum $i 
# 10 is the number of C2F6 molecules
		for { set j 0 } { $j < 10 } { incr j } { 
			set sel [atomselect top "segname PFE and resid $j"]
			$sel frame $i
			set center [measure center $sel weight mass]
			lappend coord($j) $center
		}
	}
	puts "coodinate extraction finished!"

# This loop to calculate the values 
	set num_points [expr $num_frames / 1]
	for { set i 5 } { $i < $num_points } { incr i 1 } { 
# 10 is the number of C2F6 molecules
		for { set j 0 } { $j < 10 } { incr j } { 
			#puts "$i [lindex $coord($j) $i]  [lindex $coord($j) [expr $i-5]]"
			set distance [vecdist [lindex $coord($j) $i]  [lindex $coord($j) [expr $i-5]]]	
 			lappend dis$j $distance
		}
	}

# To print out all the values 
	foreach framenumb $framenum p0 $dis0 p1 $dis1 p2 $dis2 p3 $dis3 p4 $dis4 p5 $dis5 p6 $dis6 p7 $dis7 p8 $dis8 p9 $dis9 {
		puts $fileId "$framenumb $p0 $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 "
		flush $fileId
	} 

}
