#This script fits all frames in a simulation trajectory to the first frame
#according to the {sel_text},
#and calculates an average structure of all frames, 
#
#Replace the sel_text to what you want to select

proc average {filename} {
set sel_text "protein and (name N CA C)"

#assuming the trajectory of interest is molecule top

set selall [atomselect top protein]
set num_frame [molinfo top get numframes]

#align all frames with respect to 1st frame to remove translation and rotation
#puts "Entering alignment..."
set sel1 [atomselect top $sel_text]
set ref1 [atomselect top $sel_text frame 0]
        for {set frame 0} {$frame < $num_frame} {incr frame} {
                $sel1 frame $frame
                $selall frame $frame
                set trans_mat [measure fit $sel1 $ref1]
                $selall move $trans_mat
        }
$sel1 delete
$ref1 delete
puts "fitting complete"

#Calculating average structure
puts "entering averaging..."
    set factor [expr 1./$num_frame]

	for {set i 0} {$i < $num_frame} {incr i} {
	$selall frame $i
	set all_coor($i) [$selall get {x y z}]
	}
	$selall delete
       puts "all coordinates input"

    	set ave_coor {}
    	set len [llength $all_coor(0)]

    		for {set j 0} {$j < $len} {incr j} {
			set b [veczero]
			for {set i 0} {$i < $num_frame} {incr i} {
	    		set b [vecadd $b [lindex $all_coor($i) $j]]
			}
			lappend ave_coor [vecscale $b $factor]
			#puts "The $j th atom done"
    		}
puts "average calculated"

#If one wants to write a pdb file for the averaged atructure, uncomment the 
#next 3 lines. 

set sel2 [atomselect top protein frame 0]
$sel2 set {x y z} $ave_coor
$sel2 writepdb $filename.pdb
puts "finished"
}


