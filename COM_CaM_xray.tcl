proc center_of_mass { filename } { 

set fileId [open $filename a]
set numframe [molinfo top get numframes]

set domain1 [atomselect top "protein and resid 4 to 61"]
set domain2 [atomselect top "protein and resid 90 to 147"]

	for { set i 0 } { $i < $numframe } { incr i } {
$domain1 frame $i
$domain2 frame $i
set center1 [measure center $domain1 weight mass]
set center2 [measure center $domain2 weight mass]

			set 1_x [lindex $center1 0]
			set 1_y [lindex $center1 1]
			set 1_z [lindex $center1 2]

			set 2_x [lindex $center2 0]
			set 2_y [lindex $center2 1]
			set 2_z [lindex $center2 2]

set distance [expr {pow(pow($1_x-$2_x,2)+pow($1_y-$2_y,2)+pow($1_z-$2_z,2),0.5)}]
puts "$i $distance"
puts $fileId "$i $distance"
}
close $fileId
}

