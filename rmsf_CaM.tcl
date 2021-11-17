#calculates an average structure of all frames, and then calculate the 
#root-mean-square fluctuation for atoms defined by {fluc_text}

proc rmsf { startframe endframe  filename } {
	
#------------------------------------------------------------------------
set file $filename 
#-------------------------------------------------------------------------

set fluc_text "name CA"

set selall [atomselect top all]
set num_frame [molinfo top get numframes]

if { $endframe > $num_frame }  {
	set endframe $num_frame
} 
 
set frame_count [expr $endframe - $startframe ]
puts "TOTAL $frame_count frames "
#align all frames with origin to remove translation
# Use center_of_mass to translate the frames
puts "translate according to center of mass"

for { set frame $startframe } { $frame < $endframe } { incr frame } {
	$selall frame $frame 
	set trans_mat [vecscale -1.0 [measure center $selall weight mass]]
	$selall moveby $trans_mat
}

#Calculating average structure
    set fluc [atomselect top $fluc_text]
    set factor [expr 1.0/$frame_count]

	for {set i $startframe } {$i < $endframe } {incr i} {
	$fluc frame $i
	set all_coor($i) [$fluc get {x y z}]
	}

 
    	set len [llength $all_coor($startframe)]

    		for {set j 0} {$j<$len} {incr j} {
			set b [veczero]
			for {set i $startframe } {$i < $endframe } {incr i} {
	    		set b [vecadd $b [lindex $all_coor($i) $j]]
			}
			set ave_vec [vecscale $b $factor]
			puts "$j  $ave_vec "
			set sum 0.0
			for {set i $startframe } {$i < $endframe } {incr i} {
	    		set sum [expr ($sum + [veclength2 [vecsub [lindex $all_coor($i) $j] $ave_vec]])]
			}
			set sum [expr (sqrt ($sum*$factor))]
			set fileId [open $file a]
			puts $fileId "[expr ($j+1)] [format "%.3f" $sum]"
			close $fileId
puts "residue $j $sum"
    		}

}
