#This script fits all frames in a simulation trajectory to the 1st frame
#according to the {fit_text},
#calculates an average structure of all frames, and then calculate the 
#root-mean-square fluctuation for atoms defined by {fluc_text}
#
proc rmsf { startframe endframe seg filename } {
	
#------------------------------------------------------------------------
set file $filename 
#-------------------------------------------------------------------------

#Replace the fit_text and fluc_text to what you want to select

set fit_text "segname $seg and name CA"
set fluc_text "segname $seg and name CA"

#assuming the trajectory of interest is molecule 0

set selall [atomselect top all]
set num_frame [molinfo top get numframes]

if { $endframe > $num_frame }  {
	set endframe $num_frame
} 
 

set frame_count [expr $endframe - $startframe ]
puts "TOTAL $frame_count frames "
#align all frames with respect to the 1st frame to remove translation and rotation

set sel1 [atomselect top $fit_text]
set ref1 [atomselect top $fit_text frame 0]

# Use center_of_mass to translate the frames
puts "translate according to center of mass"
set protein [atomselect top "segname $seg"]

for { set frame $startframe } { $frame < $endframe } { incr frame } {
	$protein frame $frame
	$selall frame $frame 
	set trans_mat [vecscale -1.0 [measure center $protein weight mass]]
	$selall moveby $trans_mat
}
#Calculating average structure
    set fluc [atomselect top $fluc_text]
    set factor [expr 1./$frame_count]

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
