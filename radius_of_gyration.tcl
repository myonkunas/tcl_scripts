proc gyr_radius {filename} {
set sel [atomselect top "protein"]
set num [molinfo top get numframes]
	for {set frame 0} {$frame < $num} {incr frame} {
		$sel frame $frame
# gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  		set com [measure center $sel weight mass]
  		set sum 0
  		foreach coord [$sel get {x y z}] {
    		set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
#  		return [expr sqrt($sum / ([$sel num] + 0.0))]
set ROG [expr sqrt($sum / ([$sel num] + 0.0))]

set fileId [open $filename w]
puts $fileId " $frame \t $ROG"
close $fileId
}
}
}