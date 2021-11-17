proc distance {resid1 resid2 filename} {
set atom1 [atomselect top "protein and resid $resid1 and name CA"]
set atom2 [atomselect top "protein and resid $resid2 and name CA"]

set num [molinfo top get numframes]
	for {set frame 0} {$frame < $num} {incr frame} {
	$atom1 $frame
	$atom2 $frame
	set dist [vecdist $atom1 $atom2]
puts "$dist"
set fileId [open $filename a]
puts $fileId "$frame [format "%2.3f" $dist] 
close $fileId
}
}