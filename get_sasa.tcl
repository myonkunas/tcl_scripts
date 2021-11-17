proc SASA {residue filename} {
set sasa_res [atomselect top "protein and sidechain and resid $residue"]
set numframe [molinfo top get numframes]
for { set i 0 } { $i < $numframe } { incr i } {
	$sasa_res frame $i 
	set sa [measure sasa 1.4 $sasa_res]
set fileId [open $filename a]
puts $fileId "$i $sa"

close $fileId

	}
}