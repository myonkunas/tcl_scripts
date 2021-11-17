proc get_pdbs {filename} {
 
set sel [atomselect top all]

	set num [molinfo top get numframes]
	for {set frame 0} {$frame < $num} {incr frame} {
		$sel frame $frame
		$sel writepdb $filename.$frame.pdb
        
}
}

