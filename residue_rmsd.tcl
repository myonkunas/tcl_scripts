proc rmsd_residue_over_time {{mol top} res} {
    
    # use frame 0 for the reference
    set reference [atomselect $mol "protein and name C CA N" frame 0]
    # the frame being compared
    set compare [atomselect $mol "protein and name C CA N"]
    #make a selection with all atoms
    set all [atomselect top all]
    #get the number of frames
    set num_steps [molinfo $mol get numframes]
    #open file for writing
    set fil [open residue_rmsd.dat w]
    
    foreach r $res {
	set rmsd($r) 0
    }
    
    #loop over all frames in the trajectory
    for {set frame 0} {$frame < $num_steps} {incr frame} {
	puts "Calculating rmsd for frame $frame ..."
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$all move $trans_mat
	# compute the RMSD
	
	#loop through all residues
	foreach r $res {
	    set ref [atomselect $mol "name CA and protein and resid $r" frame 0]
	    set comp [atomselect $mol "name CA and protein and resid $r" frame $frame]
	    set rmsd($r) [expr $rmsd($r) + [measure rmsd $comp $ref]]
	    $comp delete
	    $ref delete
	   
	}
    }
    set ave 0
	foreach r $res {
	    set rmsd($r) [expr $rmsd($r)/$num_steps]
	    # print the RMSD
	    puts "RMSD of residue $r is $rmsd($r)"
	    puts $fil " $r \t $rmsd($r)"
	    set res_b [atomselect $mol "resid $r"] 
            $res_b set beta $rmsd($r)
            $res_b delete
	    set ave [expr $ave + $rmsd($r)]
	    
	}
    set ave [expr $ave/[llength $res]]
    puts " Average rmsd per residue:   $ave"
    close $fil
}

