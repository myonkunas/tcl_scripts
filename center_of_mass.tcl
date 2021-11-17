# This script is used to plot the movement of center of mass during 
# the simulation. 

proc com_protein { filename } { 
	set fileId [open $filename a 0644]
	set numframe [molinfo top get numframes]
	set pro [atomselect top "segname SEG1 SEG2 "]
	
	for { set i 0 } { $i < $numframe } { incr i } {
		$pro frame $i
		set center [measure center $pro weight mass]
		puts " $i  $center "
		puts $fileId "$i $center "
	}
	close $fileId
}


proc com_lipid { filename } {  
        set fileId [open $filename a 0644]
        set numframe [molinfo top get numframes]
        set lip [atomselect top "segname LIP"]

        for { set i 0 } { $i < $numframe } { incr i } {
                $lip frame $i
                set center [measure center $lip weight mass]
                puts " $i  $center "
                puts $fileId "$i $center"
        }
	close $fileId
}

