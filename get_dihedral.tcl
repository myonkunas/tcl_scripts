#Yonkunas 10/26/04
#Get dihedrals of residues starting with start and ending with end
#for every frame in simulation

proc get_dihedral { start end filename } { 

set num [molinfo top get numframes]
set fileId [open $filename a]

for { set j $start } { $j <= $end } { incr j } { 
	for { set i 0 } { $i < $num } { incr i + 20 } {
		set atom [atomselect top "resid $j and name CA" frame $i]
		$atom frame $i
		puts $fileId "$j  [$atom get phi] [$atom get psi]"
	}
	}
	close $fileId
}
