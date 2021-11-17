proc get_halothane_interaction {halothane filename} {

set selection [atomselect top "name CA and same residue as protein and within 3 of segname $halothane"]
set num_frame [molinfo top get numframes]

for {set i 0} {$i < $num_frame} {incr i} {
$selection frame $i
set res [$selection get resname]
set id  [$selection get resid]
set fileId [open $filename a]
puts $fileId "$res $id"

close $fileId
	}
}
