proc count {filename} {
set N [molinfo top get numframes]
set fileid [open $filename w+ 0777]

for { set i 0 } {$i < $N} {incr i 50} {

set all [atomselect top all]
set center [measure center $all weight mass]
set new_vector [vecscale -1.0 $center]
$all moveby $newvector
set sel [atomselect top "name OH2 and (same residue as x<5 and x>-5 and y<5 and y>-5 and z<7 and z>-5)" frame $i]
set count [$sel num]
set ID [$sel get resid]
puts -nonewline $fileid $i
puts -nonewline $fileid " "
puts -nonewline $fileid $count
puts -nonewline $fileid " "
puts $fileid $ID
}
close $fileid
}