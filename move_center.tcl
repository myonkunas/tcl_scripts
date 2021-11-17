proc center {molid} {
set selall [atomselect $molid all]
set numframe [molinfo top get numframes]
for { set i 0 } { $i < $numframe } { incr i } {
	$selall frame $i 
	set trans_mat [vecscale -1.0 [measure center $selall weight mass]]
	$selall moveby $trans_mat
}
}