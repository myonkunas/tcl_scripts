proc move {mol} {
set selall [atomselect $mol all]
set num_frame [molinfo top get numframes]
for { set frame 0 } { $frame < $num_frame } { incr frame } {
	$selall frame $frame 
	set trans_mat [vecscale -1.0 [measure center $selall weight mass]]
	$selall moveby $trans_mat
}
}