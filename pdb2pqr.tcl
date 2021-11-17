# By zhanwu Liu, Jun 18, 2003
# edited by Yonkunas May 2005 
# to convert psf/pdb file format into PQR format, 
# so it can be used in APBS to calculate electrostatic potential.

# Assume molecule already loaded in top. contains both psf and pdb
# e.g. contains charge information

proc pqr {filename} {
set all [atomselect top all]
set num [$all num]
set frame [molinfo top get numframes]
for {set i 0 } { $i < $frame } {incr i } {
	for {set j 0 } { $j < $num } { incr j } {
		$all frame $i
		set sel [atomselect top "index $j"]
		set charge1 [$sel get charge]
		set radius1 [$sel get radius]
		$sel set occupancy $charge1
		$sel set beta $radius1
puts "index $j"
}
$all writepdb $filename.$i.pqr
}
}