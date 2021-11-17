set water [atomselect top "{same residue as within 1 of protein } and name OH2"]
set file [open waterlist.txt a]
puts here
foreach lip [$water get {segname resid} ] {
	puts $file "delatom $lip"
}
close $file
set water [atomselect top "name OH2 and {x > 42 or x < -42 or y > 42 or y < -42}"]

set file [open waterlist.txt a]
foreach lip [$water get {segname resid} ] {
	puts $file "delatom $lip"
}
close $file

