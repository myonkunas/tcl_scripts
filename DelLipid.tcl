set lipid [atomselect top "{same residue as within 1 of protein } and name P"]
set file [open lipidlist.txt a]
puts here
foreach lip [$lipid get {segname resid} ] {
	puts $file "delatom $lip"
}
close $file

set lipid [atomselect top "name P and {x > 42 or x < -42 or y > 42 or y < -42}"]

set file [open lipidlist.txt a]
foreach lip [$lipid get {segname resid} ] {
	puts $file "delatom $lip"
}
close $file

