proc get_residues {residue filename} {

set selection [atomselect top "name C and resname $residue"]

set b [llength [$selection get segname]]
for { set a 0 } {$a < $b} {incr a} {

set fileId [open $filename a]
puts $fileId "[lindex [$selection get resname] $a] [lindex [$selection get resid] $a] [lindex [$selection get segname] $a]"
close $fileId
	}
}