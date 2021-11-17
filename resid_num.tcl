#Yonkunas
#renumber residues
proc numbering {chain start end shift} {
set sel [atomselect top "chain $chain and resid $start to $end"]
set residlist {}

for {set i $start} {$i<=$end} {incr i} {
 set number [[atomselect top "chain $chain and resid $i"] num]
 set residue [atomselect top "chain $chain and resid $i"]
 for {set j 1} {$j<=$number} {incr j} { 
                set resid [expr {$shift + $i}]
                lappend residlist $resid
        }
}
$sel set resid $residlist 

}