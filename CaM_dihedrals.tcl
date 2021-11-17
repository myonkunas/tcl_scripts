#dihedrals.tcl
#Calculates the dihedral angle of any four atoms
#Similiar to example in VMD-L Mailing list

#Input the index of the four atoms and there mol ids and the outputfile name
proc dihedrals {index1 index2 index3 index4 outputfilename}  {

set mol 0

#label the dihedral angles
label add Dihedrals $mol/$index1 $mol/$index2 $mol/$index3 $mol/$index4
set dihedral [label graph Dihedrals 0]

#open the outputfile for writing
set fileid [open $outputfilename w+ 0777]

set n [molinfo top get numframes]

for {set i 0} {$i<$n} {incr i} {

#put in the dihedral angle
puts $fileid [lindex $dihedral $i]
}
#clear label for next run
label delete Dihedrals

#close outputfile
close $fileid

}
