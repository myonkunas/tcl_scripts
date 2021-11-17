#Michael Yonkunas Spring 2004
#Water counting script written for KSI
#Count the number of waters at the interface according to the size limits
#change the increment to adjust number of frames between data points (100 over 1000frames = 10 data points)
#stride data accordingly (substrate data strided by 10 = 5ps/frame)

proc count {increment filename} {

set N [molinfo top get numframes]
set fileid [open $filename w+ 0777]

for { set i 0 } {$i < $N} {incr i $increment} {

set sel [atomselect top "name OH2 and (same residue as x<5 and x>-5 and y<5 and y>-5 and z<7 and z>-5)" frame $i]
set count [$sel num]
puts -nonewline $fileid $i
puts -nonewline $fileid " "
puts $fileid $count
}
close $fileid
}