# This script is used to give a plot of H-bond
# Zhanwu Liu, Oct 02, 2002 
## set distance between donor-acceptor to 3.5, the corresponding distance between hydrogen
## and acceptor is about 2.5 (N-H bond length is 0.98 A)

set distance 3.5
## set angle to 20 degree, the <donor-acceptor-Hydrogen angle, the smaller, the stringent
set angle 30

set donor [atomselect top "name NE1 and segname SEG1 and resid 15"]
#set donor [atomselect top "name HE1 and segname SEG1"]
set acceptor2 [atomselect top {name "O.*" and resname TIP3}] 
set acceptor1 [atomselect top {name "O.*" and segname LIP }]

##set count every $step frames
set step 5
set filename "wat_PFE_SEG1_15.dat"



set fileId [open $filename a 0644]

set framenum [molinfo top get numframes]

for { set i 0 } { $i < $framenum } { incr i $step }  {
	$donor frame $i
	$acceptor1 frame $i
	$acceptor2 frame $i
	set lipidnum [measure hbonds $distance $angle $donor $acceptor1]
	set waternum [measure hbonds $distance $angle $donor $acceptor2]
	
	set lipidbondnum [llength [lindex $lipidnum 0]]
	set hindex [lindex $lipidnum 2]
	set waterbondnum [llength [lindex $waternum 0]]
        set total  [expr $lipidbondnum + $waterbondnum]
	#puts $fileId "$i $lipidbondnum"
	puts $fileId "$i $waterbondnum"
	#puts $fileId "$i  $lipidbondnum  $waterbondnum $total "
	
}
	
close $fileId
			
