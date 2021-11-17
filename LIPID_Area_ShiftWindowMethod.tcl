#LIPID area measurement by "shifting window" method
#Is this kind of "bootstrap"?

#LOGIC:
#When there is protein in the membrane, it is difficult to measure the area accurately.
#Usual method includes measure the total and deduct the area occupied by protein, however,
# area of protein is also difficult to define. 
# The shifting window method use a window (let's say 15x15 squaure) of fixed size, shift
# in the membrane plane, and count how many P atoms in the window. This method should be 
# accurate, provided enough sampling.

# By Zhanwu Liu 
# Version 052704

proc LipidArea { } {
	# Define the scanning region
	set xmin -50
	set xmax +50
	set ymin -50
	set ymax +50
	
	set totalframe [molinfo top get numframes]
	
	#Define the stepsize of each movement of window, note that the calculation is N^2, 
	#so do not choose too small value. (small value of course gives better sampling)
	set stepsize 2  
	
	# the dimension of the box (side of square)
	set boxdimension 10 
	set halfdimension [expr $boxdimension /2.0 ]
	
	# The area definately occupied by protein, use radius value, e.g. within XX A of 
	# origin is occupied by protein, this should save a lot of time 
	
	set proteinoccupied 30

	# When the square center is in this range, definately overlap with protein
	set excludecenter [expr $proteinoccupied + $halfdimension]
	
	#maximam step number in the counting
	set maxstep [expr [expr [expr $xmax - $xmin] - $boxdimension]/$stepsize ]
	
	# The center of the starting box 
	set startcenterx [expr $xmin + $halfdimension]
	set startcentery [expr $ymin + $halfdimension]
	

	puts "exlude center $excludecenter; startcenterx $startcenterx"
	#Counting of valid windows and P atoms
	set windowcount 0.0
	set pcount 0.0

	# Loop through the area
for {set frame 1 } {$frame < $totalframe} { incr frame  1 } {
	set fileID [open "temp_lipidareafilewindow10.txt" a]
	for { set i 0 } { $i < $maxstep} { incr i } {
		for {set j 0 } { $j < $maxstep } { incr j} {
			
			set centerx [expr $startcenterx + $i * $stepsize ]	
			set centery [expr $startcentery + $j * $stepsize ]	
			
			# if area in the exclude region, get out of this loop
			if {[expr $centerx * $centerx + $centery * $centery] < $excludecenter*$excludecenter} { continue }

			set windowminx [expr $centerx - $halfdimension]
			set windowmaxx [expr $centerx + $halfdimension]
			set windowminy [expr $centery - $halfdimension]
			set windowmaxy [expr $centery + $halfdimension]

			set proteinatom [atomselect top "protein and x > $windowminx and x < $windowmaxx and y > $windowminy and y < $windowmaxy" frame $frame]
			# If there is protein atom in this area, get of of this loop

			if {[$proteinatom num]} { continue }

			set patom [atomselect top "name P1 and x > $windowminx and x < $windowmaxx and y > $windowminy and y < $windowmaxy" frame $frame]
			set windowcount [expr $windowcount + 1]
			set pcount [expr $pcount + [$patom num]]
			$patom delete 
			$proteinatom delete
		}
	}
			puts $fileID " Pcount = $pcount; Windows = $windowcount; area per lipid [expr ($windowcount*2*$boxdimension*$boxdimension)/$pcount ]"
			puts  " Pcount = $pcount; Windows = $windowcount; area per lipid [expr ($windowcount*2*$boxdimension*$boxdimension)/$pcount ]"
			close $fileID
	# divide 2 because of bilayer 
}
	set area [expr ($windowcount*2*$boxdimension*$boxdimension)/$pcount ]
	
	puts "$pcount P atoms, $windowcount Windows counted. The area per lipid is $area angstrom^2"

}
