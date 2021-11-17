# Lipid bilayer analysis script 
# Zhanwu Liu, Department of Anesthesiology, University of Pittsburgh
# Apr 03, 2002
# This script will contain these functions:
# 1. atom density profile 
# 2. electrostatic potential
# 3. radial distribution function 
# 4. head group orientation order parameter
# 5. area per lipid
# 6. Volume per lipid
# 7. Scd (variation of order parameters )
# 8. interfacial width ( 90% water -> 10 % water )
# 9. orientation time correlation function
# 10. Diffusion constant
# 11. Memebrane thickness analysis 
#	thickness firstframe lastframe filename 

#################################################################################################
######################## Membrane thickness profile around Protein ##############################
#/* Jun 13, 2002 by Zhanwu Liu */
# To calculate the z-plot according to distance with protein, useful in analyzing the interaction
# between protein and lipid bilayer.

# This procedure is used to calculate the value, three inputs: first frame, last frame, outfilename.
proc Dpp_profile { firstframe lastframe filename } {
	set fileID [open $filename a 0644]
	# Maxim distance between Lipid and protein
	set max_dis 60
	# The origin of protein, in our calculation, y is membrane normal 
	set x0 0 
	set z0 0
	# Step size of distance, in A unit
	set step 1
	# Count and count2 are count for $normal> 0 and $normal< 0
	for { set dis $step } { $dis < $max_dis } { incr dis $step } {
		set count($dis) 0
		set count2($dis) 0
		set data($dis) 0 
		set data2($dis) 0
		}

	set totalframe [molinfo top get numframes]
	if { $lastframe > $totalframe } { 
		set lastframe $totalframe
		}
	for { set i $firstframe } { $i <= $lastframe } { incr i 1 } { 
		set last_dis 0
		puts "Processing Frame $i"
		for { set dis  $step } { $dis < $max_dis } { incr dis $step} {
			set atoms [atomselect top "name P and sqrt(x^2+z^2)< $dis and sqrt(x^2 + y^2)>= $last_dis"] 
			set last_dis $dis
			foreach normal [$atoms get y] {
				if { $normal > 0 } { 
					incr count($dis) 1
					set data($dis) [expr $data($dis) + $normal]
					} else { 
						incr count2($dis) 1					
						set data2($dis) [expr $data2($dis) + $normal]
						}
				}
			}
		}
	
	for { set dis $step } { $dis < $max_dis } { incr dis $step } { 
		if { $count($dis) == 0 } {
			set count($dis) 1
			}
		if { $count2($dis) == 0 } {
			set count2($dis) 1 
			}
		puts $fileID "$dis  [expr $data($dis) / $count($dis)]  [expr $data2($dis) / $count2($dis)]"
	 }

close $fileID  
}


######################## End of Membrane thickness profile around Protein #######################
#################################################################################################

##################################################################################################
########################## ATOM DENSITY PROFILE ##################################################
# /* Apr 03, 2002 by Zhanwu Liu */

# This procedure is used to calculated the atom density profile along lipid bilayer normal 
proc atomdensity { firstframe lastframe axis filename } {
	set fileID  [open $filename a 0644]
	set totalframe [molinfo top get numframes]
	if {  $lastframe > $totalframe } {
		set lastframe $totalframe
		}
# number of points in normal axis
	set pointnum 100
	set maxz  40.0
	set minz -40.0
# The interested part, such as choline, or H2O
	set selectword "name P"
# calculate some useful values
	set dz [expr ($maxz - $minz )/$pointnum ]
	
	puts "dz $dz lastframe $lastframe "
# cycle for each frame
	for { set j 0 } { $j < $pointnum } { incr j } { 
		set lower [expr $minz + $dz*$j]
		set higher [expr $lower + $dz ]
		puts "lower $lower higher $higher "
		set a 0.0
		for { set i $firstframe } { $i < $lastframe } { incr i } { 
			set atoms [atomselect top "$axis < $higher and $axis > $lower and $selectword" frame $i]	
			foreach coor [$atoms get $axis] atommass [$atoms get mass ] {
				set a [expr $a+ $atommass]
				}
		}		
		puts "$lower  $a"
		puts $fileID "$lower  $a"
	}
	close $fileID
}


######################### END OF ATOM DENSITY PROFILE ############################################
##################################################################################################

##################################################################################################
############################  THICKNESS CALCULATION ##############################################
# /* Mar 25, 2002 by Zhanwu Liu */

# This procedure is to calculate the thickness of lipid bilayer
# by the P atom
proc thickness { firstframe lastframe filename } { 
	set fileID [open $filename a 0644]
	set totalframe [molinfo top get numframes]
	if { $firstframe < 0 || $lastframe > $totalframe } { 
		puts "Wrong frame number"
		exit
		} else {
		for { set i $firstframe } { $i < $lastframe } { incr i 1 } {
			set strings [thickness_frame top z $i]
			puts  $fileID  $strings 
			}
		}
	close $fileID
}
	

# This procedure is used to calculate the thickness of lipid bilayer
# in one frame, normalaxis normal to the bilayer plane
# Zhanwu, 11-05-2002, modify to return these values:
# Thickness, average z of two layers, SD of two layers.
# SD = sqrt ( <z*z> - <z>*<z>)
proc thickness_frame { molid normalaxis framenum } {
set p1 [atomselect $molid "name P and $normalaxis > 0 "]
set p2  [atomselect $molid "name P and $normalaxis < 0 "]

$p1 frame $framenum
$p2 frame $framenum

set p1a 0
set p2a 0
set z1square   0
set z2square   0

foreach zp1 [ $p1 get $normalaxis ] {
	set p1a [expr $p1a +  $zp1 ]
	set z1square [expr $z1square+ [expr $zp1 * $zp1]]
}
foreach zp2 [ $p2 get $normalaxis ] {
        set p2a [expr $p2a + $zp2]
	set z2square [expr $z2square+ [expr $zp2 * $zp2]]
}   

set z1 [expr $p1a/[$p1 num]]
set z2 [expr $p2a/[$p2 num ]]

set thick    [expr $z1 - $z2]

set sd1 [expr sqrt([expr $z1square/[$p1 num] - $z1 * $z1])]
set sd2 [expr sqrt([expr $z2square/[$p2 num] - $z2 * $z2])]
#puts stdout "thickness of the membrane is $thick"
# Output the thickness as well as the average value in both side
set returnstring "$framenum   $thick  $z1   $z2 $sd1 $sd2"
return $returnstring
}

####################### END of THICKNESS CALCULATION ##################################################
#######################################################################################################

# This procedure is used to calculate the thickness of lipid bilayer
# in one frame, normalaxis normal to the bilayer plane
proc thickness_atom { atom } {
set normalaxis y
set p1 [atomselect top "name $atom  and $normalaxis > 0 "]
set p2  [atomselect top "name $atom  and $normalaxis < 0 "]


set p1a 0
set p2a 0

foreach zp1 [ $p1 get $normalaxis ] {
        set p1a [expr $p1a +  $zp1 ]
}
foreach zp2 [ $p2 get $normalaxis ] {
        set p2a [expr $p2a + $zp2]
}

set z1 [expr $p1a/[$p1 num]]
set z2 [expr $p2a/[$p2 num ]]

set thick    [expr $z1 - $z2]
#puts stdout "thickness of the membrane is $thick"
# Output the thickness as well as the average value in both side, and also the standard deviation
set returnstring "$thick  $z1   $z2 "
puts $returnstring
}

####################### END of THICKNESS CALCULATION ##################################################
#######################################################################################################


#######################################################################################################
####################### BEGIN OF P ATOM DISTRIBUTION CALCULATION ######################################
proc p_distribution {firstframe lastframe filename} {

set framenum [molinfo top get numframes]
	if { $lastframe > $framenum } { 
		puts " Last frame changed to $framenum " 
		set lastframe $framenum
	}


# select phosphate atom as an example, change if needed.
set Patom [atomselect top "name P"]

# The number of frames
set framecount [expr $lastframe - $firstframe]
# The number of atoms
set atomcount [$Patom num]  
# Total count divide factor equals 1.0
set factor [expr 1.0/[expr $framecount*$atomcount]]

#initialize the list used to store counting data
# stepsize 0.1 A 
for { set i -80 } { $i < 80 } { incr i 1 } { 
	set count($i) 0 
}


for { set i $firstframe } { $i < $lastframe } { incr i 1 } {
	$Patom frame $i
	foreach z [$Patom get z ] {
		set round [expr round([expr $z*2])]
		set count($round) [expr $count($round)+1]
	}		
}

set fileId [open $filename a 0644] 
for { set i -80 } { $i < 80 } { incr i 1 } { 
	set j [expr $i/2.0]
	set outdata [expr $count($i)*$factor]
	puts $fileId "$j $outdata"
#	puts "$i $count($i)"
}	

close $fileId
# End of program
}

####################################################################################
####################################################################################




