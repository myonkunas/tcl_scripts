#This script gets the theta euler angle for a vector along a bond
#Yonkunas May 31, 2005
proc theta {filename} {
set num [molinfo top get numframes]
for {set i 0} {$i<$num} {incr i} {

#select the atoms in the vector
	set atom1 [atomselect top "name C1 and segname HAL1"]
	set atom2 [atomselect top "name C2 and segname HAL1"]
#go to the correct frame
$atom1 frame $i
$atom2 frame $i

#get the coordinates so we can do operations on them
	set vec1 [lindex [$atom1 get {x y z}] 0]
	set vec2 [lindex [$atom2 get {x y z}] 0]

#get the normalized length of the vector
	set vec [vecsub $vec1 $vec2]
	set vecn [vecnorm $vec]

#project to the z axis to define the angle
	set out1 [vecdot $vecn { 0 0 1 }]
puts $out1

	set fileId [open $filename a]
	puts $fileId "$i [format "%2.3f" $out1]"
	close $fileId 

} 
}

#----------------------------------------------------------------------------------

#get the phi euler angle, the projection of the vector in the x-y plane
proc phi {filename} {
set num [molinfo top get numframes]
for {set i 0} {$i<$num} {incr i} {

#select the atoms in the vector
	set atom3 [atomselect top "name C1 and segname HAL1"]
	set atom4 [atomselect top "name C2 and segname HAL1"]
#go to the correct frame
$atom3 frame $i
$atom4 frame $i

#get the coordinates so we can do operations on them
	set vec3 [lindex [$atom3 get {x y }] 0]
	set vec4 [lindex [$atom4 get {x y }] 0]

#get the normalized length of the vector
	set vec2 [vecsub $vec3 $vec4]
	set vecn2 [vecnorm $vec2]

#project to the x axis to define the angle
	set out2 [vecdot $vecn2 { 1 0 }]
puts $out2

	set fileId [open $filename a]
	puts $fileId "$i [format "%2.3f" $out2]"
	close $fileId 
} 
}