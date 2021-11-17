#VMD  --- start of VMD description block
#Name:
# Trajectory path
#Synopsis:
# Draws the path of the center of mass of a selection through an animation
#Version:
# 1.1
#Uses VMD version:
# 1.8
#Ease of use:
# 2
#Procedures:
# <li>trajectory_path selection {color blue} -- follows the center of
# mass of the given selection.  Color is a solid color, or "scale" for color scale.
#Description:
# For each step in the animation, the center of mass of the selection is
# calculated.  A new graphics molecule is created containing lines connected
# successive coordinates.  The color is a solid color (default is blue) or
# they are mapped to the color scale from lowest (= first trajectory frame)
# to highest (= last trajectory frame).
#Example:
# <pre>
# set water [atomselect top "resid 5243"]
# trajectory_path $water scale
#Files: 
# <a href="trajectory_path.vmd">trajectory_path</a>
#Author: 
# Andrew Dalke &lt;dalke@ks.uiuc.edu&gt;
#\VMD  --- end of block

proc trajectory {selection {color blue} } {
 # save the current selection frame number
 set sel_frame [$selection frame]
 # make the list of coordinates
 set num_frames [molinfo [$selection molindex] get numframes]
 set coords {}
 #Open file to write data 
 #set fileid [open $filename a]
 
 for {set i 0} {$i < 1307} {incr i} {
   $selection frame $i
   # compute the center of mass and save it on the list
   lappend coords [measure center $selection weight mass]
   set coords_2 [measure center $selection weight mass]
#   puts -nonewline $fileid $i
#   puts -nonewline $fileid "	"
#   puts $fileid [format "%2.3s" $coords_2] 
#puts $fileid $coords_2
 }
# close $fileid
 ##### now make the graphics
 # make a new molecule for this selection
 #mol new graphics [$selection text]
 set gr_mol [molinfo top]
 set coords [lassign $coords prev]

 # use the color scale?
 if {$color == "scale"} {
   set count 0
   incr num_frames
   foreach coord $coords {
     #set color [expr 34 + int(32 * ($count + 0.0) / ($num_frames + 0.0))]
     set color [expr 17 + int(1024 * ($count + 0.0) / ($num_frames + 0.0))]
     graphics $gr_mol color $color
     graphics $gr_mol line $prev $coord
     set prev $coord
     incr count
   }
 } else {
   # constant color
   graphics $gr_mol color $color
   foreach coord $coords {
     graphics $gr_mol line $prev $coord
     set prev $coord
   }
 }

 # return the selection to its original state
 $selection frame $sel_frame
}

