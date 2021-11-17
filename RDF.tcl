#!/usr/bin/wish
# Written by Zhanwu Liu at Feb 12, 2002 for RDF calculation
# Department of Anesthesiology, University of Pittsburgh
# Known bugs: Cutoff, mbin x y z value should include digital point, 
# otherwise it will consider as integer operation, so cause trouble.
# Mar 21 2002, these bugs are removed, I just realized that tcl could cast 
# the type

set rdf  [toplevel .rdfcalc]
wm title $rdf "RDF calculation"

# List of variables used:
# Frame: startnum, endnum, skipnum, framenum
# Atoms: Segname_from, Resname_from, Name_from
# Atoms: Segname_to, Resname_to, Name_to
# Regional: cutoff, mbin
# Box: box_x, box_y, box_z 
# Control commands: Start Quit Cancel
# Data: set as global, rdfdata store data, sets is the frame index 
# rdfdata is used to store the accumulation number of atoms
# which rdffunc is used to store the g(r) value
set rdfdata(0) 0
set rdffunc(0) 0
set M_PI
set box_x 0
set box_y 0
set box_z 0

#number of trajectory frames
set framenum [molinfo top get numframes]

#construct the user interface
#create a frame for trajectory frame numbers
frame $rdf.top -borderwidth 10 -relief raised
pack $rdf.top -side top -fill x

# Create the region showing total frame number 
label $rdf.top.l -text "Total $framenum Frames "  -relief sunken
pack $rdf.top.l -side top 

foreach key { start end skip } {
frame $rdf.top.$key -borderwidth 10
label $rdf.top.$key.l -text "$key"
entry $rdf.top.$key.frame -width 8 -relief  sunken -background pink 
pack $rdf.top.$key -side left
pack $rdf.top.$key.l $rdf.top.$key.frame -side top

}

#frame for Atom information
frame $rdf.atominfo -borderwidth 10 
pack $rdf.atominfo -side top
foreach atom  { from to } {
	frame $rdf.atominfo.$atom -borderwidth 10
	label $rdf.atominfo.$atom.l  -text " $atom atom "
	pack $rdf.atominfo.$atom.l -side top -fill x
	foreach key { segname resname name } {
		frame $rdf.atominfo.$atom.$key 
		label $rdf.atominfo.$atom.$key.l -text "$key"
		entry $rdf.atominfo.$atom.$key.key -relief sunken -width 8 -background yellow
		pack $rdf.atominfo.$atom.$key.l $rdf.atominfo.$atom.$key.key -side left
		pack $rdf.atominfo.$atom.$key -side top
		}
	pack $rdf.atominfo.$atom -side left 
} 

#Frame for data infomation 
# All these data are accessed by frame_RDF function
frame $rdf.data -borderwidth 10
frame $rdf.data.cutoff
	label $rdf.data.cutoff.l -text "Cutoff"
	entry $rdf.data.cutoff.num -width 8 -background black -foreground white
pack $rdf.data.cutoff.l $rdf.data.cutoff.num -side top

frame $rdf.data.mbin 
	label $rdf.data.mbin.l -text "mbin"
	entry $rdf.data.mbin.num -width 8 -background black -foreground white
pack $rdf.data.mbin.l $rdf.data.mbin.num -side top
pack $rdf.data -side top
pack $rdf.data.cutoff $rdf.data.mbin -side left


# Frame for boundary information
# Not implemented, so just disabled the entry part
frame $rdf.boundary -borderwidth 10
pack $rdf.boundary -side top
label $rdf.boundary.l -text "Box Dimension"
#label $rdf.boundary.l2 -text " PBC NOT IMPLEMENTED " -foreground red
pack $rdf.boundary.l -side top -fill x
foreach dimension { x y z } {
	frame $rdf.boundary.$dimension
	label $rdf.boundary.$dimension.l -text "$dimension" 
	entry $rdf.boundary.$dimension.num  -width 8 -background gold 
	pack $rdf.boundary.$dimension.l $rdf.boundary.$dimension.num -side top
	pack $rdf.boundary.$dimension -side left
}
	

# Frame for command button
frame $rdf.cmdbutton 
button $rdf.cmdbutton.run -text "Run RDF" -command Run 
button $rdf.cmdbutton.quit -text "Quit" -command "quitRDF"
button $rdf.cmdbutton.save -text "Save RDF" -command "printrdf"
pack  $rdf.cmdbutton.run $rdf.cmdbutton.save -side left
pack $rdf.cmdbutton.quit -side right
pack $rdf.cmdbutton -side top
# Frame for file
frame $rdf.file
label $rdf.file.label -text "Input the file name"
entry $rdf.file.name -width 20
pack  $rdf.file.label $rdf.file.name -side top
pack $rdf.file -side top

#========================================================================
#===================== frame_RDF procedure ==============================

# procedure frame_RDF used to calculate RDF value in one frame 
# Do not process PBC now!
proc frame_RDF  { fromatomlist toatomlist cutoff dr } { 
# Frame for file
global rdfdata
global sets
# id key is used to recognize different molecules, maybe segname 
# if you use psfgen to replicate molecules. Or it might
# be resid
#set idkey "resid"
set idkey "segname"

#move to the frame 
set fromcoorlist [$fromatomlist get { x y z }] 
set fromidlist  [$fromatomlist get $idkey ]
set tocoorlist [ $toatomlist get { x y z } ]
set toidlist [$toatomlist get $idkey]

#Accumulate the atomnumber of in certain bins. m is used to differentiate
# atoms from same molecules so that it would not be counted in the calculation.
#if you want to turn off this feature, simply remove this line if { $m == $n }

foreach coordfrom $fromcoorlist m $fromidlist {
	foreach coordto $tocoorlist n $toidlist { 
		if { $m == $n } continue
		set difference [vecsub $coordfrom $coordto]

		set distance [veclength $difference]
		if { $distance < $cutoff } {
		       set group [expr $distance / $dr] 	
		       set groupint [expr int ($group) ]
			set rdfdata($groupint) [expr $rdfdata($groupint)+1]
		    } else continue
                } 
	} 
}	
#============================================================================


#============================================================================
#==================== Run procedure =========================================

# This is the main program 
# 
proc Run { } {
global rdf
global rdfdata
#sets is the number of frames calculated so far, is used to calculate
#average value 
set sets 0
global M_PI
set startnum [$rdf.top.start.frame get]
set endnum [$rdf.top.end.frame get]
set cutoff [expr double ([$rdf.data.cutoff.num get])]
set mbin [$rdf.data.mbin.num get ]
set skipnum [$rdf.top.skip.frame get]
    if { $skipnum <= 0 } {
	puts "Skip should larger than 0 "
        return
	}
global framenum
set dr [expr $cutoff /  $mbin ]
global rdffunc
global box_x
global box_y
global box_z
set box_x [expr double ([$rdf.boundary.x.num get])]
set box_y [expr double ([$rdf.boundary.y.num get])]
set box_z [expr double ([$rdf.boundary.z.num get])]



#build atom selection keywords for Fromatom
set fromstring " "
foreach key { segname resname name } {
set frominfo [$rdf.atominfo.from.$key.key get ]
if { $frominfo != "*" } {
	append fromstring "$key $frominfo "
  if { $key != "name" } { 
     append fromstring " and "
	}
    }
}


#build atom selection keywords for Toatom 
set tostring " "
foreach key { segname resname name } {
set toinfo [$rdf.atominfo.to.$key.key get ]
if { $toinfo != "*" } {
        append tostring "$key $toinfo "
  if { $key != "name" } { 
     append tostring " and "
        }
    }
}

#build selectionlist
set fromatomlist [atomselect top "$fromstring"]
set toatomlist [atomselect top "$tostring" ]
#calculate the density of to_atoms
set volume [expr $box_x * [expr $box_y * $box_z]]
set density [expr [expr [$toatomlist num]/$volume ]*[$fromatomlist num]]
puts "Density $density"

# If the number of frame is wrong?
# else initialize the store space for rdfdata, and calculate
if { [expr { $startnum > $framenum}] ||[expr {$endnum >=$framenum }] || [expr { $startnum > $endnum}] } {
puts "Wrong frame number "
} else { 
	for { set i 0 } { $i < $mbin } { incr i 1 } {
      set rdfdata($i) 0}

for { set i $startnum } { $i <= $endnum } { incr i $skipnum } {

# Call frame_RDF to calculate the numbers
# First get the information of selections of each frame
puts "Now processing Frame number $i, please be patient "
set t1 [clock seconds]
#rapidwrap $box_x $box_y $box_z $i
$fromatomlist frame $i
$toatomlist frame $i
frame_RDF $fromatomlist $toatomlist $cutoff $dr
set t2 [clock seconds]
puts "This frame cost [expr $t2 - $t1] seconds"
incr sets 1
}

# calculate the constant value for volume calculation
# deltaV = 4/3*M_PI*$dr*$dr*$dr +4*M_PI*$dr*r1*r2
set const1 [expr [expr 1.3333333* $M_PI] * [expr pow( $dr, 3)]] 
set const2 [expr [expr 4 * $M_PI] * $dr]
puts "dr $dr  const1 $const1  const2 $const2"

#processing of data, RDF calculation
puts "RDF function value"
for { set i 0 } { $i < $mbin } { incr i 1 } {
   set r1 [expr $dr*$i]
   set r2 [expr $r1+$dr]
   set deltaV [expr $const1+[expr $const2*[expr $r1*$r2]]]
   set rdffunc($i) [expr [expr $rdfdata($i)/$density]/[expr  $deltaV * $sets]]
   puts "[expr $i*$dr]     $rdffunc($i)"
}

# Draw on VMD display
draw color red
draw line {0 0 0 } {20 0 0 }
draw line {0 0 0 } { 0 20 0 }
draw line {0 0 0 } { 0 0 20 }
draw color white
for { set j 2 } { $j < $mbin } { incr j 1 } {
  set first "[expr [expr $j-1]*$dr]  $rdffunc([expr $j-1]) 0 "
  set second "[expr $j*$dr]  $rdffunc($j) 0"
  draw line $first $second
}
}
}
#=====================================================================================


#=====================================================================================
#======================== quitrdf procedure ==========================================
# To aviod quit by mistake! Prefer to cancel the quit command
proc quitRDF { } {
global rdf
frame $rdf.quitrdf
label $rdf.quitrdf.quitmessage -text "Really want to quit?" -foreground red
button $rdf.quitrdf.confirmquit -text "Quit RDF" -command "destroy $rdf " 
button $rdf.quitrdf.cancel -text "Cancel" -command "destroy $rdf.quitrdf"
bind $rdf.quitrdf.cancel <Return>  "destroy $rdf.quitrdf"
focus $rdf.quitrdf.cancel
pack $rdf.quitrdf.quitmessage -side top
pack $rdf.quitrdf.confirmquit  $rdf.quitrdf.cancel -side left 
pack $rdf.quitrdf -side top
}
#=====================================================================================


#=====================================================================================
#========================== printrdf procedure =======================================
# Print out rdf value or store to a file 
proc printrdf { } {
global rdf
set cutoff [expr double ([$rdf.data.cutoff.num get])]
set mbin [$rdf.data.mbin.num get ]
set dr  [expr $cutoff /  $mbin ]
global rdf 
set mbin [$rdf.data.mbin.num get ]
global rdffunc 

set filename [$rdf.file.name get]
set fileId [open $filename w 0755 ]

for { set i 0 } { $i < $mbin } { incr i 1 } {
puts $fileId "[expr $i*$dr]     $rdffunc($i)"
}    
close $fileId
}

#=======================================================================================
