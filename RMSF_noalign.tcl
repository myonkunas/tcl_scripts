#To determine the structural stability, 
#the RMSF of Ca atoms from their average values over the duration of
#the simulation can be calculated using following procedures.
#By selecting Starting and End, the domains involved in the calculation may vary.
#The averaged strucure of pentermer over simulation time was calculated
#by "gOpenMol"


set num [molinfo top get numframes]


#------------ output filename --------------------------------------------------
set file "PFE_rmsf_noalign_6-8s.txt"
#---------------------------------------------------------------------------------

#sel1, ref1,sel3 and ref3 are for alignment only!!!!!!!!!!!

#input psf and dcd files to vmd, select atoms
  set sel1 [atomselect top "{segname SEG1} or {segname SEG2} and {name CA} and {resid 1 to 15}"]

#use a pdb file for the reference
  set ref1 [atomselect 0 "{segname SEG1} or {segname SEG2} and {name CA} and {resid 1 to 15}"]
 
 set selall [atomselect top all]

#align all frames with respect to reference to remove translation and rotation
#
#	set num [molinfo 0 get numframes]
#	for {set frame 0} {$frame < $num} {incr frame} {
#		$sel1 frame $frame
#		$selall frame $frame
#		set trans_mat [measure fit $sel1 $ref1]
#		$selall move $trans_mat
#	}
puts "Structural alignment done"


#------- starting residue----
set start 1
#---------------------------

#------- end residue---------
set end 16
#-----------------------------



for {set residue $start} {$residue < $end}  {incr residue 1} {
	set rmsfSEG1 0
	set rmsfSEG2 0

	set ref1 [atomselect 0 "{segname SEG1} and {name CA} and {resid $residue}"]
	set xref1 [$ref1 get {x}]
	set yref1 [$ref1 get {y}]
	set zref1 [$ref1 get {z}]

	set ref2 [atomselect 0 "{segname SEG2} and {name CA} and {resid $residue}"]
	set xref2 [$ref2 get {x}]
	set yref2 [$ref2 get {y}]
	set zref2 [$ref2 get {z}]
	
	#--------------------------- selection -----------------------------
	 set sel1 [atomselect top "{segname SEG1} and {name CA} and {resid $residue}"]
	 set sel2 [atomselect top "{segname SEG2} and {name CA} and {resid $residue}"]
	#-------------------------------------------------------------------

	for {set frame 0} {$frame < $num} {incr frame} {       
	
		$sel1 frame $frame
		set xs1 [$sel1 get {x}]
		set ys1 [$sel1 get {y}]
		set zs1 [$sel1 get {z}]
	set dR1 [expr ($xs1 - $xref1) * ($xs1 - $xref1) + ($ys1 - $yref1) * ($ys1 - $yref1) + ($zs1 - $zref1) * ($zs1 - $zref1)]
	set rmsfSEG1 [expr $rmsfSEG1 + $dR1]
				
		$sel2 frame $frame
		set xs2 [$sel2 get {x}]
		set ys2 [$sel2 get {y}]
		set zs2 [$sel2 get {z}]
	set dR2 [expr ($xs2 - $xref2) * ($xs2 - $xref2) + ($ys2 - $yref2) * ($ys2 - $yref2) + ($zs2 - $zref2) * ($zs2 - $zref2)]
	set rmsfSEG2 [expr $rmsfSEG2 + $dR2] 	
 
	}

	set rmsfSEG1 [expr sqrt( $rmsfSEG1 / $num  ) ]
	set rmsfSEG2 [expr sqrt( $rmsfSEG2 / $num  ) ]
	

        puts  "residue $residue [format "%.3f" $rmsfSEG1] [format "%.3f" $rmsfSEG2]"

	set fileId [open $file a]
	puts $fileId "$residue [format "%.3f" $rmsfSEG1] [format "%.3f" $rmsfSEG2]"
	close $fileId
}
