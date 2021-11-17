#To calculate autocorrelation of vectors over simulation time
#Let a unit vector at time i be Ui and the same vector at time j be Uj
#then cos(theta_ij)=Ui*Uj and autocorrelation over time is 
#    Ci(T/2)={SUM(i=1 to T/2)[3cos(theta_ij)-1]/2}/[T/2], where T is total
# simulation time.

#Theoretically, the order parameter = S2=lim(T to infinite)Ci(T/2)
#but when T is finite, it can be estimated as
# S2 = {SUM(i=1toT)SUM(j=1toT)[3cos(theta_ij)-1]/2}/T/T
#------------------------------------------------------------------------
set file "/monet/people/tang/GRAMICIDIN/MDCALC3/PROCESS/MD3SEG1_10_20order.txt"
#set fileId [open $file a]
#-------------------------------------------------------------------------
#sel1, ref1,sel3 and ref3 are for alignment only!!!!!!!!!!!

#input psf and dcd files to vmd, select atoms
  set sel1 [atomselect 0 "{segname SEG1} and {name N} and {resid 1 to 15}"]

#use a pdb file for the reference
  set ref1 [atomselect 1 "{segname SEG1} and {name N} and {resid 1 to 15}"]

set selall [atomselect 0 all]

#align all frames with respect to reference to remove translation and rotation


	set num [molinfo 0 get numframes]
	for {set frame 0} {$frame < $num} {incr frame} {
		$sel1 frame $frame
		$selall frame $frame
		set trans_mat [measure fit $sel1 $ref1]
		$selall move $trans_mat
	}
puts "Structural alignment done"
#Loop for all residues

for {set residue 1} {$residue <= 2} {incr residue} {

  set SA1 [atomselect 0 "{segname SEG1} and {name N} and {resid $residue}"]
  set SB1 [atomselect 0 "{segname SEG1} and {name HN} and {resid $residue}"]

#initializing order
set order 0

  for {set j 0} {$j < 2} {incr j} {
puts "Entering residue $residue loop $j"
    	#get the correct frame
	$SA1 frame $j
	$SB1 frame $j
	#get coordinates
	set JA1 [$SA1 get {x y z}]
	set JB1 [$SB1 get {x y z}]
	#assign vectors
#        set JA1 [lvarpop JA1]
#        set JB1 [lvarpop JB1]
	#cal culate vector from N-H
#	set diffJ [vecsub $JB1 $JA1]
#       set lengthJ [veclength $diffJ]
#set JA1 [lvarpop [$SA1 get {x y z}]]
#set JB1 [lvarpop [$SB1 get {x y z}]]
# set lengthJ [veclength [vecsub $JB1 $JA1]] 

		for {set i 0} {$i < 2} {incr i} {

		#get the correct frame
		$SA1 frame $i
		$SB1 frame $i

        	#compute the vector difference

		set IA1 [$SA1 get {x y z}]
		set IB1 [$SB1 get {x y z}]
#        	set IA1 [lvarpop IA1]
#       	set IB1 [lvarpop IB1]

#		set diffI [vecsub $IB1 $IA1]
#	   	set lengthI [veclength $diffI]

#set IA1 [lvarpop [$SA1 get {x y z}]]
#set IB1 [lvarpop [$SB1 get {x y z}]]
#set lengthI [veclength [vecsub $IB1 $IA1]] 

#set dotIJ [vecdot{vecsub $IB1 $IA1}{vecsub $JB1 $JA1}]
set kos1 [expr ([vecdot[vecsub $IB1 $IA1][vecsub $JB1 $JA1]] / ([veclength [vecsub $IB1 $IA1]] *[veclength [vecsub $JB1 $JA1]] ))]

# 	   	set dotIJ [vecdot $diffI $diffJ]
#            	set kos1 [expr ($dotIJ / ($lengthI * $lengthJ))]
   		set order1 [ expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            	set order [expr ($order + $order1)]
		#puts  "framei $i $order1"
   		}
#puts  "frame $j $order"
	}
#	set order [expr (0.0 + $order /($num*$num))]
	set order [expr (0.0 + $order /(1*1))]
#puts  "residue $residue $order"

	set fileId [open $file a]
	puts $fileId "$residue [format "%.3f" $order]"
	close $fileId
}


