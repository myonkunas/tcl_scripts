#To calculate autocorrelation of vectors over simulation time
#Let a unit vector at time i be Ui and the same vector at time j be Uj
#then cos(theta_ij)=Ui*Uj and autocorrelation over time is 
#    Ci(T/2)={SUM(i=1 to T/2)[3cos(theta_ij)-1]/2}/[T/2], where T is total
# simulation time.

#Theoretically, the order parameter = S2=lim(T to infinite)Ci(T/2)
#but when T is finite, it can be estimated as
# S2 = {SUM(i=1toT)SUM(j=1toT)[3cos(theta_ij)-1]/2}/T/T
#------------------------------------------------------------------------
set file "~/GRAMICIDIN/ORDER_CAL/MD3SEG1_0_200.txt"
#-------------------------------------------------------------------------

#change your {fit_text} here
set fit_text "name N or name CA or name C"

#sel1 and ref1 are for alignment only!!!!!!!!!!!
#turn off display for speed

display update off

#assuming the molecule of interest is on the top
#input psf and dcd files to vmd, select atoms 
  set sel1 [atomselect top $fit_text]

#use the first frame as the reference
  set ref1 [atomselect top $fit_text frame 0]

set selall [atomselect top all]

#align all frames with respect to reference to remove translation and rotation

	set num [molinfo top get numframes]
	for {set frame 0} {$frame < $num} {incr frame} {
		$sel1 frame $frame
		$selall frame $frame
		set trans_mat [measure fit $sel1 $ref1]
		$selall move $trans_mat
	}
puts "Structural alignment done"
#Loop for all residues

set factor [expr 1.0/(200*200)]
#set factor [expr 1.0/($num*$num)]
for {set residue 1} {$residue <= 1} {incr residue} {

  set AI [atomselect top "segname SEG1 and name N and resid $residue"]
  set BI [atomselect top "segname SEG1 and name HN and resid $residue"]
  set AJ [atomselect top "segname SEG1 and name N and resid $residue"]
  set BJ [atomselect top "segname SEG1 and name HN and resid $residue"]


#initializing order
set order 0

  for {set j 0} {$j < 200} {incr j} {

puts "Entering residue $residue loop $j"
	#go to the correct frame
	$AJ frame $j
	$BJ frame $j
	#get vector coordinates
  	#set AJco [lindex [$AJ get {x y z}] 0]
	#set BJco [lindex [$BJ get {x y z}] 0]
	set AJco [$AJ get {x y z}]
	set BJco [$BJ get {x y z}]
	#normalizing difference vector
	#set subJ [vecnorm [vecsub $BJco $AJco]]

		for {set i 0} {$i < 200} {incr i} {
    	
		$AI frame $i
		$BI frame $i

		set AIco [$AI get {x y z}]
		set BIco [$BI get {x y z}]

  		#set AIco [lindex [$AI get {x y z}] 0]
		#set BIco [lindex [$BI get {x y z}] 0]
		#set subI [vecnorm [vecsub $BIco $AIco]]

		#set kos1 [vecdot $subI $subJ]
		#set order1 [expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            	#set order [expr ($order + $order1)]
		}
	}
	#set order [expr (0.0 + $order*$factor)]

	#set fileId [open $file a]
	#puts $fileId "$residue [format "%.3f" $order]"
	#close $fileId
}

display update on

