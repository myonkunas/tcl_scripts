#Script to calculate autocorrelation function
#Verion 03052001, by Yan Xu
#
#To calculate autocorrelation of vectors over simulation time
#Let a unit vector at time i be Ui and the same vector at time j be Uj
#then cos(theta_ij)=Ui*Uj and autocorrelation over time is 
#    Ci(T/2)={SUM(i=1 to T/2)[3cos(theta_ij)-1]/2}/[T/2], where T is total
# simulation time.

#Theoretically, the order parameter = S2=lim(T to infinite)Ci(T/2)
#but when T is finite, it can be estimated as
# S2 = {SUM(i=1toT)SUM(j=1toT)[3cos(theta_ij)-1]/2}/T/T
#------------------------------------------------------------------------
#-------------------------------------------------------------------------

#change your {fit_text} here
#-------------------------------------------
set fit_text "name P or name O11 or name O12"
#-------------------------------------------


#change selection here for A, B, C, and D vector arrays. 
#These are used for order and auto correlation.
#-----------------------------------------------
#set textA "name C21 or name C22 or name C31 or name C32"
#set textB "name C23 or name C24 or name C33 or name C34"
# first half is interfacial, second is not
#set interf_resid " resid 25 31 34 35 36 46 48 66 68 69 0 3 8 40 42 49 54 57 59 60 65"
set textA "name C21 or name C22 or name C23 or name C24 or name C25 or name C26 or name C27 or name C28 or name C29 or name C210 or name C211 or name C212 "
set textB "name C23 or name C24 or name C25 or name C26 or name C27 or name C28 or name C29 or name C210 or name C211 or name C212 or name C213 or name C214 "
#set textC "name C31 or name C32 or name C33 or name C34 or name C35 or name C36 or name C37 or name C38 or name C39 or name C310 or name C311 or name C312"
#set textD "name C33 or name C34 or name C35 or name C36 or name C37 or name C38 or name C39 or name C310 or name C311 or name C312 or name C313 or name C314"


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
set numhalf [expr ($num / 2)]
set factor [expr (1.0 / $numhalf)]

	for {set frame 0} {$frame < $num} {incr frame} {
		$sel1 frame $frame
		$selall frame $frame
		set trans_mat [measure fit $sel1 $ref1]
		$selall move $trans_mat
	}
puts "Structural alignment done"


  set A [atomselect top $textA]
  set B [atomselect top $textB]
#  set C [atomselect top $textC]
#  set D [atomselect top $textD]

#Because get {xyz} command in vmd is very slow, have it done once and save 
#all coordinates in all_coor(i).

	for {set i 0} {$i < $num} {incr i} {

	#go to the correct frame
        $A frame $i
	$B frame $i
	#$C frame $i
	#$D frame $i
        set all_coorA($i) [$A get {x y z}]
	set all_coorB($i) [$B get {x y z}]
        #set all_coorC($i) [$C get {x y z}]
	#set all_coorD($i) [$D get {x y z}]
	}

	set len [llength $all_coorA(0)]

puts "Done inputting $num frames of $len coordinates."

$sel1 delete
$selall delete
$ref1 delete
$A delete
$B delete
#$C delete
#$D delete


#------------------------------------------------------------------------
set file1 "PFE_interfacial_lip_order.dat"
#set file2 "C:MD100InterfaceLIP_S2_chain2.txt"
set fileId [open $file1 a]
puts $fileId "Data in the sequence of: C21-C23, C22-C24, ... and C212-C214"
close $fileId
#set fileId [open $file2 a]
#puts $fileId "Data in the sequence of: C33-C35, C34-C36 ... C312-C314"
#close $fileId
#-------------------------------------------------------------------------


#The maximum time length for order calculation is numframes,
#but one can choose shorter time.  If this is desired, enter your length here
#-------------------------------------------

#-------------------------------------------

display update off



set factor [expr 1.0/($num*$num)]



#k is the index for residues
	for {set k 38} {$k < $len} {incr k} {

		#initializing order
		set order1 0.0
		#set order2 0.0

#i and j are indices for running frames
  		for {set j 0} {$j < $num} {incr j} {

			#get vector coordinates

  			set AJco [lindex $all_coorA($j) $k]
			set BJco [lindex $all_coorB($j) $k]
  			#set CJco [lindex $all_coorC($j) $k]
			#set DJco [lindex $all_coorD($j) $k]

			#normalizing difference vector
			set subJ1 [vecnorm [vecsub $BJco $AJco]]
			#set subJ2 [vecnorm [vecsub $DJco $CJco]]


				for {set i 0} {$i < $num} {incr i} {
    	  				set AIco [lindex $all_coorA($i) $k]
					set BIco [lindex $all_coorB($i) $k]
					set subI1 [vecnorm [vecsub $BIco $AIco]]
   	  				#set CIco [lindex $all_coorC($i) $k]
					#set DIco [lindex $all_coorD($i) $k]
					#set subI2 [vecnorm [vecsub $DIco $CIco]]

					set kos1 [vecdot $subI1 $subJ1]
					set orderr1 [expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            			set order1 [expr ($order1 + $orderr1)]
					#set kos2 [vecdot $subI2 $subJ2]
					#set orderr2 [expr (0.5 * (3.0 * $kos2 * $kos2 - 1.0)) ]
            			#set order2 [expr ($order2 + $orderr2)]
				}
    		}
		set order1 [expr ($order1*$factor)]
		#set order2 [expr ($order2*$factor)]
		#write out
		set fileId [open $file1 a]
		puts $fileId "[expr (1 + $k)] [format "%.3f" $order1]"
		close $fileId
		#set fileId [open $file2 a]
		#puts $fileId "[expr (1 + $k)] [format "%.3f" $order2]"
		#close $fileId

		puts "Lipid [expr (1 + $k)] done"
	}







display update on


