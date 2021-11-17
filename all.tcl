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
set file "MD8DLETRPCbCg_autocorr.txt"
#-------------------------------------------------------------------------

#change your {fit_text} here
#-------------------------------------------
set fit_text "name N or name CA or name C"
#-------------------------------------------


#change selection here for A and B vector arrays. 
#These are used for order and auto correlation.
#-----------------------------------------------
set textA "(resname DLE or resname TRP) and name CB"
set textB "(resname DLE or resname TRP) and name CG"

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

#Because get {xyz} command in vmd is very slow, have it done once and save 
#all coordinates in all_coor(i).

	for {set i 0} {$i < $num} {incr i} {

	#go to the correct frame
        $A frame $i
	$B frame $i
        set all_coorA($i) [$A get {x y z}]
	set all_coorB($i) [$B get {x y z}]
	}

	set len [llength $all_coorA(0)]

puts "Done inputting $num frames of $len coordinates."

$sel1 delete
$selall delete
$ref1 delete
$A delete
$B delete

#k is the index for residues
	for {set k 0} {$k < $len} {incr k} {


#i and j are indices for running frames
  		for {set j 0} {$j < $numhalf} {incr j} {

		#initializing order
		set order 0.0

			for {set i 0} {$i < $numhalf} {incr i} {
			set n [expr ($i + $j)]
			#get vector coordinates

  			set ANco [lindex $all_coorA($n) $k]
			set BNco [lindex $all_coorB($n) $k]
			#normalizing difference vector
			set subN [vecnorm [vecsub $BNco $ANco]]

    	  		set AIco [lindex $all_coorA($i) $k]
			set BIco [lindex $all_coorB($i) $k]
			set subI [vecnorm [vecsub $BIco $AIco]]

			set kos1 [vecdot $subI $subN]
			set orderr [expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            		set order [expr ($order + $orderr)]
			}
			set order [expr ($order * $factor)]

#write out
			set fileId [open $file a]
			puts $fileId "$j [format "%.3f" $order]"
			close $fileId
    		}

		puts "Residue [expr (1 + $k)] done"
	}
puts "autocorr done, entering order parameter calculation..."

#------------------------------------------------------------------------
set file "MD8DLETRPCbCg_2368ps.txt"
#-------------------------------------------------------------------------


#The maximum time length for order calculation is numframes,
#but one can choose shorter time.  If this is desired, enter your length here
#-------------------------------------------
set num1 2368
#-------------------------------------------


#If num1<numframes, one can use num1 as sliding time window to sample multiple
#sets of order.  Here we do 10 sets.

set m [expr ($num - $num1)/9]

set begin 0
set factor [expr 1.0/($num1*$num1)]

for {set n 0} {$n < 10} {incr n} {


#k is the index for residues
	for {set k 0} {$k < $len} {incr k} {

		#initializing order
		set order 0.0

#i and j are indices for running frames
  		for {set j $begin} {$j < [expr ($begin+$num1)]} {incr j} {

			#get vector coordinates

  			set AJco [lindex $all_coorA($j) $k]
			set BJco [lindex $all_coorB($j) $k]
			#normalizing difference vector
			set subJ [vecnorm [vecsub $BJco $AJco]]

				for {set i $begin} {$i < [expr ($begin+$num1)]} {incr i} {
    	  				set AIco [lindex $all_coorA($i) $k]
					set BIco [lindex $all_coorB($i) $k]
					set subI [vecnorm [vecsub $BIco $AIco]]

					set kos1 [vecdot $subI $subJ]
					set orderr [expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            			set order [expr ($order + $orderr)]
				}
    		}
		set order [expr ($order*$factor)]
		#write out
set fileId [open $file a]
		puts $fileId "[expr (1 + $k)] [format "%.3f" $order]"
		close $fileId

		puts "Residue [expr (1 + $k)] done"
	}

set begin [expr ($begin + $m)]
set fileId [open $file a]
puts $fileId ""
close $fileId

if {$m < 10} {

	puts "There is no need to repeat.  Done!"
	break
	} else {
	continue
	}
}

puts "Full order done.  Entering 200 ps order calculations..."
#------------------------------------------------------------------------
set file "MD8DLETRPCbCg_order_200ps.txt"
#-------------------------------------------------------------------------


#The maximum time length for order calculation is numframes,
#but one can choose shorter time.  If this is desired, enter your length here
#-------------------------------------------
set num1 200
#-------------------------------------------


#If num1<numframes, one can use num1 as sliding time window to sample multiple
#sets of order.  Here we do 10 sets.

set m [expr ($num - $num1)/9]

set begin 0
set factor [expr 1.0/($num1*$num1)]

for {set n 0} {$n < 10} {incr n} {


#k is the index for residues
	for {set k 0} {$k < $len} {incr k} {

		#initializing order
		set order 0.0

#i and j are indices for running frames
  		for {set j $begin} {$j < [expr ($begin+$num1)]} {incr j} {

			#get vector coordinates

  			set AJco [lindex $all_coorA($j) $k]
			set BJco [lindex $all_coorB($j) $k]
			#normalizing difference vector
			set subJ [vecnorm [vecsub $BJco $AJco]]

				for {set i $begin} {$i < [expr ($begin+$num1)]} {incr i} {
    	  				set AIco [lindex $all_coorA($i) $k]
					set BIco [lindex $all_coorB($i) $k]
					set subI [vecnorm [vecsub $BIco $AIco]]

					set kos1 [vecdot $subI $subJ]
					set orderr [expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            			set order [expr ($order + $orderr)]
				}
    		}
		set order [expr ($order*$factor)]
		#write out
set fileId [open $file a]
		puts $fileId "[expr (1 + $k)] [format "%.3f" $order]"
		close $fileId

		puts "Residue [expr (1 + $k)] done"
	}

set begin [expr ($begin + $m)]
set fileId [open $file a]
puts $fileId ""
close $fileId

if {$m < 10} {

	puts "There is no need to repeat.  Done!"
	break
	} else {
	continue
	}

}

display update on


