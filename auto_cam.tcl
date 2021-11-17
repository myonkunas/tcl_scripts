#Script to calculate autocorrelation function
#Verion 03052001, by Yan Xu
#Modified 10/21/04 Yonkunas
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

proc auto {residue1 residue2 file} {

#change your {fit_text} here
#-------------------------------------------
set fit_text "name N or name CA or name C"
#-------------------------------------------

#sel1 and ref1 are for alignment only!!!!!!!!!!!
#turn off display for speed

#display update off

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
set all [atomselect top all]

set A [atomselect top "protein and (resid $residue1 to $residue2) and name N"]
set B [atomselect top "protein and (resid $residue1 to $residue2) and name HN"]

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
$all delete
puts "Done inputting $num frames of $len coordinates."


#k is the index for residues
	for {set k 0} {$k < $len} {incr k} {
                        set fileId [open $file.$k.txt a]
                        puts $fileId " " 
                        close $fileId


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


			set fileId [open $file.$k.txt a]
			puts $fileId "$j [format "%2.3f" $order]"
			close $fileId
		
    		}

		puts "Residue [expr (1 + $k)] done"
	}
	}
#display update on

