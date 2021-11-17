### This script contains many scripts to analyze the simulation of gA channel.
### Procedures:
#   rmsf startframe endframe  filename :run rmsf using center-of-mass fitting.
#   rmsd pdbmolId filename : measure the rmsd of top molecule comparing to pdb
#   S2 filename : calculate the S2, some parameters may need be changed
#   autocorr file : calculate autocorrelation function

########################################################################
################## rmsf using center-of-mass fitting ###################
#This script fits all frames in a simulation trajectory to the 1st frame
#according to the {fit_text},
#calculates an average structure of all frames, and then calculate the 
#root-mean-square fluctuation for atoms defined by {fluc_text}
#
set num 0
set all_coorA(0) 0
set all_coorB(0) 0



proc rmsf { startframe endframe  filename } {
	
#------------------------------------------------------------------------
set file $filename 
#-------------------------------------------------------------------------

#Replace the fit_text and fluc_text to what you want to select

set fit_text "name CA or name N or name C"
set fluc_text "name CA"


set selall [atomselect top all]
set num_frame [molinfo top get numframes]

# change $endframe,so a correct frame count could be obtained.
if { $endframe > $num_frame }  {
	set endframe $num_frame
} 
 

set frame_count [expr $endframe - $startframe ]
puts "TOTAL $frame_count frames "
#align all frames with respect to the 1st frame to remove translation and rotation

set sel1 [atomselect top $fit_text]
# First frame as reference, not useful when using center of mass fitting.
set ref1 [atomselect top $fit_text frame 0]

# Use center_of_mass to translate the frames
puts "translate according to center of mass"
set protein [atomselect top "segname SEG1 SEG2"]
# move each frame according to center of mass
for { set frame $startframe } { $frame < $endframe } { incr frame } {
	$protein frame $frame
	$selall frame $frame 
	set trans_mat [vecscale -1.0 [measure center $protein weight mass]]
	$selall moveby $trans_mat
}

### measure fit will dampen the fluctuation
#        for {set frame $startframe } {$frame < $endframe } {incr frame} {
#                $sel1 frame $frame
#                $selall frame $frame
#                set trans_mat [measure fit $sel1 $ref1]
#                $selall move $trans_mat
#        }

#Calculating average structure
    set fluc [atomselect top $fluc_text]
    set factor [expr 1./$frame_count]

	for {set i $startframe } {$i < $endframe } {incr i} {
	$fluc frame $i
	set all_coor($i) [$fluc get {x y z}]
	}

 
# how many atoms (set of coordinates ) in each frame
    	set len [llength $all_coor($startframe)]

    		for {set j 0} {$j<$len} {incr j} {
			set b [veczero]
			for {set i $startframe } {$i < $endframe } {incr i} {
	    		set b [vecadd $b [lindex $all_coor($i) $j]]
			}
			set ave_vec [vecscale $b $factor]
			puts "$j  $ave_vec "
			set sum 0.0
			for {set i $startframe } {$i < $endframe } {incr i} {
	    		set sum [expr ($sum + [veclength2 [vecsub [lindex $all_coor($i) $j] $ave_vec]])]
			}
			set sum [expr (sqrt ($sum*$factor))]
			set fileId [open $file a]
			puts $fileId "[expr ($j+1)] [format "%.3f" $sum]"
			close $fileId
puts "residue $j $sum"
    		}

}
############# End of rmsf procedure  #######################################
############################################################################

############################################################################
################# RMSD calculation #########################################
####Calculate the RMSD of c-alpha, comparing to the initial docked structure, eg#   a pdb file or the first frame of any dcd file
#Aligns the molecule before computing its RMSD; 
#Prints the RMSD of the protein atoms between each timestep.

#sel1, ref1,sel3 and ref3 are for alignment only!!!!!!!!!!!

# calculation assume dcd file to be calculated is in the top, and the molId of 
# reference pdb file was input as a parameter.  If the same molId as top mole
# was input, the first frame will be used

proc rmsd  { pdbmolId filename } { 

# Note: in the selection, must use segname, because in DMPC lipid there is also
# an atom called "N" 
set sel1 [atomselect top "{segname SEG1} and {name CA N O C} and {resid 1 to 15}"]
 set sel2 [atomselect top "{segname SEG2} and {name CA N O C} and {resid 1 to 15}"]

  #use a pdb file for the reference
 set ref1 [atomselect $pdbmolId "{segname SEG1} and {name CA N O C} and {resid 1 to 15}" frame 0]
 set ref2 [atomselect $pdbmolId "{segname SEG2} and {name CA N O C} and {resid 1 to 15} " frame 0]


#align all frames with respect to reference to remove translation and rotation

# open the file for writing
        set fileId [open $filename a]
        set num [molinfo top get numframes]
        for {set frame 0} {$frame < $num} {incr frame} {
                $sel1 frame $frame
                $sel2 frame $frame
                set trans_mat1 [measure fit $sel1 $ref1]
                set trans_mat2 [measure fit $sel2 $ref2]
                $sel1 move $trans_mat1
                $sel2 move $trans_mat2

        #compute the RMSD
        set rmsd1 [measure rmsd $sel1 $ref1]
        set rmsd2 [measure rmsd $sel2 $ref2]

        puts $fileId "$frame [format "%2.3f" $rmsd1] [format "%2.3f" $rmsd2]"
	}
#" close fileId to flush
        close $fileId 
}

################### End of rmsd procedure ###############################
#########################################################################

#########################################################################
#################### Vector Loader, a procedure called by "auto" or "S2" 
#################### To extract the coordinates. 
proc VectorLoader { } { 
#Script to load vectors for auto.tcl and/or S2.tcl calculations
#Verion 03052001, by Yan Xu
#Because of memory leak in VMD, coordinates should be loaded into arrays to 
#save memory and time
#-------------------------------------------------------------------------

#change your {fit_text} here for structual fitting
global num
global all_coorA
global all_coorB
#-------------------------------------------
set fit_text "segname SEG1 SEG2 and name N CA C "
#-------------------------------------------

#change selection here for A and B vector arrays. 
#These are used for order and auto correlation.
#-----------------------------------------------
set textA "protein and (segname SEG1 or segname SEG2) and name N"
set textB "protein and (segname SEG1 or segname SEG2) and name HN"

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

display update on
}
 
####################### End of VectorLoader ######################
##################################################################


##################################################################
################ S2 calculation ##################################
#Script to calculate order parameter
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
proc  S2 { filename } { 
run VectorLoader
global num
global all_coorA
global all_coorB
#------------------------------------------------------------------------
set file $filename
#-------------------------------------------------------------------------


#The maximum time length for order calculation is numframes,
#but one can choose shorter time.  If this is desired, enter your length here
#-------------------------------------------
set num1 1000
#-------------------------------------------

display update off

#If num1<numframes, one can use num1 as sliding time window to sample multiple
#sets of order.  Here we do 10 sets.
# Zhanwu , Jul 27, 8 iterations 
set m [expr ($num - $num1)/7]

set begin 0
set factor [expr 1.0/($num1*$num1)]

for {set n 0} {$n < 8 } {incr n} {


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
}
######################## End of S2 procedure #######################
####################################################################


####################################################################
##################### Auto correlation function ####################
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

proc autocorr { filename }{ 
#------------------------------------------------------------------------
set file $filename 
#-------------------------------------------------------------------------

#change your {fit_text} here
#-------------------------------------------
set fit_text "segname SEG1 SEG2 and name N CA C"
#-------------------------------------------

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

#only consider residue 10 of seg 1 and res15 of seg2, so directly use the atom number

  set A [atomselect top "protein and (segname SEG1 or segname SEG2) and name N"]
#set A [atomselect top "index 137 518"]
#set B [atomselect top "index 138 519"]
  set B [atomselect top "protein and (segname SEG1 or segname SEG2) and name HN"]

#Because get {xyz} command in vmd is very slow, have it done once and save 
#all coordinates in all_coor(i).

	for {set i 0} {$i < $num} {incr i} {

	#go to the correct frame
      $A frame $i
	$B frame $i
      set all_coorA($i) [$A get {x y z}]
	set all_coorB($i) [$B get {x y z}]
	}

$selall delete

	set len [llength $all_coorA(0)]

puts "Done inputting $num frames of $len coordinates."


#k is the index for residues
	for {set k 0} {$k < $len} {incr k} {
                        set fileId [open $file a]
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

#write out
			set fileId [open $file a]
			puts $fileId "$j [format "%.3f" $order]"
			close $fileId
    		}

		puts "Residue [expr (1 + $k)] done"
	}

display update on
}
####################### End of autocorr #####################################
#############################################################################

