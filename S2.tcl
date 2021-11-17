#!!!!!RUN VectorLoader.tcl before running this script!!!!!!
#!!!!!RUN VectorLoader.tcl before running this script!!!!!!
#!!!!!RUN VectorLoader.tcl before running this script!!!!!!
#!!!!!RUN VectorLoader.tcl before running this script!!!!!!
#!!!!!RUN VectorLoader.tcl before running this script!!!!!!
#
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
#------------------------------------------------------------------------
set file "control_MD_order8ns.txt"
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

