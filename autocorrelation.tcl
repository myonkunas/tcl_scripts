#To calculate autocorrelation of vectors over simulation time
#Let a unit vector at time i be Ui and the same vector at time j be Uj
#then cos(theta_ij)=Ui*Uj and autocorrelation over time is 
#    Ci(T/2)={SUM(i=1 to T/2)[3cos(theta_ij)-1]/2}/[T/2], where T is total
# simulation time.


#Theoretically, the order parameter = S2=lim(T to infinite)Ci(T/2)

#but when T is finite, it can be estimated as
# S2 = {SUM(i=1toT)SUM(j=1toT)[3cos(theta_ij)-1]/2}/T/T

#sel1, ref1,sel3 and ref3 are for alignment only!!!!!!!!!!!

#input psf and dcd files to vmd, select atoms
  set sel1 [atomselect 0 "segname SEG1 and name N and resid 1 to 15"]

#use a pdb file for the reference
  set ref1 [atomselect 1 "segname SEG1 and name N and resid 1 to 15"]

# N-H vectors defined here
# set rn1 [[atomselect 1 "segname SEG1 and name N and resid 1"] get {x y z}]
#  set rn1 [lvarpop rn1]
#  set rh1 [[atomselect 1 "segname SEG1 and name HN and resid 1"] get {x y z}]
#  set rh1 [lvarpop rh1]
  set SA1 [atomselect 0 "segname SEG1 and name N and resid 1"]
  set SB1 [atomselect 0 "segname SEG1 and name HN and resid 1"]
#  set diff1 [vecsub $rh1 $rn1]
#  set D1 [veclength $diff1]

set num [molinfo 0 get numframes]
set order 0

#set file "/monet/people/tang/GRAMICIDIN/MDCALC3/PROCESS/order.txt"
set file "order.txt"
set fileId [open $file a]

for {set j 1} {$j < $num/2} {incr j} {
   for {set i 0} {$i < $num/2} {incr i} {
	#get the correct frame
	$sel1 frame [expr $i+$j]

	#compute the transformation
	set trans_mat2 [measure fit $sel1 $ref1]

  #get the correct frame and do the alignment
	$SA1 frame [expr $i+$j]
	$SA1 move $trans_mat2
	$SB1 frame [expr $i+$j]
	$SB1 move $trans_mat2

	set JA1 [$SA1 get {x y z}]
	set JB1 [$SB1 get {x y z}]
        set JA1 [lvarpop JA1]
        set JB1 [lvarpop JB1]
	set diffJ [vecsub $JB1 $JA1]
        set lengthJ [veclength $diffJ]

      $sel1 frame $i

	#compute the transformation
	set trans_mat1 [measure fit $sel1 $ref1]

	#get the correct frame and do the alignment
	$SA1 frame $i
	$SA1 move $trans_mat1
	$SB1 frame $i
	$SB1 move $trans_mat1

        #compute the vector difference

	set IA1 [$SA1 get {x y z}]
	set IB1 [$SB1 get {x y z}]
        set IA1 [lvarpop IA1]
        set IB1 [lvarpop IB1]
		set diffI [vecsub $IB1 $IA1]
	   	set lengthI [veclength $diffI]
 	   	set dotIJ [vecdot $diffI $diffJ]

            set kos1 [expr $dotIJ / ($lengthI * $lengthJ)]
   		set order1 [ expr (0.5 * (3.0 * $kos1 * $kos1 - 1.0)) ]
            set order [expr $order + $order1]
        puts  "frame $i $order1"
   }

        set order [expr 0.0 + $order /($num*0.5)]
	puts $fileId "$order"
}
	close $fileId
