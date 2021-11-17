set totalframenum [molinfo top get numframes]
set fileId [open OaroundNa-11.txt a ]

set oxygen [atomselect top {name "O.*" and within 3.5 of segname ION}]
set sodium [atomselect top "segname ION"]

for {set framenum 0 } { $framenum < $totalframenum } { incr framenum } { 
	$oxygen frame $framenum
	$sodium frame $framenum
	$oxygen update 
	
	puts $fileId "[$sodium get z] [$oxygen get {segname resid name}]"
}

close $fileId
