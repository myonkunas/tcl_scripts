# This script come from 
# http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/1549.html
# by From: Justin Gullingsrud (justin_at_ks.uiuc.edu)
# Date: Wed Jun 04 2003 - 12:29:43 CDT

# Here's a little script in the same format as Jan's script that does the
# same thing using the label command. It uses the current timestep by
# default, instead of the first timestep. It also blows away any existing
# Dihedral labels; this is because of a deficiency in the scripting
# interface.

# calculate the X1 and X2 dihedral angle of gA, input segment name
proc ga_trp_x { } {
    foreach segment { SEG1 } {
	foreach res { 9 11 13 15 } {
		set a1 [[atomselect top "segname $segment and resid $res and name N "] list] 
		set a2 [[atomselect top "segname $segment and resid $res and name CA "] list]
		set a3 [[atomselect top "segname $segment and resid $res and name CB "] list]
		set a4 [[atomselect top "segname $segment and resid $res and name CG"] list]
		set filename "$segment\_$res\_X1.txt"
		all_dihed_angle $a1 $a2 $a3 $a4 $filename
	}
	foreach res { 9 11 13 15 } {
		set a1 [[atomselect top "segname $segment and resid $res and name CA"] list] 
		set a2 [[atomselect top "segname $segment and resid $res and name CB "] list]
		set a3 [[atomselect top "segname $segment and resid $res and name CG "] list]
		set a4 [[atomselect top "segname $segment and resid $res and name CD1"] list]
		set filename "$segment\_$res\_X2.txt"
		all_dihed_angle $a1 $a2 $a3 $a4 $filename
	}
   }
}



proc ga_8_10_x { } {
    foreach segment { SEG1 } {
	foreach res {  10 } {
		set a1 [[atomselect top "segname $segment and resid $res and name N "] list] 
		set a2 [[atomselect top "segname $segment and resid $res and name CA "] list]
		set a3 [[atomselect top "segname $segment and resid $res and name CB "] list]
		set a4 [[atomselect top "segname $segment and resid $res and name CG"] list]
		set filename "$segment\_$res\_X1.txt"
		all_dihed_angle $a1 $a2 $a3 $a4 $filename
puts "here"
	}
	foreach res {  10 } {
		set a1 [[atomselect top "segname $segment and resid $res and name CA"] list] 
		set a2 [[atomselect top "segname $segment and resid $res and name CB "] list]
		set a3 [[atomselect top "segname $segment and resid $res and name CG "] list]
		set a4 [[atomselect top "segname $segment and resid $res and name CD2"] list]
		set filename "$segment\_$res\_X2.txt"
		all_dihed_angle $a1 $a2 $a3 $a4 $filename
	}
   }
}


proc my_dihed_angle { a1 a2 a3 a4 { frame now }} {
  # Delete all existing Dihdral labels so that the one we add has index 0.
  label delete Dihedrals

  # Use the top molecule
  set molid [molinfo top]

  # Override the current frame
  set curframe [molinfo top get frame]
  molinfo top set frame $frame

  # Add the label
  label add Dihedrals $molid/$a1 $molid/$a2 $molid/$a3 $molid/$a4

  # Get its value
  set curval [lindex [lindex [label list Dihedrals] 0] 4]

  # Restore the old frame
  molinfo top set frame $curframe

  return $curval
}

# Note that if you just want the value of the dihedral for all timesteps,
# there's a much faster way to do it: the "label graph" command. Here's
# a script that returns the values for all timesteps:

proc all_dihed_angle { a1 a2 a3 a4  filename } {
  # Delete all existing Dihdral labels so that the one we add has index 0.
  label delete Dihedrals

  # Use the top molecule
  set molid [molinfo top]

  # Add the dihedral monitor
  label add Dihedrals $molid/$a1 $molid/$a2 $molid/$a3 $molid/$a4

  # Get all values
  # and write to a file, added by Zhanwu
  return [label graph Dihedrals 0  $filename]
}

#Cheers,
#Justin 
