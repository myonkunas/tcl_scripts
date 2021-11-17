# BigDCD
# Justin Gullingsrud
# vmd@ks.uiuc.edu
# 15 May 2002

# Purpose: Use this script to analyze one or more DCD files that don't fit into 
# memory.  The script will arrage for your analysis function to be called
# each time a frame is loaded, then delete the frame from memory.
# The analysis script must accept one argument; BigDCD will keep track of how
# many timesteps have been loaded and call your script with that number.
# 
# How to include this function: either source the script directory, or 
# (better) place the script in one of the directories in your auto_path 
# variable and include "package require bigdcd" in your script.
#
# Example 1: 
# This computes the center of mass for each frame in the DCD file.
#
# proc mycenter { frame } {
#   global all
#   puts "$frame: [measure center $all weight mass]"
# }
# mol load psf alanin.psf
# set all [atomselect top all]
# $all global
# bigdcd mycenter alanin.dcd
#

# Example 2:
# This computes the RMS distance between each frame in a DCD file and a 
# reference pdb file.  
#
# proc myrmsd { frame } {
#   global ref sel all
#   $all move [measure fit $sel $ref]
#   puts "$frame: [measure rmsd $sel $ref]"
# }
# mol load psf protein.psf
# set all [atomselect top all]
# set ref [atomselect top "name CA" frame 0]
# set sel [atomselect top "name CA"]
# animate read pdb protein.pdb
# bigdcd myrmsd eq01.dcd  eq02.dcd eq03.dcd

proc bigdcd { script args } {
  global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
  
  set bigdcd_frame 0
  set bigdcd_firstframe [molinfo top get numframes]
  set bigdcd_proc $script
  
  uplevel #0 trace variable vmd_frame w bigdcd_callback
  foreach dcd $args {
    animate read dcd $dcd waitfor 0
  }
}

proc bigdcd_callback { name1 name2 op } {
  global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
 
  # If we're out of frames, we're also done
  set thisframe $vmd_frame($name2)
  if { $thisframe < $bigdcd_firstframe } {
    bigdcd_done
    return
  }
 
  incr bigdcd_frame
  if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } { 
    puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
    bigdcd_done
    return
  }
  animate delete beg $thisframe end $thisframe 
  return $msg
}

proc bigdcd_done { } {
  puts "bigdcd_done"
  after idle uplevel #0 trace vdelete vmd_frame w bigdcd_callback
}
  
