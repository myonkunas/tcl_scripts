# hole: a script for running the HOLE program written by Oliver Smart
# URL: http://www.bip.bham.ac.uk/osmart/hole/top.html
#
# Here's example usage:
#
#  source hole.tcl
#  mol load pdb pull06-pore-0ns.pdb
#  set sel [atomselect top all]
#  Hole::runhole $sel
#
# Justin Gullingsrud
# jgulling@mccammon.ucsd.edu
# 25 October 2003


package provide hole 1.0

namespace eval Hole {
  # Customize the following lines to set the paths to the hole executable and
  # the radius file.
  variable holebin /usr/local/HOLE/exe_linux/
  variable holerad /usr/local/HOLE/rad/simple.rad

  # Customize the following to set default values.  See the runhole 
  # comments for what these parameters do.
  variable cvect [list 0. 0. 1.]
  variable cpoint [list 0. 0. 0.]
  variable sample 0.05
  variable endrad 10.
}

# Routine for calling hole.  Pass molid, frame, and optional keyword
# arguments:
#   -cvect {x y z}    vector parallel to channel axis
#   -cpoint {x y z}   a point in the pore
#   -sample s         distance between samples
# Returns a list whose elements are of the form {z r resname resid},
# where z is the channel coordinate, r is the radius at that coordinate,
# and resname and resid give the identity of the residue nearest the center
# of the pore at that coordinate.
proc Hole::runhole {holeoutput sel args} {
  variable holebin
  variable holerad
  variable cvect
  variable cpoint
  variable sample
  variable endrad

  # parse options
  foreach { opt val } $args {
    switch $opt {
      -cvect { 
        if { [llength $cvect] != 3 } {
          error "-cvect must have three elements"
        }
        set cvect $val 
      }
      -cpoint { 
        if { [llength $cpoint] != 3 } {
          error "-cpoint must have three elements"
        }
        set cpoint $val 
      }
      -sample { set sample $val }
      -endrad { set endrad $val }
      default {
        error "Unknown option '$opt'"
      }
    }
  }
  if { [$sel num] < 3 } {
    error "Not enough atoms ([$sel num]) found in selection."
  }
  # write coordinates to files
  set pdb tmpholeinputfiles.pdb
  $sel writepdb $pdb

  # construct HOLE input string
  set str "\ncoord $pdb\n"
  append str "radius $holerad\n"
  append str "cvect $cvect\n"
  append str "cpoint $cpoint\n"
  append str "sample $sample\n"
  append str "endrad $endrad\n"
  append str "SPHPDB $holeoutput.sph\n"

  # Call HOLE and collect output
  puts "Calling HOLE..."
  flush stdout
  set result [exec "$holebin\hole" "<< $str"]


  set rawdata_z [list]
  set rawdata_r [list]
  set rawdata_resname [list]
  set rawdata_resid [list]


  set fileid [open $holeoutput.xls a+ 0777]

  set lines [split $result \n]
  set n [llength $lines]
  for { set i 0 } { $i < $n } { incr i } {
    set line [lindex $lines $i]
    if { [string first highest $line] != -1 } {
      incr i
      foreach { at point x y z } [lindex $lines $i] { break }
      incr i
      set line [lindex $lines $i]
      set r [string trim [string range $line 22 30]]
      set aname [string trim [string range $line 31 37]]
      set resname [string trim [string range $line 38 40]]
      set resid [string trim [string range $line 44 end]]
      puts $fileid "$z	$r	$resname $resid"
    }
  }
    close $fileid
    file delete tmpholeinputfiles.pdb
    puts "Calling HOLE Sph_Process"
    catch {exec $holebin\sph_process -dotden 15 -color 2 3 -sos $holeoutput.sph $holeoutput.sos} fid
    set fileid2 [open $holeoutput.sph_process_log.txt a+ 0777]     
    puts $fileid2 $fid
    close $fileid2
 
    puts "Calling HOLE Sos_Triangle" 
    catch {exec $holebin\sos_triangle -s -v < $holeoutput.sos > $holeoutput.vmd} fid2
    set fileid3 [open $holeoutput.sos_triangle_log.txt a+ 0777]
    puts $fileid3 $fid2
    close $fileid3

    puts "Finished"
}


