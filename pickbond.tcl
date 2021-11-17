
# Pickbond
# Author: Justin Gullingsrud
# Date: 11 January, 2001
# Compatible VMD versions: Post-1.6; i.e, anything built after the above 
#    date ;-)  You can rewrite the delbond and addbond procs below to reflect 
#    the current state of VMD's bond changing API; it's undergone a bit of flux.


# This is a script to allow users to add or delete bonds by picking atoms
# with the mouse.  Type "pickbond add" at the text console, and now picking
# pairs of atoms will result a new bond being created between them.
# "pickbond del" will cause bonds between picked pairs of atoms to be 
# deleted.  "pickbond stop" ceases all pickbond behavior.


proc delbond { molid atom1 atom2 } {
  if { $atom1 > $atom2 } {
    set tmp $atom1
    set atom1 $atom2
    set atom2 $tmp
  }
  set sel [atomselect $molid "index $atom1 $atom2"]
  lassign [$sel getbonds] bond1 bond2
  
  set id [lsearch -exact $bond1 $atom2]
  if { $id != -1 } {
    set bond1 [lreplace $bond1 $id $id ]
  }
  set id [lsearch -exact $bond2 $atom1]
  if { $id != -1 } {
    set bond2 [lreplace $bond2 $id $id ]
  }
  $sel setbonds [list $bond1 $bond2] 
  $sel delete
}

proc addbond { molid atom1 atom2 } {
  if {$atom1 == $atom2 } {
    #error "cannot bond an atom to itself"
    return
  }
  if { $atom1 > $atom2 } {
    set tmp $atom1
    set atom1 $atom2
    set atom2 $tmp
  }
  set sel [atomselect $molid "index $atom1 $atom2"]
  lassign [$sel getbonds] bond1 bond2

  set id [lsearch -exact $bond1 $atom2]
  if { $id == -1 } {
    lappend bond1 $atom2
  }
  set id [lsearch -exact $bond2 $atom1]
  if { $id == -1 } {
    lappend bond2 $atom1
  }
  $sel setbonds [list $bond1 $bond2] 
  $sel delete
}

# set this proc to trace on vmd_pick_atom

set pickbond_molid -1
set pickbond_atom1 -1
set pickbond_proc addbond

proc pickbond_trace { name1 name2 op } {
  global pickbond_proc pickbond_molid pickbond_atom1 vmd_pick_mol vmd_pick_atom
 
  if { $pickbond_molid == -1 } {
    # this is the first atom
    set pickbond_molid $vmd_pick_mol
    set pickbond_atom1 $vmd_pick_atom
  } else {
    # we have two atoms: call the proc
    catch {$pickbond_proc $pickbond_molid $pickbond_atom1 $vmd_pick_atom}
    set pickbond_molid -1
    set pickbond_atom1 -1
  }
}

      
#
# pickbond [add | del | stop ] 
#
proc pickbond { op } {
  global pickbond_molid pickbond_atom1 pickbond_proc

  set pickbond_molid -1
  set pickbond_atom1 -1

  uplevel {trace vdelete vmd_pick_atom w pickbond_trace}

  if { $op == "add" } {
    set pickbond_proc addbond
    uplevel {trace variable vmd_pick_atom w pickbond_trace}
  } elseif { $op == "del" } {
    set pickbond_proc delbond
    uplevel {trace variable vmd_pick_atom w pickbond_trace}
  } elseif { $op == "stop" } {
    puts "No more pickbond"
  } else {
    error "Unknown option: choose one of add, del, or stop"
  }
}


