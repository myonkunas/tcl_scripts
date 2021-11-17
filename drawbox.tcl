

proc box_molecule {molid} {
      # get the min and max values for each of the directions
      # (I'm not sure if this is the best way ... )
      set sel [atomselect top all]

      set coords [lsort -real [$sel get x]]
      set minx [lindex $coords 0]
      set maxx [lindex [lsort -real -decreasing $coords] 0]

      set coords [lsort -real [$sel get y]]
      set miny [lindex $coords 0]
      set maxy [lindex [lsort -real -decreasing $coords] 0]

      set coords [lsort -real [$sel get z]]
      set minz [lindex $coords 0]
      set maxz [lindex [lsort -real -decreasing $coords] 0]
puts "$minx $miny $minz $maxx $maxy $maxz"
      # and draw the lines
      draw materials off
      draw color red
      draw line "$minx $miny $minz" "$maxx $miny $minz"
      draw line "$minx $miny $minz" "$minx $maxy $minz"
      draw line "$minx $miny $minz" "$minx $miny $maxz"

      draw line "$maxx $miny $minz" "$maxx $maxy $minz"
      draw line "$maxx $miny $minz" "$maxx $miny $maxz"

      draw line "$minx $maxy $minz" "$maxx $maxy $minz"
      draw line "$minx $maxy $minz" "$minx $maxy $maxz"

      draw line "$minx $miny $maxz" "$maxx $miny $maxz"
      draw line "$minx $miny $maxz" "$minx $maxy $maxz"

      draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
      draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
      draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
}


proc mybox {x y z} { 
	
	set minx [expr [expr 0 - $x ]/2]
	set maxx [expr $x/2]
	set miny [expr [expr 0 - $y ]/2]
	set maxy [expr $y/2]
	set minz [expr [expr 0 - $z ]/2]
	set maxz [expr $z/2]

	# and draw the lines
      draw materials off
      draw color yellow
      draw line "$minx $miny $minz" "$maxx $miny $minz"
      draw line "$minx $miny $minz" "$minx $maxy $minz"
      draw line "$minx $miny $minz" "$minx $miny $maxz"

      draw line "$maxx $miny $minz" "$maxx $maxy $minz"
      draw line "$maxx $miny $minz" "$maxx $miny $maxz"

      draw line "$minx $maxy $minz" "$maxx $maxy $minz"
      draw line "$minx $maxy $minz" "$minx $maxy $maxz"

      draw line "$minx $miny $maxz" "$maxx $miny $maxz"
      draw line "$minx $miny $maxz" "$minx $maxy $maxz"

      draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
      draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
      draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
}

