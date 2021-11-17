
# This script is used to check the geometries after each simulation 

puts " CHECK BY LOOKING AT FIRST" 
proc check_halothane {molnum } { 
	for { set i 0 } { $i < $molnum } { incr i } {
		set resid [expr round(rand()*216)]
		set 

