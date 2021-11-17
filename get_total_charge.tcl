proc get_total_charge {{molid }} {
	eval "vecadd [[atomselect $molid all] get charge]"
}
