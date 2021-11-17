#This script fits all num_pdb to the first pdb
#according to the {fit_text},
#and calculates an average structure of all frames, 
#


#Replace the fit_text to what you want to select
set fit_text "backbone and resid 4 to 19"

#change number of pdb files read in
set num_pdb 30

for {set i 0} {$i < $num_pdb} {incr i} {
  set j [expr $i+1]
  mol load pdb "low_$j.pdb"
  mol delrep 0 top
  mol representation lines
  mol color name
mol selection {backbone and name C or name N or name CA}
mol material Opaque
mol addrep top
  set selall($i) [atomselect $i all]
  set fit($i) [atomselect $i $fit_text]
}



#align all structures with respect to 1st structure to remove translation and rotation
        for {set i 0} {$i < $num_pdb} {incr i} {             
                set trans_mat [measure fit $fit($i) $fit(0)]
                $selall($i) move $trans_mat
        }




