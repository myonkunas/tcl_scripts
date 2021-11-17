# need psfgen module and topology
package require psfgen
topology top_all27_prot_lipid.inp

# load structures
resetpsf
readpsf  AChR_A7_362POPC_ionized_081304.psf
coordpdb AChR_A7_362POPC_ionized_081304.pdb
readpsf  PoreWater_2.psf
coordpdb PoreWater_2.pdb

# write structure
writepsf AChR_A7_362POPC_ionized_water_081304.psf
writepdb AChR_A7_362POPC_ionized_water_081304.pdb

#quit
