mol new d95n_0_working.psf
mol addfile d95n_0_working_min-0.dcd waitfor all
set a [atomselect top protein]
set nf [molinfo top get numframes]
set res [ expr $nf - 1 ]
animate write pdb {d95n_0_working.pdb} beg $res end $res sel $a
exit