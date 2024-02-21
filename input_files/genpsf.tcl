package require psfgen	 
topology top_all27_prot_lipid_P.inp
pdbalias residue HIS HSE	 
pdbalias atom ILE CD1 CD	 
segment A {pdb d95n_0.pdb}
coordpdb d95n_0.pdb A
guesscoord	 
writepdb d95n_0_working.pdb
writepsf d95n_0_working.psf