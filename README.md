<h1 align="center">MC_side_chain_conf</h1>
<h2 align="center">

## Description
  
**Monte Carlo simulation for protein's rotamers**
  
  The algorithm finds the optimal conformation of rotamers using Monte Carlo simulation. At each step, a side-bond around which the rotation will occur is randomly chosen. If, after changing the protein geometry, the total energy is less than the best energy, the best energy variable is overwritten.
  
  The cycle continues until a certain number of steps are passed or all available change attempts are made that do not result in a better conformation.
The output includes a log file and a file with bond numbers and rotation angles in degrees, as well as the best energy.
  
  Initial data are presented as a config file and a protein in pdb format. 

## How to run
  - **Protein choice** 
  
    Place the protein's pdb file in the pdb_files directory
  - **Parametres choice**
    
    Enter the parameters in the config file
  - **Running code**
  
    Run run_MC.py

## Project puckages

```
abc
random
configparser
pathlib
collections
math
numpy
```

## Future scope

- Add energy function
- Add aminoacids: PRO, GLN, HIS
- Add function to analyze pdb file and refill missing atoms
