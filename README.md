# SimulationStats
A Java suite of post-analysis tools for molecular modeling simulations

The current state of the software includes a "core" framework that handles basic information  
of atoms and simulation boxes. Currently, the only file format supported is Tinker formatted .arc and .xyz files.  
It has one calculation implemented thus far, which is a RDF calculator.

New in recent update:
Contains a .pdb to lammps .dat file builder. Can take in a pdb formatted file constructed from Packmol and convert it into a lammps input file. The code asks the user to provide "blueprint" files, which act as a way to assign bonds and angles to different molecules. After reading in the system file and any pdb "blueprint" files, then it will ask the user for the mass and charge of atoms contained in the system.
