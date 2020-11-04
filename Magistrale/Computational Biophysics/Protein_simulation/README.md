# Molecular Dynamics of a Lisozyme proteins
This simulation is made using GROMACS. We follow the instruction on http://www.mdtutorials.com/gmx/lysozyme/

The .mdp are necessary for the simulation.
Those are the shell commands:

    gmx pdb2gmx -f 1aki.pdb -o 1aki_processed.gro -water spce
This command creates files that Gromacs can handle: a topology file, a positions restrained file and a post-processed structure file.
We choose SPC/E model for water molecules.

    gmx editconf -f 1aki_processed.gro -o 1aki_newbox.gro -c -d 1.0 -bt cubic
This inserts the molecule inside the correct geometry. We can choose among different ones.

    gmx solvate -cp 1aki_newbox.gro -cs spc216.gro -o 1aki_solv.gro -p topol.top
It adds solvate to the topology; spc216.gro is a good solvent configuration for our water model.

    gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

This command neutralizes the charge using Na+ and Cl- ions (in this case 8 Cl- are needed).

    gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em
    
This command creates the microcanonical equilibrium for the system by minimizing energy.
In this case, the output is a .tpr file and not a .gro!
If you want to know the output just use the following command:

    gmx energy -f em.edr -o potential.xvg
    
More information about output and so on can be found on http://manual.gromacs.org/documentation/current/onlinehelp/gmx-energy.html.
    
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -v -deffnm nvt

This command creates an NVT ensemble by coupling the system to a thermostat for 100 ps.
A similar command equilibrates also the system pressure, reaching an NPT ensemble.
Such MD runs require several minutes of execution.

If you want to visualize the proteins or the system you can edit your .gro file into a .pdb using the editconf command:

    gmx ediconf -f 1AKI_solv_ions.gro -o 1AKI_solv_ions.pdb

Now we are ready for the MD simulation:

    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
    gmx mdrun -v -deffnm md

This code makes the protein evolve for 1 ns. Be aware that it can lasts hours on your PC! Use "-v" tag to see the run progress.

    gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
    gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns

The first line corrects any positional error due to periodic boundary conditions. Next line calculates the root mean squared of the backbones at each frame of the simulation. Other useful commands for statistical analysis of the trajectory can be found at http://manual.gromacs.org/documentation/5.1/user-guide/cmdline.html#commands-by-topic.
