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

If you want to visualize the proteins or the system you can edit your .gro file into a .pdb using the editconf command:

    gmx ediconf -f 1AKI_solv_ions.gro -o 1AKI_solv_ions.pdb


Be aware that the mdrun codes can lasts hours!
