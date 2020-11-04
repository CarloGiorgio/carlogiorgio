# Molecular Dynamics of a Lisozyme proteins
This simulation is made using GROMACS. We follows the instruction on http://www.mdtutorials.com/gmx/lysozyme/

The .mdp are necessary for the simulation.
Those are the shell commands:

    gmx pdb2gmx -f 1aki.pdb -o 1aki_processed.gro -water spce
This command create a file that Gromacs can hand

    gmx editconf -f 1aki_processed.gro -o 1aki_newbox.gro -c -d 1.0 -bt cubic
This insert the molecule inside the correct geometry. We can choose among different ones.

    gmx solvate -cp 1aki_newbox.gro -cs spc216.gro -o 1aki_solv.gro -p topol.top
It adds solvate spc216 is the type of model of water used.

    gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

This command neutralize the charge using Na and CL

    gmx grompp -f em.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em
    
This command create the microcanonical equilibrium for the system
In this case, the output is a .tpr file and not a .gro!
If you whant to know the output just use the following command:
    
    
More inforrmation about output and so on can be found on http://manual.gromacs.org/documentation/current/onlinehelp/gmx-energy.html.
    
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -v -deffnm nvt

If you want to visualize the proteins or the system you can edit your .gro file into a .pdb using the editconf command:

    gmx ediconf -f 1AKI_solv_ions.gro -o 1AKI_solv_ions.pdb


Be aware that the mdrun codes can lasts hours!
