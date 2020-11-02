# Molecular Dynamics of a Lisozyme proteins
This simulation is made using GROMACS. We follows the instruction on http://www.mdtutorials.com/gmx/lysozyme/

The .mdp are necessary for the simulation.
Those are the shell commands:

    gmx pdb2gmx -f 1aki.pdb -o 1aki_processed.gro -water spce
    gmx editconf -f 1aki_processed.gro -o 1aki_newbox.gro -c -d 1.0 -bt cubic
    gmx solvate -cp 1aki_newbox.gro -cs spc216.gro -o 1aki_solv.gro -p topol.top
    gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
    gmx grompp -f em.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -v -deffnm nvt

Be aware that the mdrun codes can lasts hours!
