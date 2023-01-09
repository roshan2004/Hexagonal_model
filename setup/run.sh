#!/usr/bin/bash
gmx grompp -p system.top -f martini_em.mdp -c molecule.gro -o 1-min.tpr -maxwarn 2
gmx mdrun -v -deffnm 1-min -nt 4
gmx grompp -p system.top -f martini_eq.mdp -c 1-min.gro -o 2-eq.tpr -maxwarn 2
gmx mdrun -v -deffnm 2-eq -nt 4
gmx grompp -p system.top -f martini_run.mdp -c 2-eq.gro -o 3-run.tpr -maxwarn 2
gmx mdrun -v -deffnm 3-run -nt 4
