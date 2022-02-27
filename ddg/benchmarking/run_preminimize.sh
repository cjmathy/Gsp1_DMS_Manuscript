#!/bin/bash

# This script runs preminimization for the cartesian_ddg protocol on a set of
# cleaned PDB files.
#
# date: 2020-07-22
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

#ROSETTABIN="/wynton/home/kortemme/cjmathy/Rosetta/main/source/bin"
#ROSETTADB="/wynton/home/kortemme/cjmathy/Rosetta/main/database"
ROSETTABIN="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/source/bin"
ROSETTADB="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/database"

PDB=$(sed "$1q;d" pdbs_to_minimize.txt)
PDBDIR="pdbs_renumbered/"
OUTDIR="pdbs_minimized/"

mkdir -p $OUTDIR

$ROSETTABIN/relax.default.linuxgccrelease \
    -s $PDBDIR$PDB \
    -out:path:all $OUTDIR \
    -database $ROSETTADB \
    -use_input_sc \
    -constrain_relax_to_start_coords \
    -ignore_unrecognized_res \
    -nstruct 20 \
    -relax:cartesian \
    -relax:coord_constrain_sidechains  \
    -relax:min_type lbfgs_armijo_nonmonotone \
    -score:weights ref2015_cart \
    -relax:script cart2.script

