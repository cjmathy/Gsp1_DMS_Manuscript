#!/bin/bash

# This script runs preminimization for the cartesian_ddg protocol on a
# cleaned Gsp1 structure
#
# date: 2020-07-24
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

ROSETTABIN="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/source/bin"
ROSETTADB="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/database"
PDB="3m1i_cleaned.pdb"
OUTDIR="preminimization_out/"

mkdir -p $OUTDIR

$ROSETTABIN/relax.default.linuxgccrelease \
    -s $PDB \
    -out:path:all $OUTDIR \
    -database $ROSETTADB \
    -in:file:movemap gsp1.movemap \
    -use_input_sc \
    -constrain_relax_to_start_coords \
    -ignore_unrecognized_res \
    -nstruct 20 \
    -relax:cartesian \
    -relax:coord_constrain_sidechains  \
    -relax:min_type lbfgs_armijo_nonmonotone \
    -score:weights ref2015_cart \
    -relax:script cart2.script > $OUTDIR"minimization.log"


