#!/bin/bash

# This script runs cartesian_ddg on a set of preminimized pdbs using
# mutfiles to specify the point mutations being made
#
# date: 2020-07-24
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

#ROSETTABIN="/wynton/home/kortemme/cjmathy/Rosetta/main/source/bin"
#ROSETTADB="/wynton/home/kortemme/cjmathy/Rosetta/main/database"
ROSETTABIN="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/source/bin"
ROSETTADB="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/database"

MUTFILE=$(find mutfiles/ -type f -iname "${1}_*")
PDB_ID=${MUTFILE: -8:-4}
PDB="pdbs_for_ddg/"$PDB_ID".pdb"
OUTDIR="ddg_out"
mkdir -p $OUTDIR

$ROSETTABIN/cartesian_ddg.linuxgccrelease \
    -database $ROSETTADB \
    -s $PDB \
    -out:path:all $OUTDIR \
    -ddg:mut_file $MUTFILE \
    -ddg:output_dir $OUTDIR \
    -ddg:iterations 3 \
    -ddg::cartesian \
    -ddg::dump_pdbs false \
    -ddg::bbnbrs 1 \
    -fa_max_dis 9.0 \
    -score:weights ref2015_cart

