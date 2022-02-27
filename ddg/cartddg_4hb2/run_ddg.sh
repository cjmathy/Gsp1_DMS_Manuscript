#!/bin/bash

# This script runs cartesian_ddg on a preminimized Gsp1 pdb using
# mutfiles to specify the point mutations being made
#
# date: 2021-03-04
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

ROSETTABIN="/wynton/home/kortemme/cjmathy/Rosetta/main/source/bin"
ROSETTADB="/wynton/home/kortemme/cjmathy/Rosetta/main/database"
#ROSETTABIN="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/source/bin"
#ROSETTADB="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/database"

MUTFILE=$(find mutfiles/ -type f -iname "${1}_*")
OUTDIR="ddg_out/"
mkdir -p $OUTDIR

$ROSETTABIN/cartesian_ddg.linuxgccrelease \
    -database $ROSETTADB \
    -s 4hb2_cleaned_0007_trunc.pdb \
    -out:path:all $OUTDIR \
    -ddg:mut_file $MUTFILE \
    -ddg:output_dir $OUTDIR \
    -ddg:iterations 5 \
    -ddg::cartesian \
    -ddg::dump_pdbs true \
    -ddg::bbnbrs 1 \
    -fa_max_dis 9.0 \
    -score:weights ref2015_cart

