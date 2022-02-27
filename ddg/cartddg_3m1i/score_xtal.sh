#!/bin/bash

# This script scores the input cleaned Gsp1 structure
# date: 2020-07-24
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

ROSETTABIN="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/source/bin"
ROSETTADB="/kortemmelab/home/cjmathy/rosetta/Rosetta/main/database"
PDB="3m1i_cleaned.pdb"
OUTDIR="preminimization_out/"

mkdir -p $OUTDIR

$ROSETTABIN/score_jd2.default.linuxgccrelease \
    -s $PDB \
    -database $ROSETTADB \
    -ignore_unrecognized_res \
    -score:weights ref2015_cart \
    -out:file:scorefile score_3m1i_cleaned.sc

