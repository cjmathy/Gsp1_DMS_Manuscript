#!/bin/bash

# This script scores the input PDBs for cartesian_ddg
# (the original, cleaned PDB and the minimized model)
#
# date: 2021-04-20
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

ROSETTABIN="/Users/cjmathy/Rosetta/main/source/bin"
ROSETTADB="/Users/cjmathy/Rosetta/main/database"

$ROSETTABIN/score_jd2.macosclangrelease \
    -database $ROSETTADB \
    -s 3gj0_cleaned.pdb \
    -fa_max_dis 9.0 \
    -score:weights ref2015_cart

$ROSETTABIN/score_jd2.macosclangrelease \
    -database $ROSETTADB \
    -s preminimization_out/3gj0_cleaned_0009.pdb \
    -fa_max_dis 9.0 \
    -score:weights ref2015_cart
