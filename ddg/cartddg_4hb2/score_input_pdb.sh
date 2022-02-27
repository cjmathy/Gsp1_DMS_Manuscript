#!/bin/bash

# This script scores the input PDB for cartesian_ddg
# mutfiles to specify the point mutations being made
#
# date: 2021-03-04
# author: Chris Mathy
# email: chris.mathy@ucsf.edu
# email: cjmathy@gmail.com

ROSETTABIN="/Users/cjmathy/Rosetta/main/source/bin"
ROSETTADB="/Users/cjmathy/Rosetta/main/database"

$ROSETTABIN/score_jd2.macosclangrelease \
    -database $ROSETTADB \
    -s 4hb2_cleaned_0007_trunc.pdb \
    -fa_max_dis 9.0 \
    -score:weights ref2015_cart

