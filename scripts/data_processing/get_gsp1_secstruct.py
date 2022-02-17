#!/usr/bin/env python

# this script should be in the scripts/ subdirectory of the project
# last updated 2020-09-29

'''This script loads in each Ran/Gsp1 crystal structure and writes out a 
datafile with secondary structure information
'''

import os, sys, glob
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import seaborn as sns
import pyrosetta
from pyrosetta import *
from pyrosetta import PyMOLMover

from pyrosetta_utils import pose_from_chains, pose_aa_only

from Bio.PDB import *
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

pyrosetta.init()
pymol = PyMOLMover()

# read in data on PDB files, compiled by hand
df = pd.read_csv('../../data/pdb_structures_info.txt', sep='\t')

# get sequences for S. cerevisiae Gsp1 and human Ran (from UNIPROT)
gsp1_record = SeqIO.read('../../data/gsp1.fasta', 'fasta')
ran_record = SeqIO.read('../../data/ran.fasta', 'fasta')

# perform an alignment, so we can map the Ran structure to Gsp1 sequence
alignment = pairwise2.align.globalds(gsp1_record.seq, ran_record.seq,
                                     blosum62, -11, -1) # default param on BLAST
seqs_aligned = {}
seqs_aligned['gsp1'] = alignment[0][0]
seqs_aligned['ran'] = alignment[0][1]

# set directory of PDBs to include  
pdb_dir = '../../data/pdbs_ran/'

# container for df's that will be bound together 
list_of_dfs = []

for pdb in glob.glob(os.path.join(pdb_dir, '*.pdb')): 
        
    pose = pyrosetta.pose_from_pdb(pdb)
    pose_aa = pose_aa_only(pose) # remove Mg and nucleotide

    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
    DSSP.apply(pose_aa)    # populates the pose's Pose.secstruct

    ss = pose_aa.secstruct()
    seq = pose_aa.sequence()
    nres = pose_aa.total_residue()

    ss_df = pd.DataFrame(
        data={'amino_acid': list(seq),
              'secondary_structure': list(ss),
              'pdb_id': [pdb.split('/')[-1].split('.')[0]]*len(seq)}
    )

    alignment = pairwise2.align.localds(
        gsp1_record.seq, seq, blosum62, -11, -1) # default param on BLAST
    aligned_seq = alignment[0][1]
    aa_positions = [i for i, aa in enumerate(aligned_seq) if aa is not '-']
    ss_df['position'] = aa_positions
    
    list_of_dfs.append(ss_df)

# bind together all dataframes and write to file
ss_df = (pd.concat(list_of_dfs)
         [['pdb_id','amino_acid','position','secondary_structure']]
        )
ss_df.to_csv('../../data/gsp1_secondary_structure_long.csv', index=False)

# pivot the dataframe to wide format so that each row is a position
# and each column is a different PDB
ss_df_wide = ss_df.pivot(index='position', columns='pdb_id', values='secondary_structure')
temp_df = pd.DataFrame(data={'amino_acid': list(gsp1_record.seq),
                             'position': list(range(1,220))})
ss_df_wide = temp_df.merge(ss_df_wide, on='position', how='left')

# compute a "consensus" secondary structure assignment by taking the 
# mode across all PDBs. Any NaNs assign as loops (these are just at
# the termini in this case anyway).
ss_df_wide['consensus_ss'] = ss_df_wide.iloc[:,2:].mode(axis='columns')[0]
ss_df_wide.consensus_ss.fillna('L', inplace=True)
ss_df_wide.to_csv('../../data/gsp1_secondary_structure_wide.csv', index=False)
