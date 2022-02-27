#!/usr/bin/env python

"""script to merge various ddg benchmark datasets

Date: 2020-07-17
Author: Chris Mathy
email: cjmathy@gmail.com
"""

import pandas as pd
import numpy as np
import re

# function to convert 'A G 44 S' format (chain aa_from position aa_to) to G44S
# for kellogg/curated_protherm datasets
def convert_mutation_format(s):
    vec = re.split(r' |_', s)
    if len(vec) == 4:
        return ''.join(vec[1:])
    if len(vec) == 8:
        return ''.join(vec[1:4])+'/'+''.join(vec[5:8])

# full datasets
# pd.read_csv('kellogg.csv', skiprows=21)                   # 1210 rows
# pd.read_csv('curatedprotherm.csv', skiprows=21)           # 2971 rows
# pd.read_csv('ref2015.Lizddg.txt', delim_whitespace=True)  # 1181 rows
# pd.read_csv('balanced_benchmark.csv')                     # 768 rows    

# Kellogg and Protherm can be merged on Mutations and PDBFileID, they have 801 in common
# pd.merge(protherm, kellogg, on=['PDBFileID','Mutations'])

# read in liz kellogg's dataset and ignore multiple mutations
kellogg = (pd.read_csv('kellogg.csv', skiprows=21)
           .loc[lambda df: df.DSSPSimpleTypes.apply(len) == 1]  # ignore records w/ >1 mutation
           .assign(chain = lambda df: df.Mutations.apply(
               lambda row: row[0]))
           .assign(mutation_short = lambda df: df.Mutations.apply(
               lambda row: convert_mutation_format(row)))
           .rename(columns={'#RecordID':'kellogg_ID','PDBFileID':'pdb', 'Mutations':'mutation',
                            'DDG':'ddg_expt', 'IndividualDDGs':'ddg_expt_individ'})
           [['pdb','mutation','chain','mutation_short','kellogg_ID','ddg_expt','ddg_expt_individ']]
)

# read in the protherm dataset and ignore multiple mutations   
protherm = (pd.read_csv('curatedprotherm.csv', skiprows=21)
            .loc[lambda df: df.DSSPSimpleTypes.apply(len) == 1]  # ignore records w/ >1 mutation
            .assign(chain = lambda df: df.Mutations.apply(
                lambda row: row[0]))
            .assign(mutation_short = lambda df: df.Mutations.apply(
                lambda row: convert_mutation_format(row)))
            .rename(columns={'RecordID':'protherm_ID','PDBFileID':'pdb', 'Mutations':'mutation',
                            'DDG':'ddg_expt', 'IndividualDDGs':'ddg_expt_individ'})
            [['pdb','mutation','chain','mutation_short','protherm_ID','ddg_expt','ddg_expt_individ']]
)

# read in hahnbeom park's dataset
# NOTE we don't process it here, as it uses the renumbered pdb numbering
hp = (pd.read_csv('ref2015.Lizddg.txt', delim_whitespace=True)
      .rename(columns={"#": 'pdb', 'Expt.':'ddg_expt', 'Calc.':'ddg_pred'})
      .assign(mutation_short = lambda df: df.pdb.str.split('_',expand=True)[1])
      .assign(pdb = lambda df: df.pdb.str.split('_',expand=True)[0].str.upper())
      [['pdb', 'mutation_short', 'ddg_expt', 'ddg_pred']]
     )

# read in brandon frenz's dataset
frenz = (pd.read_csv('balanced_benchmark.csv')
         .assign(mutation_short = lambda df: df[['Wild Type', 'Residue Number', 'Mutation']].apply(
             lambda row: ''.join(row.values), axis=1))
         .rename(columns={'Protherm ID':'protherm_ID','PDB ID':'pdb','Experimental DDG':'ddg_expt'})
         [['protherm_ID','pdb', 'mutation_short', 'ddg_expt']]
)

df = pd.merge(protherm, kellogg,
              on=['pdb','chain','mutation_short','mutation'],
              how='outer', suffixes=('_protherm','_kellogg'))

(pd.merge(df, frenz, on = ['protherm_ID','pdb','mutation_short'], how='outer')
 .rename(columns={'ddg_expt':'ddg_expt_frenz'})
 .rename_axis('record_ID')
 .to_csv('benchmark_sets_merged_new.csv')
)
