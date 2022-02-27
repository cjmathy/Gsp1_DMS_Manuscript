#!/usr/bin/env python

"""This script reads in the minimized score data
and scores the initial structures to allow for comparison.
Some bulk analysis is run to confirm the structures are 
converged, and the top files are copied to another location
ddg benchmarking run, and a file is written listing these
filenames

Date: 2020-07-24
Author: Chris Mathy
email: cjmathy@gmail.com
"""

import os
import shutil
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pyrosetta
pyrosetta.init()

# read in the minimized score file
sc = (pd.read_csv('pdbs_minimized/score.sc', delim_whitespace=True, skiprows=1)
      .assign(pdb = lambda df: df.description.apply(lambda row: row.split('_')[0]))
      .drop(['SCORE:'], axis=1)
     )

scorefxn = pyrosetta.create_score_function('~/rosetta/Rosetta/main/database/scoring/weights/ref2015_cart.wts')
# score all of the unminized pdbs using PyRosetta
def score_cart(pdb):
    return scorefxn(pyrosetta.pose_from_pdb(pdb))

pdbs = glob.glob('pdbs_renumbered/*.pdb')
init_score = (
    pd.DataFrame(data = {'file':pdbs})
    .assign(pdb = lambda df: df.file.apply(lambda x: re.split(r'/|\.', x)[1]))
    .assign(init_score = lambda df: df.file.apply(lambda pdb: score_cart(pdb)))
    [['pdb','init_score']]
)

# merge in the initial scores
sc = pd.merge(sc, init_score, on='pdb')

# computing change in score for each structure
# shows that every structure was improved in terms of it's energy
sc['delta_score'] = sc.total_score - sc.init_score
print('smallest change in score: {}'.format(max(sc.delta_score)))
plt.hist(sc.delta_score, range=(-1000,0), bins=20)


# overall, we see that our models are well converged, as evidenced by a low
# coefficient of variation across the top5 scoring structures

sc['rank'] = (
    sc.groupby('pdb')
    ['total_score'].rank('first', ascending=True)
    .astype('int64'))

sc_top5 = sc[sc['rank'] <= 5].sort_values(['pdb','rank'])

(sc_top5[['pdb','rank','total_score']]
 .pivot(index='pdb',columns='rank',values='total_score')
 .assign(mean = lambda df: df.max(axis=1),
         sd = lambda df: df.std(axis=1),
         cv = lambda df: df['sd']/df['mean'])
.sort_values('cv', ascending=True)
)

# we can look for some bad models/convergence by plotting the
# mean and standard deviation of the best 5 models

(sc_top5.groupby('pdb')
 .agg(['mean','std'])
 .reset_index()
 ['total_score']
 .pipe(lambda df: sns.scatterplot(data=df, x='mean', y='std'))
)

def print_score_info(sc, pdb):
    total = list(sc[sc.pdb == pdb].total_score)
    delta = list(sc[sc.pdb == pdb].delta_score)
    init = round(list(sc[sc.pdb == pdb].init_score)[0])
    best = min(total)
    mean = round(np.mean(total), 1)
    sd = round(np.std(total), 1)
    delta = round(np.mean(delta), 1)
    print('''{}:
    initial score = {}
    best model score = {}
    best 5 models mean score = {}
    best 5 models mean sd = {}
    best 5 models mean Δscore = {}
    '''.format(pdb, init, best, mean, sd, delta))
    return

print('The only bad scoring minimized model (PDB 1BF4) is nonetheless minimized and well-converged')
print('Note that this model has a DNA molecule, perhaps explaining the relatively high score')
print_score_info(sc, '1BF4')

print('The PDB with the highest standard deviation (PDB 1RHG) had a very bad starting score')
print('Note that this model had long unmodeled loops in the initial xtal structure')
print('The RMSD of the best and 5th best 1RHG models using PyMOL rms_cur is 1.369, and most of the deviation is in the loops')
print_score_info(sc, '1RHG')


#Take the best model from each and write it to a file
best_models = sc[sc['rank'] == 1]['description']

os.makedirs('pdbs_for_ddg', exist_ok=True)

for model in best_models:
    shutil.copy('pdbs_minimized/'+model+'.pdb', 'pdbs_for_ddg/'+model[:4]+'.pdb')

with open('best_minimized_models.txt', 'w') as f:
    for model in list(best_models):
        f.write('pdbs_minimized/'+model+'.pdb\n')
