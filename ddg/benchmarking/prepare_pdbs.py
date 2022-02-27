#!/usr/bin/env python

"""This script cleans, renumbers, and prepares mutfiles for the
ddg benchmarking run.

Date: 2020-07-23
Author: Chris Mathy
email: cjmathy@gmail.com
"""

import os
import numpy as np
import pandas as pd


def clean_pdb(pdbfile, outdir):
    fname = pdbfile.split('/')[-1]

    # read in PDB, excluding HETATOM lines
    with open(pdbfile, 'r') as f:
        lines = [l for l in f if (l.startswith('ATOM') or l.startswith('TER'))]

    # write out this cleaned version
    with open(os.path.join(outdir,fname), 'w') as g:
        for line in lines:
            g.write(line)

    return lines


def renumber_pdb(pdbid, pdblines, df, renumdir, mutdir):

    df = df.assign(
        chain=lambda df:df.mutation.apply(lambda mut: mut.split(' ')[0]),
        aa_from=lambda df:df.mutation.apply(lambda mut: mut.split(' ')[1]),
        seqnum=lambda df:df.mutation.apply(lambda mut: mut.split(' ')[2]),
        aa_to=lambda df:df.mutation.apply(lambda mut: mut.split(' ')[3]))

    chains_with_mutations = set(df.chain.values)
    seqnums_with_mutations = set(df.seqnum.values)

    curr_resnum = 1
    curr_seqnum = pdblines[0][22:26].strip()
    renumbered_lines = []

    for line in pdblines:
        chain = line[21]

        seqnum = line[22:26].strip()

        if curr_seqnum != seqnum:
            curr_seqnum = seqnum
            curr_resnum+=1

        if (not(line.startswith('ATOM') or line.startswith('TER'))):
            print('\n--THIS LINE OF {} IS NEITHER ATOM/TER--'.format(pdbid))
            print(line)

        else:
            # update the line to reflect the renumbering
            renumbered_lines.append(line[:22]+str(curr_resnum).rjust(4)+line[26:])

            # if this residue needs to be mutated, make the mutfile(s)
            if seqnum in seqnums_with_mutations:
                make_mutfile(df[df.seqnum==seqnum], curr_resnum, mutdir)
                seqnums_with_mutations.remove(seqnum)

    # write the renumbered pdb to a file
    with open(os.path.join(renumdir,pdbid+'.pdb'), 'w') as f:
        for line in renumbered_lines:
            f.write(line)

    return


def make_mutfile(df, res, mutdir):

    for _, row in df.iterrows():
        pdb = row['pdb']
        rid = str(row['record_ID'])
        mut = row['mutation']
        mut = ' '.join((mut[2], str(res), mut[-1]))

        with open(os.path.join(mutdir,rid+'_'+pdb+'.mut'), 'w') as f:
            f.write('total 1\n')
            f.write('1\n')
            f.write(mut+'\n')

    return


cleandir = 'pdbs_cleaned'
renumdir = 'pdbs_renumbered'
mutdir = 'mutfiles'

for newdir in (cleandir, renumdir, mutdir):
    os.makedirs(newdir, exist_ok=True)

dset = pd.read_csv('datasets/benchmark_sets_merged.csv')

for pdbid in np.unique(dset.pdb.values):

    # get info for that pdb
    pdbfile = os.path.join('ddg/input/pdbs', pdbid+'.pdb')

    mutations_df = dset[dset.pdb == pdbid][['record_ID','pdb','mutation']]

    # clean that pdb
    pdblines = clean_pdb(pdbfile, cleandir)

    # renumber and simultaneously write out mutfiles using the unique record_ID
    renumber_pdb(pdbid, pdblines, mutations_df, renumdir, mutdir)


