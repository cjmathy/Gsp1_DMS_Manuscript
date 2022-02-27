#!/usr/bin/env python

"""This script cleans, renumbers, and prepares mutfiles for the
ddg benchmarking run of 4hb2.pdb

Date: 2021-03-04
Author: Chris Mathy
email: cjmathy@gmail.com
"""

import os
import numpy as np
import pandas as pd


def make_mutfiles(n, aa_from, counter):

    for aa_to in aa1:
        counter+=1
        fname = str(counter)+'_'+aa_from+str(n)+aa_to+'.mut'
        with open('mutfiles/'+fname, 'w') as f:
            f.write('total 1\n')
            f.write('1\n')
            f.write(aa_from+' '+str(n)+' '+aa_to+'\n')
    return



aa1 = ['W','Y','F','M','L','I','V','A','C','G',
       'P','T','S','Q','N','E','D','H','R','K']
aa3 = ['TRP', 'TYR', 'PHE', 'MET', 'LEU',
       'ILE', 'VAL', 'ALA', 'CYS', 'GLY',
       'PRO', 'THR', 'SER', 'GLN', 'ASN',
       'GLU', 'ASP', 'HIS', 'ARG', 'LYS']
aa3to1 = dict(zip(aa3,aa1))

# read in pdb, ignoring REMARKS, waters and the END line
with open('4hb2.pdb', 'r') as f:
    lines = [l for l in f if (
        l.startswith('ATOM') or
        l.startswith('TER') or
        l[17:20]=='GNP' or
        l[18:20]=='MG')]

# now all lines have at least 22 characters, so remove chains B and C
lines = [l for l in lines if l[21]=='A']

# mutfile directory
os.makedirs('mutfiles/', exist_ok=True)

# renumber pdb and make mutfiles
curr_resnum = 0
curr_seqnum = 9
mutfile_counter=0
renumbered_lines = []

numbering_map = []

for l in lines:
    seqnum = l[22:26].strip()
    aa = l[17:20]

    # if at next residue, update resnum counter
    if seqnum != curr_seqnum:
        curr_seqnum = seqnum
        curr_resnum+=1
        numbering_map.append((curr_resnum, seqnum))
        
        if (aa in aa3):
            make_mutfiles(curr_resnum, aa3to1[aa], mutfile_counter)
            mutfile_counter+=20

    # update the line to reflect the renumbering
    renumbered_lines.append(l[:22]+str(curr_resnum).rjust(4)+l[26:])

# write cleaned/renumbered pdb
with open('4hb2_cleaned.pdb','w') as f:
    for l in renumbered_lines:
        f.write(l)

with open('residue_number_map.txt', 'w') as f:
    f.write('resnum,seqnum\n')
    for tup in numbering_map:
        f.write(str(tup[0])+', '+str(tup[1]))
        f.write('\n')
    
