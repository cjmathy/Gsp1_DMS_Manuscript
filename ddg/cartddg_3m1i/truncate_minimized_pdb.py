#!/usr/bin/env python

"""This script cleans, renumbers, and prepares mutfiles for the
ddg benchmarking run of 3m1i.pdb

Date: 2020-07-24
Author: Chris Mathy
email: cjmathy@gmail.com
"""

import os



with open('preminimization_out/3m1i_cleaned_0012.pdb', 'r') as f:
    lines = [l for l in f if (
        l.startswith('ATOM') or
        l.startswith('TER') or
        l[17:20]=='GTP' or
        l[17:19]=='MG')]

# write 

lines_new = []

for l in lines:
    if not (int(l[23:26]) > 174) and (int(l[23:26]) < 200):
        lines_new.append(l)
    elif (int(l[23:26]) == 200):
        lines_new.append(l[:23]+' 175')

for l in lines_new:
    if int(l[23:26]) >= 200:

with open('3m1i_cleaned_0012_trunc.pdb','w') as g:
    for l in renumbered_lines:
        g.write(l)

