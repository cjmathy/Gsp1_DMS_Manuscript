from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import *

# set up for illustrator to read figures
import matplotlib
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = 'Arial'

import os
os.chdir('/Users/cjmathy/gdrive/gsp1_dms')

# get list of toxic positions
toxics = list(
    pd.read_csv('data/toxic_positions.csv')
    .query('position != 220') 
    .position
    )
n_toxics=len(toxics)

# get residues within 4 angstroms of nucleotide from PDB 3m1iA
with open('data/residues_contacting_nucleotide.txt', 'r') as f:
    contacting = [int(line.strip('\n')) for line in f.readlines() if '#' not in line]
n_contact = len(contacting)

# get gtpase regions
gtpase_regions = [
    *range(19, 26+1),    # P-loop 
    *range(34, 47+1),    # Switch I (Canonical Ras superfamily)
    *range(67, 81+1),    # Switch II (Canonical Ras superfamily)
    *range(124, 127+1),  # N/TKxD motif
    *range(152, 154+1)   # SAK motif
]
# n_gr_intersect_contacts = len(set(gtpase_regions).intersection(set(contacting)))
n_gr_union_contacts = len(set(gtpase_regions).union(set(contacting)))

# get core positions
cores = list(
    pd.read_csv('data/burial.csv')
    .query('sasa_group_category=="structure core"')
    .position
)
n_core = len(cores)


### Make Figure 2A
fig = plt.figure(figsize=(2, 2), facecolor='w')
ax = fig.add_axes([0, 0, 1, 1])
c = venn2(
    [set(toxics), set(gtpase_regions).union(set(contacting))],
    ('Toxic/GOF\npositions ({})'.format(n_toxics), 
     'Active site\npositions ({})'.format(n_gr_union_contacts))
)

c.get_patch_by_id('10').set_color('#CD2029')    # red
c.get_patch_by_id('10').set_edgecolor('black')
c.get_patch_by_id('10').set_alpha(0.6)
c.get_patch_by_id('01').set_color('#FFFFFF')    # white
c.get_patch_by_id('01').set_edgecolor('black')   
c.get_patch_by_id('01').set_alpha(1)
c.get_patch_by_id('11').set_color('#4169E1')    # royal blue
c.get_patch_by_id('11').set_edgecolor('black')
c.get_patch_by_id('11').set_alpha(0.6)

lblA = c.get_label_by_id("A")
lblA.set_position((-0.7, 0.05))
lblB = c.get_label_by_id("B")
lblB.set_position((0.7, 0.05))


plt.savefig('figures/Fig2/Fig2B_venn.pdf', dpi=300, bbox_inches='tight')




# What are the positions in the GTPase region that contact the nucleotide 
# but are not T/GOF

(pd.read_csv('data/Gsp1_fitness_scores.csv')
 .query('position in [41,124,154]')
 .groupby(['position', 'bin'])
 .size()
)

# # Could include a 3-category venn diagram:
# fig = plt.figure(facecolor='w')
# c = venn3(
#     [set(toxics), set(gtpase_regions), set(contacting)],
#     ('Toxic/GOFs',
#      'Positions in the GTPase regions',
#      'Residues contacting the nucleotide'))

# def genbin(l, n, bs=''):
#     if len(bs) == n:
#         l.append(bs)
#     else:
#         genbin(l, n, bs + '0')
#         genbin(l, n, bs + '1')
#     return l

# for i in genbin([], 3):
#     if i != '000':
#         c.get_patch_by_id(i).set_edgecolor('black')
# plt.savefig('figures/Fig2/Fig2A_venn3.pdf', dpi=300, bbox_inches='tight')


## Make Figure 2B
fig = plt.figure(figsize=(2, 2), facecolor='w')
ax = fig.add_axes([0, 0, 1, 1])
c = venn2(
    [set(toxics), set(cores)],
    ('Toxic/GOF\npositions ({})'.format(n_toxics), 
     'Core\npositions ({})'.format(n_core))
)

c.get_patch_by_id('10').set_color('#FFFFFF')    # white
c.get_patch_by_id('10').set_edgecolor('black')
c.get_patch_by_id('10').set_alpha(1)
c.get_patch_by_id('01').set_color('#FFA54F')    # tan1
c.get_patch_by_id('01').set_edgecolor('black')   
c.get_patch_by_id('01').set_alpha(0.6)
c.get_patch_by_id('11').set_color('#CD2029')    # red
c.get_patch_by_id('11').set_edgecolor('black')
c.get_patch_by_id('11').set_alpha(0.6)

lblA = c.get_label_by_id("A")
lblA.set_position((-0.7, 0.05))
lblB = c.get_label_by_id("B")
lblB.set_position((0.7, 0.05))

plt.savefig('figures/Fig2/Fig2C_venn.pdf', dpi=300, bbox_inches='tight')








