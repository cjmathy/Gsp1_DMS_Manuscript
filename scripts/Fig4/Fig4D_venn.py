import os
os.chdir('/Users/cjmathy/gdrive/gsp1_dms')

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import *

# set up for illustrator to read figures
import matplotlib
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = 'Arial'

# load HRas labels from GOF DMS study and SCA
df = pd.read_csv('data/Ras_data/Ras_sca_gof_table.csv')

tox = df.loc[df.tox].pos.values
sca = df.loc[df.sca].pos.values
gof_all = df.loc[df.gof_all].pos.values
gof_2plus = df.loc[df.gof_2plus].pos.values

# plot GOF venn
fig = plt.figure(figsize=(2, 2), facecolor='w')
ax = fig.add_axes([0, 0, 1, 1])
c = venn2(
    [set(tox), set(gof_2plus)],
    ('Gsp1 toxic/GOF\npositions ({})'.format(len(tox)), 
     'HRas GOF\n positions ({})'.format(len(gof_2plus)))
)

c.get_patch_by_id('10').set_color('#cd2028')    # red
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

plt.title('GOF positions among 156 aligned sites in Gsp1 and HRas')

plt.savefig('figures/Fig4/Fig4D_Venn_gof.pdf', dpi=300, bbox_inches='tight')

# plot SCA venn
fig = plt.figure(figsize=(2, 2), facecolor='w')
ax = fig.add_axes([0, 0, 1, 1])
c = venn2(
    [set(tox), set(sca)],
    ('Gsp1 toxic/GOF\npositions ({})'.format(len(tox)), 
     'HRas sector\n positions (SCA) ({})'.format(len(sca)))
)

c.get_patch_by_id('10').set_color('#cd2028')    # red
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

plt.savefig('figures/Fig4/Fig4E_Venn_sca.pdf', dpi=300, bbox_inches='tight')