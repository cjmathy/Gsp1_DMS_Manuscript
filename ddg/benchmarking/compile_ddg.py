import os
import glob
import re
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# sns.set(style="whitegrid", palette="pastel", color_codes=True, font='Helvetica')

# Set font sizes
SMALL_SIZE = 12 #6
MEDIUM_SIZE = 13 #7
BIG_SIZE = 14 #8

plt.rcParams.update({
    'legend.loc': 'center left',
#     'font.family': 'Helvetica',
    'font.size': BIG_SIZE,         # controls default text sizes
    'axes.titlesize': SMALL_SIZE,  # fontsize of the axes title
    'axes.labelsize': MEDIUM_SIZE, # fontsize of the x and y labels
    'xtick.labelsize': SMALL_SIZE, # fontsize of the tick labels
    'ytick.labelsize': SMALL_SIZE, # fontsize of the tick labels
    'legend.fontsize': SMALL_SIZE, # legend fontsize
    'figure.titlesize': BIG_SIZE,   # fontsize of the figure title
    'figure.figsize': (6,6),
    'figure.dpi': 120
})

def read_ddg(ddgfile):
    with open(ddgfile,'r') as f:
        lines = f.readlines()
        rows = []
        for row, line in enumerate(lines):            
            l = line.split()
            d = {}
            d['round'] = l[1][5]
            d['type'] = l[2][:-1]
            d['total_score'] = l[3]
            for i in np.arange(4, len(l), 2):
                d[l[i][:-1]]=l[i+1]
            rows.append(pd.DataFrame(data=d,index=[row]))
        df = pd.concat(rows)
        df['file'] = ddgfile
    return df

dfs = []
for ddgfile in glob.glob('ddg_out_second_run/*.ddg'):
    dfs.append(read_ddg(ddgfile))

df = (pd.concat(dfs)
      .reset_index()
      .drop('index', 1)
      .assign(total_score=lambda df: df.total_score.astype(float),
              record_ID=lambda df: df.file.apply(lambda row: row.split('/')[1].split('_')[0]),
              calc=lambda df: df.type.apply(lambda row: row.split('_')[0])
             )
      .assign(record_ID=lambda df: df.record_ID.astype('int64'))
      [['record_ID','calc','total_score']]
     )


df_ddg = (df
    .pivot_table(index='record_ID', columns='calc', values='total_score')
    .assign(ddg_calc = lambda df: df.MUT-df.WT)
    .reset_index()
    [['record_ID','ddg_calc']]
)

dataset = pd.read_csv('datasets/benchmark_sets_merged.csv')

merged_df = pd.merge(dataset, df_ddg, on='record_ID')
merged_df.head()
merged_df.to_csv('ddg_benchmark_dataset.csv', index=False)

protherm_corr = 'pearson corr: = {}'.format(round(merged_df['ddg_expt_protherm'].corr(merged_df['ddg_calc']), 3))
kellogg_corr = 'pearson corr: = {}'.format(round(merged_df['ddg_expt_kellogg'].corr(merged_df['ddg_calc']), 3))
frenz_corr = 'pearson corr: = {}'.format(round(merged_df['ddg_expt_frenz'].corr(merged_df['ddg_calc']), 3))

print('Protherm {}'.format(protherm_corr))
print('Protherm {}'.format(kellogg_corr))
print('Protherm {}'.format(frenz_corr))

df_to_plot = (
    merged_df
    [['record_ID','pdb','mutation', 'ddg_calc', 'ddg_expt_protherm', 'ddg_expt_kellogg', 'ddg_expt_frenz']]
    .melt(id_vars=('record_ID','pdb','mutation', 'ddg_calc'),
                value_vars=('ddg_expt_protherm', 'ddg_expt_kellogg', 'ddg_expt_frenz'),
                var_name = 'dataset', value_name='ddg_expt')
    .assign(dataset=lambda df: df.dataset.apply(lambda row: row.split('_')[2]))
)

g = sns.relplot(data=df_to_plot, x='ddg_calc', y='ddg_expt', col='dataset')
print('Protherm dataset {}'.format(protherm_corr))
print('Kellogg dataset {}'.format(kellogg_corr))
print('Protherm dataset {}'.format(frenz_corr))

corrs = [protherm_corr, kellogg_corr, frenz_corr]

for i, ax in enumerate(g.axes[0]):
    ax.axhline(0, ls='--',color='red')
    ax.axvline(0, ls='--',color='red')
    ax.text(30,15,corrs[i])
    
    
plt.show()
g.savefig('dataset_corr.png')

# want to compute a scaling factor based on the kellogg data

vals_for_regression = (df_to_plot
 [df_to_plot['dataset'] == 'kellogg']
 .dropna()
 [['ddg_calc','ddg_expt']]
)

x = vals_for_regression.ddg_calc.values
y = vals_for_regression.ddg_expt.values

popt, pcov = sp.optimize.curve_fit(lambda x, slope: x*slope, x, y)

scaling_factor = np.round(popt[0],3)

# plot a version with just kellogg
g = sns.relplot(data=df_to_plot[df_to_plot['dataset'] == 'kellogg'],
                x='ddg_calc', y='ddg_expt', col='dataset')

corrs = [kellogg_corr]

ax = g.axes[0][0]
ax.axhline(0, ls='--',color='red')
ax.axvline(0, ls='--',color='red')
ax.text(10,-4,round(merged_df['ddg_expt_kellogg'].corr(merged_df['ddg_calc']), 3))

xline = np.linspace(-10, 25, 36)
yline = scaling_factor*xline

plt.plot(xline, yline, color='green')
    
plt.show()

print('scaling factor = {}'.format(scaling_factor))

# now we want to look at how the Park et al 2016 data compares on these slices

hp = (pd.read_csv('datasets/ref2015.Lizddg.txt', delim_whitespace=True)
      .rename(columns={"#": 'pdb', 'Expt.':'ddg_expt', 'Calc.':'ddg_pred'})
      .assign(mutation_short = lambda df: df.pdb.str.split('_',expand=True)[1])
      .assign(pdb = lambda df: df.pdb.str.split('_',expand=True)[0].str.upper())
      [['pdb', 'mutation_short', 'ddg_expt', 'ddg_pred']]
     )
hp

# these mutations are in pdb numbering, so read in the mutfiles

record_IDs = []
mutations = []
pdbs = []

for mutf in glob.glob('mutfiles/*.mut'):
    record_IDs.append(mutf.split('/')[1].split('_')[0])
    pdbs.append(mutf.split('/')[1].split('_')[1][:4])

    with open(mutf, 'r') as f:
        lines = [l for l in f]
        mutations.append(lines[2].replace(' ','')[:-1])

        
        
df_resnum = pd.DataFrame(data={'record_ID': record_IDs,
                               'mutation_short':mutations,
                               'pdb':pdbs})

hp = pd.merge(hp, df_resnum, on=['pdb','mutation_short'])


# get the record ID's used for the kellogg dataset
kellogg_subset_computed = (
    merged_df
    [['record_ID', 'pdb', 'mutation_short', 'ddg_expt_kellogg', 'ddg_calc']]
    .dropna()
    .assign(record_ID=lambda df: df.record_ID.astype(str))
)
print(round(kellogg_subset_computed['ddg_expt_kellogg'].corr(kellogg_subset_computed['ddg_calc']), 3))

record_ID_used = set(kellogg_subset_computed.record_ID.values)

hp_subset = hp[hp['record_ID'].isin(record_ID_used)]
hp_subset_corr = round(hp_subset['ddg_expt'].corr(hp_subset['ddg_pred']),3)
g = sns.relplot(data=hp_subset, x='ddg_pred', y='ddg_expt')
ax = g.axes[0][0]
ax.axhline(0, ls='--',color='red')
ax.axvline(0, ls='--',color='red')
ax.text(10,-4,hp_subset_corr)

g.savefig('park_ddg_subset.png')