#!/usr/bin/env python
# coding: utf-8

import glob
import re
import numpy as np
import pandas as pd


# Function to parse the ddg output file created by cartesian_ddg
def read_ddg(ddgfile):
    with open(ddgfile,'r') as f:
        lines = f.readlines()
        rows = []
        for row, line in enumerate(lines):            
            l = line.split()
            d = {}
            d['round'] = l[1][5]
            d['type'] = l[2][:-1]
            d['total_score'] = float(l[3])
            for i in np.arange(4, len(l), 2):
                
                colname = l[i][:-1] 
                value = l[i+1]

                # set correct dtypes                
                if colname not in ['round', 'type', 'file']:
                    value = float(value)

                d[colname]=value
            rows.append(pd.DataFrame(data=d,index=[row]))
        df = pd.concat(rows)
        df['file'] = ddgfile
    return df

# convert rosetta numbering (mutation file name uses this) to
# sequence numbering
res2seq = (pd.read_csv('3gj0_residue_number_map.txt')
           .set_index('resnum').to_dict()['seqnum']
          )

# read in all results (scores for each MUT and paired WT), concatenate in a dataframe
dfs = []
for ddgfile in glob.glob('ddgs/*.ddg'):
    dfs.append(read_ddg(ddgfile))

df = (pd.concat(dfs)
      .reset_index()
      .drop('index', 1)
      .assign(total_score=lambda df: df.total_score.astype(float),
              mut=lambda df: df.file.apply(lambda row: row.split('_')[1][:-4]),
              mutation=lambda df: df.mut.apply(lambda row: row[0]+str(res2seq[float(row[1:-1])])+row[-1]),
              aa_from=lambda df: df.mutation.apply(lambda row: row[0]),
              position=lambda df: df.mutation.apply(lambda row: int(row[1:-1])),
              aa_to=lambda df: df.mutation.apply(lambda row: row[-1]),
              run_type=lambda df: df.type.apply(lambda row: row.split('_')[0])
             )
      .groupby(['mutation','aa_from','position','aa_to','run_type'])
      .agg("mean")
      .reset_index()
      .sort_values(['position','aa_to'])
     )


# compute the ddG by subtracting the average score across the five replicates,
# so (mean total_score of MUT) - (mean total_score of WT)
# then, scale the ddG it down by 0.298, as computed from the Kellogg benchmark
# (see cartddg/benchmarking/benchmarking_analysis.ipynb)
ddg = (df
    [['mutation','run_type','total_score']]
    .groupby(['mutation','run_type'])
    .agg('mean')
    .pivot_table(index='mutation', columns='run_type', values='total_score')
    .assign(ddg = lambda df: (df.MUT-df.WT)*0.298) # scaling factor included, computed from Kellogg benchmark
    .reset_index()
    .rename_axis(None, axis=1)
    [['mutation','ddg']]
)

# merge the computed ddg into the score term dataframe
df = pd.merge(df, ddg, on=['mutation'])

# write out results to a csv
df.to_csv('3gj0_gsp1_ddg_with_scores.csv', index=False, float_format='%g')

# also write a simpler file of just ddg's
df[['mutation','ddg']].drop_duplicates().to_csv('3gj0_gsp1_ddg.csv', index=False, float_format='%g')

