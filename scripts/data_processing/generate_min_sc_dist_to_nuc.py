import numpy as np
import pandas as pd
import pyrosetta
from pyrosetta import *
pyrosetta.init()


def get_heavy_atoms(res):
    return [atom for i, atom in enumerate(res.atoms()) 
            if is_heavy_atom(res.atom_name(i+1))]

def is_heavy_atom(atom_name):
    if 'H' in atom_name:
        return False
    return True

def get_sc_heavy_atoms(res):
    return [atom for i, atom in enumerate(res.atoms()) 
            if is_sc_heavy_atom(res.atom_name(i+1), res)]

def is_sc_heavy_atom(atom_name, res):
    if res.name() == 'GLY':
       if atom_name.strip() == 'CA':
           return True
    if atom_name.strip() in ['N','CA','C','O']:
        return False
    if 'H' in atom_name:
        return False
    return True

def closest_sc_heavy_atom_dist(residue, nucleotide):
    closest = 100
    for atom1 in get_sc_heavy_atoms(residue):
        for atom2 in get_heavy_atoms(nucleotide):
            dist = (atom1.xyz() - atom2.xyz()).norm()
            if dist < closest:
                closest = dist
    return closest

def get_sc_distances_to_nucleotide(protein, nucleotide):
    positions = []
    distances = []
    for residue in protein.residues:
        pos = int(protein.pdb_info().pose2pdb(residue.seqpos()).split()[0])
        if (pos < 220):
            positions.append(pos)
            distances.append(closest_sc_heavy_atom_dist(residue, nucleotide))
    return pd.DataFrame({'pos': positions,'dist': distances})

# GTP
gsp1_gtp = pose_from_pdb('/Users/cjmathy/gdrive/gsp1_dms/data/pdbs/pdbs_ran/3m1i.pdb')
GTP = gsp1_gtp.residue(201)
df_gtp = get_sc_distances_to_nucleotide(gsp1_gtp, GTP)

# GDP
ran_gdp = pose_from_pdb('/Users/cjmathy/gdrive/gsp1_dms/data/pdbs/pdbs_ran/3gj0.pdb')
GDP = ran_gdp.residue(209)
df_gdp = get_sc_distances_to_nucleotide(ran_gdp, GDP)

# Convert Ran and Gsp1 sequences
# only use sequence positions 10-180
df_aln = pd.read_csv('/Users/cjmathy/gdrive/gsp1_dms/data/Gsp1_Ran_sequence_alignment.csv')
df_aln = df_aln[df_aln['pos_Sc']!='-']
df_aln = df_aln[df_aln['pos_Hs']!='-']
df_aln = df_aln.astype({'pos_Sc': 'int', 'pos_Hs': 'int'})

df_gtp = df_gtp.astype({'pos':'Int64'})
df_gdp = df_gdp.astype({'pos':'Int64'})

df = (
    df_aln
    .merge(df_gtp, how='left', left_on='pos_Sc', right_on='pos')
    .rename(columns={'pos':'pos_gtp', 'dist':'dist_gtp'})
    .merge(df_gdp, how='left', left_on='pos_Hs', right_on='pos')
    .rename(columns={'pos':'pos_gdp', 'dist':'dist_gdp'})
    [['pos_Sc', 'dist_gtp', 'dist_gdp']]
    .rename(columns={'pos_Sc':'pos', 'dist':'dist_gdp'})
)

df.to_csv('/Users/cjmathy/gdrive/gsp1_dms/data/dist_to_nuc.txt', index=False)
