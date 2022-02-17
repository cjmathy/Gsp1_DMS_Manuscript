import pyrosetta
from pyrosetta import *
pyrosetta.init()

gsp1 = pose_from_pdb('/Users/cjmathy/gdrive/gsp1_dms/data/pdbs/pdbs_ran/3m1i.pdb')

def get_heavy_atoms(res):
    return [atom for i, atom in enumerate(res.atoms()) 
            if 'H' not in res.atom_name(i+1)]

def residues_contact(res1, res2, cutoff):
    for atom1 in get_heavy_atoms(res1):
        for atom2 in get_heavy_atoms(res2):
            if (atom1.xyz() - atom2.xyz()).norm() < cutoff:
                return True    
    return False

nucleotide = gsp1.residue(201)
magnesium = gsp1.residue(202)
contacting_residues = []

for residue in gsp1.residues:
    if residues_contact(residue, nucleotide, cutoff=4) or residues_contact(residue, magnesium, cutoff=4):
        pos = int(gsp1.pdb_info().pose2pdb(residue.seqpos()).split()[0])
        if pos < 221:  # excludes the nucleotide and Mg
            contacting_residues.append(pos)

with open('/Users/cjmathy/gdrive/gsp1_dms/data/residues_contacting_nucleotide.txt','w') as f:
    f.write('# Residues numbers are from sequence numbering, not PDB numbering.\n')
    f.write('# Distances were computed using PDB 3m1i, chain A.\n')
    f.write('# Contacting residues determined by a 4 Å cutoff\n')
    for res in contacting_residues:
        f.write(str(res)+'\n')