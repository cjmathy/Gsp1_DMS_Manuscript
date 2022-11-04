import pyrosetta
import pandas as pd

# Global variables

class AminoAcidList():
    def __init__(self):
        self.ONELETTER = ['W','Y','F','M','L','I','V','A','C','G',
                          'P','T','S','Q','N','E','D','H','R','K']
        self.oneletter = ['w','y','f','m','l','i','v','a','c','g',
                          'p','t','s','q','n','e','d','h','r','k']
        self.THREELETTER = ['TRP', 'TYR', 'PHE', 'MET', 'LEU',
                            'ILE', 'VAL', 'ALA', 'CYS', 'GLY',
                            'PRO', 'THR', 'SER', 'GLN', 'ASN',
                            'GLU', 'ASP', 'HIS', 'ARG', 'LYS']
        self.threeletter = ['Trp', 'Tyr', 'Phe', 'Met', 'Leu',
                            'Ile', 'Val', 'Ala', 'Cys', 'Gly',
                            'Pro', 'Thr', 'Ser', 'Gln', 'Asn',
                            'Glu', 'Asp', 'His', 'Arg', 'Lys']
    
    def __str__(self):
        lines = ['ONELETTER', ', '.join(self.ONELETTER),
                 'oneletter', ', '.join(self.oneletter),
                 'THREELETTER', ', '.join(self.THREELETTER),
                 'threeletter', ', '.join(self.threeletter)]
        return('\n'.join(lines))
    
    def __repr__(self):
        return(print(self))
        
    
aminoacids = AminoAcidList()

'''Functions for use in pyrosetta'''

def pose_from_chains(pose, chains):
    ''' Extracts a new pose after including only residues with certain
    chain labels

    Keywork arguments:
    pose -- a Pose object
    chains -- a list of letter corresponding to the chain in the pose to be extracted
    '''

    r = [] # residues to keep

    for res in pose.residues:
        num = res.seqpos()
        chainletter = pose.pdb_info().chain(num)

        if chainletter in chains:
            r.append(num)

    return (pyrosetta.rosetta.core.pose.Pose(pose, r[0], r[-1]))

def pose_aa_only(pose):
    '''Extracts a new pose containing only amino acid Residue objects'''
    r = [res.seqpos() for res in pose.residues if res.is_protein()]
    return (pyrosetta.rosetta.core.pose.Pose(pose, r[0], r[-1]))
