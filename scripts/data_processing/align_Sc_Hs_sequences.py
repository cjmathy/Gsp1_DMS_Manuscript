import pandas as pd
from Bio.PDB import *
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")


seq_Sc = str(SeqIO.read('/Users/cjmathy/gdrive/gsp1_dms/data/gsp1.fasta', 'fasta').seq)
seq_Hs = str(SeqIO.read('/Users/cjmathy/gdrive/gsp1_dms/data/ran.fasta', 'fasta').seq)


alignment = pairwise2.align.globalds(seq_Sc, seq_Hs, blosum62, -11, -1) # default param on BLAST

# this alignment is pretty close, but from comparing structures (e.g. 3m1i with 4ol0)
# it is pretty clear that Gsp1 10 VPT 13 aligns with Ran 8 QVQ 10. So the sequence for
# Gsp1 in this alignment should start with 'MSAPAANGE-VPTFKL' not 'MSAPAANGEVPT-FKL'

seq_Sc_aligned = 'MSAPAANGE-VPT' + alignment[0][0][13:]
seq_Hs_aligned = alignment[0][1]

df_seq = pd.DataFrame(list(zip(seq_Sc_aligned,
                               seq_Hs_aligned,
                               list(range(1,10)) + ['-'] + list(range(10,220)),
                               [1] + ['-','-','-'] + list(range(2,215)) + ['-',215,216],
                              )), columns=['aa_Sc','aa_Hs', 'pos_Sc', 'pos_Hs'])

# remove the row for the one insertion at human R
df_seq.to_csv('/Users/cjmathy/gdrive/gsp1_dms/data/Gsp1_Ran_sequence_alignment.csv',
               index=False)