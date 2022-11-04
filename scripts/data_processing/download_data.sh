#!/bin/bash
# Download Ran PDBs from PDB-REDO

wget https://pdb-redo.eu/db/1i2m/1i2m_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/1k5d/1k5d_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/1wa5/1wa5_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/2bku/2bku_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3a6p/3a6p_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3ch5/3ch5_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3ea5/3ea5_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3gj0/3gj0_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3icq/3icq_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3m1i/3m1i_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3nc1/3nc1_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3w3z/3w3z_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/3wyf/3wyf_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/4hb2/4hb2_final.pdb -P ../data/pdbs_raw/
wget https://pdb-redo.eu/db/4ol0/4ol0_final.pdb -P ../data/pdbs_raw/
wget https://files.rcsb.org/download/1a2k.pdb -P ../data/pdbs_raw/
wget https://files.rcsb.org/download/1ibr.pdb -P ../data/pdbs_raw/
wget https://files.rcsb.org/download/1qbk.pdb -P ../data/pdbs_raw/
wget https://files.rcsb.org/download/3k8y.pdb -P ../data/pdbs_raw/
wget https://files.rcsb.org/download/2rge.pdb -P ../data/pdbs_raw/

rename 's/_final//' ../data/pdbs_raw/*_final*

# Download FASTA sequence files for S. cerevisiae Gsp1 and human Ran

wget https://www.uniprot.org/uniprot/P32835.fasta -O ../data/gsp1.fasta
wget https://www.uniprot.org/uniprot/P62826.fasta -O ../data/ran.fasta

