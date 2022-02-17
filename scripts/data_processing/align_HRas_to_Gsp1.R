# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# load additional packages
library(bio3d)

# this script needs an install of the package 'muscle'
# which can be installed using conda
path_to_muscle <- '/Users/cjmathy/miniconda3/envs/gsp1_dms/bin/muscle'

## Read in the DMS data for human HRas from Badaru et al.
## this file contains the WT sequence of human HRas
HRas_index_sequence <- 
  read_tsv('data/Ras_data/HRas_seq_index.txt', col_types = cols(.default = "c")) %>% 
  pivot_longer(everything(), names_to = 'pos', values_to = 'wt') %>% 
  mutate('pos' = as.integer(pos))

## note: use PDB: 3L8Z as a reference wt human HRas structure (GNP bound)
#### Structural alignment of human HRas and Gsp1 with bio3d
pdb_dir <- 'data/Ras_data/raw_pdbs/'
files_list <- list.files(path = pdb_dir)
files_list <- files_list[grepl('pdb', files_list)]
file_paths <-  file.path(pdb_dir, files_list)
unq.ids <- c('3m1i', '1k5d', '1wq1', '3l8z')
# Download PDBs

raw.files <- get.pdb(unq.ids, path=pdb_dir)

# Extract and align sequences 

pdbs <- pdbaln(raw.files, exefile=path_to_muscle)
pdbs <- pdbaln(file_paths, exefile=path_to_muscle)
core <- core.find(pdbs)
core.inds <- print(core, vol = 0.5)
xyz <- pdbfit(pdbs, core.inds, outpath = 'data/Ras_data/Ran_HRas_aligned')

# get the sequence alignment
alignment <- pdbs$ali %>% 
    as.data.frame() %>%
    rownames_to_column(var = 'pdb_path') %>%
    as_tibble() %>% 
    filter(grepl('3l8z', pdb_path) | grepl('3m1i', pdb_path)) %>% 
    pivot_longer(-pdb_path, names_to = 'aln_num', values_to = 'aa') %>%
    mutate('aln_num' = str_sub(aln_num, start = 2),
           'pdb' = case_when(
               grepl('3l8z', pdb_path) ~ 'HRas',
               grepl('3m1i', pdb_path) ~ 'Gsp1')
          ) %>%
  mutate('aln_num' = as.integer(aln_num)) %>%
  pivot_wider(-pdb_path, names_from = pdb, values_from = aa)

## get residue numbers from PDB file
resnumbers <- pdbs$resno %>% 
    as.data.frame() %>%
    rownames_to_column(var = 'pdb_path') %>%
    as_tibble() %>% 
    filter(grepl('3l8z', pdb_path) | grepl('3m1i', pdb_path)) %>% 
    pivot_longer(-pdb_path, names_to = 'aln_num', values_to = 'res_num') %>%
  mutate('aln_num' = str_sub(aln_num, start = 2),
         pdb = case_when(
           grepl('3l8z', pdb_path) ~ 'HRas_res_num',
           grepl('3m1i', pdb_path) ~ 'Gsp1_res_num'
         )) %>%
  mutate('aln_num' = as.integer(aln_num)) %>%
  pivot_wider(-pdb_path, names_from = pdb, values_from = res_num)


# merge the alignment and residue info into an alignment index for HRas and Gsp1
# Institute a couple of fixes:
#     - Fix the Gsp1 sequence so that position 71 is a Q (there is a Q71L 
#       mutation in the xtal structure)
#     - HRas 164 RQH 166 are misaligned, they are contigious with the previous 
#       residues and should be at alignment positions 165-167, not 191-193
#     - Some gaps have an "NA" for the Gsp1 residue instead of a '-'
alignment <- inner_join(alignment, resnumbers, by = c("aln_num")) %>% 
    filter(aln_num < 168) %>%
    mutate('Gsp1' = ifelse(Gsp1_res_num==71, 'Q', Gsp1)) %>% 
    mutate('HRas' = case_when(
               aln_num == 165 ~ 'R',
               aln_num == 166 ~ 'Q',
               aln_num == 167 ~ 'H',
               TRUE ~ HRas)) %>%
    mutate('HRas_res_num' = case_when(
               aln_num == 165 ~ 164L,
               aln_num == 166 ~ 165L,
               aln_num == 167 ~ 166L,
               TRUE ~ HRas_res_num)) %>% 
    mutate('Gsp1'=ifelse(is.na(Gsp1), '-', Gsp1))

write_tsv(alignment, 'data/Ras_data/HRas_Gsp1_index.txt')


# write fast file for computing sequence identity and sequence similarity
alignment <- read_tsv('data/Ras_data/HRas_Gsp1_index.txt', col_types=cols())
Gsp1_sequence_aligned <- pull(alignment, Gsp1) %>% paste0(collapse='')
HRas_sequence_aligned <- pull(alignment, HRas) %>% paste0(collapse='')

fileConn<-file('data/Ras_data/Gsp1_HRas_aligned.fasta')
writeLines(
  c('>Gsp1',Gsp1_sequence_aligned,'>HRas',HRas_sequence_aligned),
  sep='\n',
  fileConn
)
close(fileConn)

# compute sequence identity and similarity
# http://thegrantlab.org/bio3d/reference/seqidentity.html
aln <- read.fasta('data/Ras_data/Gsp1_HRas_aligned.fasta')
aln
seqidentity(aln, normalize=TRUE, similarity=FALSE)['Gsp1','HRas']
seqidentity(aln, normalize=TRUE, similarity=TRUE)['Gsp1','HRas']
