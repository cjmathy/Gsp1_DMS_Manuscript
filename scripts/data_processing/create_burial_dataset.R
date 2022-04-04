# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R') 


# read in the file that has SASA info for all of the RAN/Gsp1 complexes
# this is taken from Perica and Mathy et al. 2021 Nature
# https://github.com/tinaperica/Gsp1_manuscript/blob/master/Data/Gsp1_interfaces_SASA_and_conservation.txt
sasa <- read_tsv('data/Gsp1_interfaces_SASA_and_conservation.txt', col_types=cols()) %>% 
  filter(protein == 'GSP1') %>% 
  select(file_path, partner, species, 'pdb_num' = resno, 'pdb_aa' = aa, yeast_seq, interface, deltarASA, deltaASA, rASAc, yeast_num, identical)


# Group residues by SASA
# the goal is to define three regions: buried, interface, surface
# get residues that are buried in every structure
# get all residues (in any of the PDB files) listed to have relative
#   surface accessible area less than 25$
# this gets you all the residues that are buried (in interior or interface) in EVERY STRUCTURE (always buried, never solvent exposed)
monomer_core_residues <- sasa %>% ##### 62 interior residues
  group_by(yeast_num) %>% 
  summarize('max_rASAc' = max(rASAc)) %>% 
  ungroup() %>% 
  filter(max_rASAc < 0.25) %>% 
  pull(yeast_num)

# get residues that form the cores of interaction interfaces
# these are residues that are part of the interface core in at least one of the complexes
interface_core_residues <- sasa %>% ### 87 residues
  filter(interface == 'core') %>% 
  pull(yeast_num) %>% unique()

# are there residue that fit into both categories? always buried and sometimes interface core?
interface_core_residues[interface_core_residues %in% monomer_core_residues]

# yes, the only residue that fits both is residue 47 - which is notably the last residue of switch I
# mostly SURFACE residues
# residues that are always either surface or interface rim, i.e. their rASAc is always above 0.25
surface_rim_residues <- sasa %>% ### 40 residues
  group_by(yeast_num) %>% 
  summarize('min_rASAc' = min(rASAc)) %>% 
  ungroup() %>% 
  filter(min_rASAc > 0.25) %>% 
  pull(yeast_num)

# which residues don't fit into any of the three categories?
mixed_sasa_residues <- sasa %>% 
  filter(! yeast_num %in% c(monomer_core_residues, interface_core_residues, surface_rim_residues)) %>% 
  arrange(yeast_num) %>% pull(yeast_num) %>% unique()

unsorted_residues <- tibble('yeast_num' = 1:219) %>% 
  filter(! yeast_num %in% c(monomer_core_residues, interface_core_residues, surface_rim_residues, mixed_sasa_residues)) %>% 
  arrange(yeast_num) %>% pull(yeast_num) %>% unique()

# add the unsorted ones to surface
surface_rim_residues <- append(surface_rim_residues, unsorted_residues)

# write out SASA group categories to a datafile
sasa_categories <- c('structure core', 'interface core', 'mixed', 'surface')

data.frame(list('position'=seq(1,220))) %>% 
  mutate('sasa_group_category' = case_when(
    position %in% monomer_core_residues ~ 'structure core',
    position %in% interface_core_residues ~ 'interface core',
    position %in% surface_rim_residues ~ 'surface',
    position %in% mixed_sasa_residues ~ 'mixed'
  )) %>% 
  mutate('sasa_group_category' = factor(sasa_group_category, sasa_categories)) %>% 
  filter(position < 220) %>%  # remove the stop codon mutations
  select(position, sasa_group_category) %>% 
  unique() %>% 
  write_csv('data/burial.csv')

