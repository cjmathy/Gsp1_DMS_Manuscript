source('scripts/config_workspace.R')


# Data S1: Gsp1 Fitness Scores
df_fit <- read_csv('data/Gsp1_fitness_scores.csv', col_types=cols())
df_reads <- read_csv('data/raw_counts.csv', col_types=cols())

inner_join(df_reads, df_fit, by=c('mutant','aa_from','position','aa_to')) %>% 
write_csv('data/DataS1_Gsp1_Reads_Scores_Bins.csv')

# Data S2: Gsp1 ddG data
df_seq <-
    read_csv('data/Gsp1_Ran_sequence_alignment.csv', col_types=cols()) %>% 
    mutate(pos_Sc = as.integer(pos_Sc), pos_Hs = as.integer(pos_Hs))
df_ddg_3m1i <-
    read_csv('data/ddg/3m1i_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Sc', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(pdb_id='3m1i', species='Sc') %>% 
    right_join(df_seq, by='pos_Sc')
df_ddg_3gj0 <-
    read_csv('data/ddg/3gj0_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Hs', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(pdb_id='3gj0', species='Hs') %>% 
    right_join(df_seq, by=c('pos_Hs'))
bind_rows(df_ddg_3m1i, df_ddg_3gj0) %>% 
    select(mutation,aa_from,aa_to,pdb_id,species,pos_Sc,aa_Sc,pos_Hs,aa_Hs,ddg) %>% 
    filter(mutation!='NA') %>% 
    write_csv('data/DataS2_Gsp1_DDG.csv')
  
# Data S2: Gsp1 benchmark data
read_csv('data/ddg/ddg_benchmark_dataset.csv', col_types=cols()) %>% 
    mutate(position = as.integer(substring(mutation_short, 2, nchar(mutation_short)-1)),
           aa_from = substring(mutation_short, 1, 1),
           aa_to = substring(mutation_short, nchar(mutation_short), nchar(mutation_short))) %>% 
    mutate(ddg_expt = case_when(
              is.na(ddg_expt_kellogg) ~ ddg_expt_protherm,
              is.na(ddg_expt_protherm) ~ ddg_expt_kellogg,
              T ~ ddg_expt_protherm)) %>% 
    mutate(ddg_calc_adj = ddg_calc*0.298) %>% 
    select(record_ID, pdb, 'mutation_full'=mutation, chain, 'mutation'=mutation_short, position, ddg_expt, ddg_calc, ddg_calc_adj) %>% 
    write_csv('data/DataS3_Benchmark_DDG.csv')
