# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')


# load fitness data
df_fit <- read_csv('data/Gsp1_fitness_scores.csv', col_types = cols())


df_seq <-
    read_csv('data/Gsp1_Ran_sequence_alignment.csv', col_types=cols()) %>% 
    mutate(pos_Sc = as.integer(pos_Sc), pos_Hs = as.integer(pos_Hs))

df_ddg_3m1i <-
    read_csv('data/ddg/3m1i_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Sc', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(struct='3m1i', species='Sc') %>% 
    right_join(df_seq, by='pos_Sc')

df_ddg_4hb2 <-
    read_csv('data/ddg/4hb2_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Hs', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(struct='4hb2', species='Hs') %>% 
    right_join(df_seq, by=c('pos_Hs'))

df_ddg_3gj0 <-
    read_csv('data/ddg/3gj0_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Hs', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(struct='3gj0', species='Hs') %>% 
    right_join(df_seq, by=c('pos_Hs'))

## ANALYSIS 1: Compare 3m1i with 4hb2
# Compare the two Gsp1 GTP structures to see whether we can use
# ddGs calculated on a human structure
# Note that we filter out positions not included in all structures,
# as well as the catalytic glutamate 71 (human 69), which is mutated 
# to an L in the 3m1i structure.
# This leaves us comparing 172 residues, im yeast numbering 11-70, 72-183
df_ddg_3m1i_4hb2 <-
    bind_rows(df_ddg_3m1i, df_ddg_4hb2) %>% 
    filter(!is.na(struct)) %>%
    filter(pos_Sc != 71) %>%
    select(struct, aa_to, pos_Sc, pos_Hs, aa_Sc, aa_Hs, ddg) %>% 
    pivot_wider(names_from=struct, values_from=ddg, names_prefix='ddg_') %>% 
    filter(!is.na(ddg_3m1i) & !is.na(ddg_4hb2))

df_ddg_3m1i_4hb2 %>% 
    ggplot(aes(x=ddg_3m1i, y=ddg_4hb2)) +
    geom_point() +
    theme_custom

## ANALYSIS 2: Compare 4hb2 with 3gj0
# Compare the two human structures to understand the difference in ddG between states
# This time we can include Q69, which is not mutated in either structure
# This leaves us comparing 173 residues, im yeast numbering 11-183
df_ddg_4hb2_3gj0 <- 
    bind_rows(df_ddg_4hb2, df_ddg_3gj0) %>% 
    filter(!is.na(struct)) %>%
    select(struct, aa_to, pos_Sc, pos_Hs, aa_Sc, aa_Hs, ddg) %>% 
    pivot_wider(names_from=struct, values_from=ddg, names_prefix='ddg_') %>% 
    filter(!is.na(ddg_4hb2) & !is.na(ddg_3gj0))
    
df_ddg_4hb2_3gj0 %>% 
    ggplot(aes(x=ddg_3gj0, y=ddg_4hb2)) +
    geom_point() +
    theme_custom


# merge all three datasets and save

bind_rows(df_ddg_3m1i, df_ddg_3gj0, df_ddg_4hb2) %>% 
    select(mutation,aa_from,aa_to,struct,species,pos_Sc,aa_Sc,pos_Hs,aa_Hs,ddg) %>% 
    left_join(
        select(df_fit, position, aa_to, 'fitness'=score, bin, fitness_low_reads_flag=low_reads_flag),
        by=c('pos_Sc'='position','aa_to'='aa_to')) %>% 
    select(mutation, aa_to, struct, pos_Sc, ddg, fitness, bin, fitness_low_reads_flag) %>% 
    filter((pos_Sc >= 11) & (pos_Sc <= 183)) %>% # shared residues across the PDBs
    pivot_wider(id_cols=c('pos_Sc','aa_to', fitness, bin, fitness_low_reads_flag), names_from='struct', values_from='ddg') %>% 
    mutate('ddG_Gsp1_GTP'=`3m1i`, 'ddG_Ran_GDP'=`3gj0`, 'ddG_Ran_GTP'=`4hb2`) %>% 
    write_csv('data/ddg/Gsp1_ddG_merged_all_states_shared_residues.csv')
  