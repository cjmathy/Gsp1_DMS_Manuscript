source('scripts/config_workspace.R')

df_aln <- read_delim('data/Ras_data/HRas_Gsp1_index.txt', delim='\t', show_col_types=F)

df_Gsp1 <-
    read_csv('data/Gsp1_fitness_scores.csv', show_col_types=F) %>% 
    left_join(select(df_aln, aln_num, 'position'=Gsp1_res_num),by='position') %>% 
    select(aln_num,'aa_from_Gsp1'=aa_from, 'position_Gsp1'='position', aa_to,
          'score_Gsp1'=score, 'bin_Gsp1'=bin,'low_reads_flag_Gsp1'=low_reads_flag)

df_HRas <-
    read_csv('data/Ras_data/HRas188_BaF3_processed.csv', show_col_types=F) %>% 
    left_join(select(df_aln, aln_num, 'position'=HRas_res_num),by='position') %>% 
    select(aln_num, 'aa_from_HRas'=aa_from, 'position_HRas'=position, aa_to,
           'score_HRas'=score)

activating_threshold <- 1.5*sd(pull(df_HRas, score_HRas))
activating_threshold

# Positions unique to one protein or the other do not appear in df
# i.e aln_num = 105,106,133 in HRas are insertions compared to Gsp1
# and aln_num = 38 in Gsp1 is not in HRas
# none of the HRas positions 105, 106, 133 have activating mutations
# also the termini that are not included in both experiments, i.e. aln_num = 1
# gets removed because it is present in the Gsp1 dataset but not in the HRas one

df_HRas %>% 
    mutate('TP'=ifelse(score_HRas > activating_threshold, T, F)) %>% 
    group_by(position_HRas, TP) %>% 
    summarise(n=n()) %>% 
    filter(TP) %>% 
    filter(n>1) %>% 
    pull(position_HRas) %>% 
    paste(collapse='+')

