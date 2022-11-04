# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# Write a tab-delimited file for the pymol script that makes Fig S7 structures

toxics <-
    read_csv('data/toxic_positions.csv', col_types=cols()) %>% 
    filter(position != 220) %>% 
    arrange(desc(n)) %>% 
    pull(position)

df_tox <- read_delim('data/toxic_positions.csv', col_types=cols(), delim=',')
df_aln <- read_delim('data/Ras_data/HRas_Gsp1_index.txt', col_types=cols(), delim='\t')

df_toxics <- 
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    filter(!low_reads_flag) %>% 
    mutate(is_toxic = ifelse((aa_to != '*') & (bin == 'toxic'), T, F)) %>% 
    group_by(position, is_toxic) %>%
    summarise(n = n(), .groups='drop') %>% 
    complete(position, is_toxic, fill=list(n=0)) %>% 
    filter(is_toxic) %>% 
    select(position, 'n_toxics'=n) %>% 
    left_join(df_aln, by=c('position'='Gsp1_res_num')) %>% 
    select('pos_HRas'=HRas_res_num, n_toxics) %>% 
    filter(!is.na(pos_HRas))

write_delim(df_toxics, 'data/Ras_data/n_toxic_for_HRas_pos.txt', delim='\t', col_names=FALSE)


### now write out activating mutations from Hidalgo et al 2022 eLife

df_HRas <-
    read_csv('data/Ras_data/HRas188_BaF3_processed.csv', show_col_types=F) %>% 
    left_join(select(df_aln, aln_num, 'position'=HRas_res_num),by='position') %>% 
    select(aln_num, 'aa_from_HRas'=aa_from, 'position_HRas'=position, aa_to,
           'score_HRas'=score)

activating_threshold <- 1.5*sd(pull(df_HRas, score_HRas))
activating_threshold


### print list of positions with one or more activating mutation in Ras
df_activating <- 
    df_HRas %>% 
    mutate('is_activating'=ifelse(score_HRas > activating_threshold, T, F)) %>%
    group_by(position_HRas, is_activating) %>% 
    summarise(n = n(), .groups='drop') %>% 
    complete(position_HRas, is_activating, fill=list(n=0)) %>% 
    filter(is_activating) %>% 
    select('pos_HRas'=position_HRas, 'n_activating'=n)
    
write_delim(df_activating, 'data/Ras_data/n_activating_for_HRas_pos.txt', delim='\t', col_names=FALSE)

# output is 4+8+12+13+14+16+18+19+20+22+23+24+34+56+58+59+60+61+63+68+73+77+79+81+83+84+85+92+116+117+118+119+120+127+136+138+144+145+146+147+148+152+155+156+157
# rewritten for pymol is 4+8+12-14+16+18-20+22-24+34+56+58-61+63+68+73+77+79+81+83-85+92+116-120+127+136+138+144-148+152+155-157

