library(tidyverse)

df_seq <- read_csv('data/Gsp1_Ran_sequence_alignment.csv', col_types=cols())

conserved_positions <-
    df_seq %>%
    filter(aa_Sc==aa_Hs) %>%
    mutate(pos_Sc = as.numeric(pos_Sc)) %>% 
    pull(pos_Sc)

# 82% identity
length(conserved_positions) / 219

df_fit <- 
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>% 
    filter(!low_reads_flag) %>% 
    mutate(is_toxic = ifelse(bin == 'toxic', T, F)) %>% 
    mutate(is_conserved = ifelse(position %in% conserved_positions, T, F))
    
df1 <- 
    df_fit %>% 
    filter(is_toxic) %>% 
    group_by(position) %>% 
    mutate(n_tox = n()) %>% 
    select(position, is_conserved, n_tox) %>% 
    unique %>% 
    filter(is_conserved)

# 136 of 179 (76%) of conserved positions have at least 1 toxic
# And so 43/179 (24%) have no toxics
nrow(df1)

# 117 of 179 (65%) of conserved positions have at least 2 toxics
df1 %>% filter(n_tox > 1) %>% nrow



### Which conserved positions have no toxics?
conserved_no_toxics <-
    df_fit %>% 
    group_by(position) %>% 
    mutate(n=n()) %>% # need this since we removed low reads mutants
    filter(!is_toxic) %>% 
    mutate(n_nottox = n()) %>% 
    select(position, is_conserved, n_nottox, n) %>% 
    unique %>% 
    filter(is_conserved, n_nottox==n) %>% 
    pull(position)
conserved_no_toxics


# Do they have many STOP-like's?

df_fit %>% 
    filter(position %in% conserved_no_toxics) %>% 
    group_by(position, bin) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from='bin', values_from='n') %>% 
    print(n=49)

# no, most are WT-like



# finally, of the 60 positions with at least 10 toxics, how many are conserved

toxic_positions <-
    read_csv('data/toxic_positions.csv', col_types=cols()) %>% 
    filter(position != '220') %>% 
    pull(position)



conserved_positions

toxic_positions

length(intersect(toxic_positions, conserved_positions))







