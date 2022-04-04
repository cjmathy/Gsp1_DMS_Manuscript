# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

df_sets <- read_csv('data/Ras_data/Ras_sca_gof_table.csv', show_col_types=F)
df_annot <- read_csv('scripts/Fig4/toxic_annotations.txt', show_col_types=F)
df <- 
    left_join(df_sets, df_annot, by=c('pos'='position')) %>%
    select(gof=-gof_all, -`C-terminal extension`) %>%  
    mutate(
        'I'=ifelse(`Active site regions`==0, F, T),
        'II'=ifelse(`Distal sites affecting switching`==0, F, T),
        'III'=ifelse(`PTM sites`==0, F, T),
        'IV'=ifelse(`Regulator interfaces`==0, F, T)
    ) %>% 
    select(pos, tox, sca, 'gof'=gof_2plus, I, II, III, IV)


df_to_plot <-
    df %>% 
    select(-sca) %>% 
    filter(tox | gof) %>% 
    mutate(set = case_when(
        tox & gof ~ 'Gsp1\ntoxic/GOF\nand HRas\nGOF',
        tox & !gof ~ 'Gsp1\ntoxic/GOF\nonly',
        !tox & gof ~ 'HRas\nGOF\nonly')) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nGOF',
        'HRas\nGOF\nonly')))  %>% 
    select(-tox, -gof) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(names_to='category', values_to='count', cols=c('I','II','III','IV')) %>% 
    mutate(category = factor(category, levels=c('I','II','III','IV')))

# make a dataframe with just the GTPase regions that are aligned with HRas
# i.e. without the Ran/Gsp1 specific extensions defined in Vetter 1999.

active_site_conserved <- c(seq(19,26), seq(41,47), seq(69,77), seq(124,127), seq(152,154))

df_HRas_canonical <-
    df %>% 
    mutate('I' = ifelse(pos %in% active_site_conserved, T, F)) %>% 
    select(-sca) %>% 
    filter(tox | gof) %>% 
    mutate(set = case_when(
        tox & gof ~ 'Gsp1\ntoxic/GOF\nand HRas\nGOF',
        tox & !gof ~ 'Gsp1\ntoxic/GOF\nonly',
        !tox & gof ~ 'HRas\nGOF\nonly')) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nGOF',
        'HRas\nGOF\nonly'))) %>% 
    select(-tox, -gof, -II, -III, -IV) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(names_to='category', values_to='count', cols='I') %>% 
    mutate(category = factor(category, levels=c('I','II','III','IV')))

df_to_plot %>% 
    ggplot(aes(x=category, y=count, fill=category)) +    
    geom_bar(stat='identity', color='black', lwd=0.2, show.legend=F) +
    scale_fill_manual(values=c('lightgray', '#CD2027', '#CD2027', '#CD2027')) +
    geom_bar(stat='identity', show.legend=F, color='black', lwd=0.2, fill='white', data=df_HRas_canonical) +
    facet_grid(rows=vars(set)) +
    theme_custom +
    theme(strip.text.y = element_text(angle = 0, margin = margin(r = 3, l = 3, unit = "pt")))

ggsave('figures/Fig4/Fig4D_vennbar_GOF.pdf', width=1.8, height=2.2)


#### SCA
df_to_plot <-
    df %>% 
    select(-gof) %>% 
    filter(tox | sca) %>% 
    mutate(set = case_when(
        tox & sca ~ 'Gsp1\ntoxic/GOF\nand HRas\nSCA',
        tox & !sca ~ 'Gsp1\ntoxic/GOF\nonly',
        !tox & sca ~ 'HRas\nSCA\nonly')) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nSCA',
        'HRas\nSCA\nonly'))) %>% 
    select(-tox, -sca) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(names_to='category', values_to='count', cols=c('I','II','III','IV')) %>% 
    mutate(category = factor(category, levels=c('I','II','III','IV')))


df_HRas_canonical <-
    df %>% 
    mutate('I' = ifelse(pos %in% active_site_conserved, T, F)) %>% 
    select(-gof) %>% 
    filter(tox | sca) %>% 
    mutate(set = case_when(
        tox & sca ~ 'Gsp1\ntoxic/GOF\nand HRas\nSCA',
        tox & !sca ~ 'Gsp1\ntoxic/GOF\nonly',
        !tox & sca ~ 'HRas\nSCA\nonly')) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nSCA',
        'HRas\nSCA\nonly'
        ))) %>% 
    select(-tox, -sca, -II, -III, -IV) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(names_to='category', values_to='count', cols='I') %>% 
    mutate(category = factor(category, levels=c('I','II','III','IV')))

df_to_plot %>% 
    ggplot(aes(x=category, y=count, fill=category)) +    
    geom_bar(stat='identity', color='black', lwd=0.2, show.legend=F) +
    scale_fill_manual(values=c('lightgray', '#CD2027', '#CD2027', '#CD2027')) +
    geom_bar(stat='identity', show.legend=F, color='black', lwd=0.2, fill='white', data=df_HRas_canonical) +
    facet_grid(rows=vars(set)) +
    theme_custom +
    theme(strip.text.y = element_text(angle = 0, margin = margin(r = 3, l = 3, unit = "pt")))

ggsave('figures/Fig4/Fig4E_vennbar_SCA.pdf', width=1.8, height=2.2)



### Analysis of positions identified by SCA but not toxics

sca_not_tox <-
    df %>% 
    filter(sca, !tox) %>% 
    pull(pos)

read_csv('data/Gsp1_fitness_scores.csv', col_types = cols()) %>% 
    filter(!low_reads_flag) %>% 
    filter(position %in% sca_not_tox) %>% 
    group_by(position, bin) %>% 
    mutate(n = n()) %>% 
    select(position, bin, n) %>% 
    unique() %>% 
    pivot_wider(id_cols=position, names_from=bin, values_from=n) %>% 
    select(toxic, `STOP-like`, intermediate, `WT-like`, beneficial) %>% 
    arrange(desc(toxic)) %>% 
    print(n=30)

