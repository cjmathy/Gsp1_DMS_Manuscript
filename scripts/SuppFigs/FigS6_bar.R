# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

df_sets <- read_csv('data/Ras_data/Ras_sca_gof_table.csv', show_col_types=F)
df_annot <- read_csv('scripts/Fig4/toxic_annotations.txt', show_col_types=F)

df_n_tox <- 
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    filter(!low_reads_flag) %>% 
    filter(aa_to != '*') %>% 
    group_by(position, bin) %>% 
    summarise(n_tox=n(), .groups='drop') %>% 
    complete(position, bin, fill = list(n_tox = 0)) %>% 
    filter(bin=='toxic') %>% 
    select(-bin)

df <- 
    left_join(df_sets, df_annot, by=c('pos'='position')) %>%
    select(gof=-gof_all, -`C-terminal extension`) %>%  
    mutate(category = ifelse(`Active site regions`==1, 'Active\nsite', 'Distal\nsites')) %>% 
    mutate(category = factor(category, levels=c('Active\nsite', 'Distal\nsites'))) %>% 
    select(pos, tox, sca, 'gof'=gof_2plus, category) %>% 
    left_join(df_n_tox, by=c('pos'='position')) %>% 
    mutate(`Number of toxics` = case_when(
        n_tox >= 10 ~ '10 or more',
        n_tox >= 5 ~ 'at least 5',
        n_tox >= 1 ~ 'at least 1',
        T ~ 'None'
    )) %>% 
    mutate(`Number of toxics` = factor(`Number of toxics`, levels=c(
        '10 or more','at least 5','at least 1','None')))

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
    group_by(set, category, `Number of toxics`) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% 
    complete(set, category, `Number of toxics`, fill = list(count = 0))

df_to_plot %>% 
    ggplot(aes(x=category, y=count, fill=`Number of toxics`)) +    
    geom_bar(stat='identity', position='stack', color='black', lwd=0.2, show.legend=T) +
    scale_fill_manual(values=c('gray40','gray70','lightgray','white')) +
    facet_grid(rows=vars(set)) +
    scale_y_continuous(limits=c(0,17), breaks=c(0,5,10,15)) +
    ylab('Position count') + xlab('') +
    theme_custom +
    theme(strip.text.y = element_text(angle = 0, size = 5, margin = margin(r = 2, l = 2, unit = "pt")))

ggsave('figures/FigS6/FigS6_vennbar_GOF.pdf', width=2.5, height=2.5)


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
        'HRas\nSCA\nonly')))  %>% 
    select(-tox, -sca) %>% 
    group_by(set, category, `Number of toxics`) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% 
    complete(set, category, `Number of toxics`, fill = list(count = 0))



df_to_plot %>% 
    ggplot(aes(x=category, y=count, fill=`Number of toxics`)) +    
    geom_bar(stat='identity', position='stack', color='black', lwd=0.2, show.legend=T) +
    scale_fill_manual(values=c('gray40','gray70','lightgray','white')) +
    facet_grid(rows=vars(set)) +
    scale_y_continuous(limits=c(0,26), breaks=c(0,10,20)) +
    ylab('') + xlab('') +
    theme_custom +
    theme(strip.text.y = element_text(angle = 0, size = 5, margin = margin(r = 2, l = 2, unit = "pt")))

ggsave('figures/FigS6/FigS6_vennbar_SCA.pdf', width=2.5, height=2.5)



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

