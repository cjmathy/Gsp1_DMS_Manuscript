# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')
library(ggrepel)

df_annot <- read_csv('scripts/Fig4/toxic_annotations.txt', show_col_types=F)

df_dists <-
    read_csv('data/dist_to_nuc.txt', show_col_types=F) %>% 
    rename('dist'=dist_gtp) %>% 
    select(-dist_gdp)


df_n_tox <- 
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    filter(!low_reads_flag) %>% 
    filter(aa_to != '*') %>% 
    group_by(position, bin) %>% 
    summarise(n_tox=n(), .groups='drop') %>% 
    complete(position, bin, fill = list(n_tox = 0)) %>% 
    filter(bin=='toxic') %>% 
    select(-bin)

df_label <-
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    group_by(mutant) %>% 
    mutate(label = str_c(aa_from, position)) %>% 
    ungroup() %>% 
    select(position, label) %>% 
    unique()

df <- 
    left_join(df_annot, df_dists, by=c('position'='pos')) %>% 
    left_join(df_n_tox, by='position') %>% 
    left_join(df_label, by='position') %>% 
    filter(!is.na(dist)) %>% 
    filter(position < 181) %>% 
    mutate(cat = ifelse(`Active site regions`==1, 'Active site', 'Distal sites')) %>% 
    mutate(cat = case_when(
        `Distal sites affecting switching`==1 ~ 'Distal sites affecting switching',
        `PTM sites` == 1 ~ 'PTM sites',
        `Regulator interfaces`==1 ~ 'Regulator interfaces', 
        `Active site regions` == 1 ~ 'Active site regions',
        T ~ 'None'
    )) %>% 
    mutate(cat = factor(cat, levels = c(
        'Active site regions',
        'Distal sites affecting switching',
        'PTM sites',
        'Regulator interfaces',
        'None'
    ))) %>% 
    mutate(label = ifelse((n_tox > 9.5) & (cat != 'Active site regions'), label, NA))

# number of toxics further than 5 angstroms from the nucleotide
df %>% filter(n_tox >=10, dist > 5.0, `Active site regions`!=1)


# distances of novel allosteric cluster
df %>% 
    filter(position %in% c(28,50,54,156,159,163)) %>% 
    select(position, dist)

df %>% 
    filter(n_tox > 4, n_tox < 10) %>% 
    filter(!cat == 'Active site regions') %>% 
    arrange(desc(dist), n_tox) %>% 
    filter(position %in% c(29,30,31,49,51,135,162)) %>% 
    print(n=22) 

    
df %>% filter(n_tox > 7, n_tox < 10, dist > 20)


df %>% filter(position %in% c(56,59,176,178))

df %>% filter


T56E



df %>% 
    ggplot(aes(x = dist, y = n_tox)) + 
    annotate('rect', xmin = 0, xmax = 5, ymin = 9.5, ymax = 20, alpha = 0.6, fill = '#c3dadb') +
    annotate('rect', xmin = 0, xmax = 5, ymin = -1, ymax = 9.5, alpha = 0.6, fill = '#c3dadb') +
    annotate('rect', xmin = 5, xmax = 31, ymin = -1, ymax = 9.5, alpha = 0.6, fill = '#c3dadb') +
    geom_vline(xintercept=5, lty='dotted') +
    geom_hline(yintercept=9.5, lty='dotted') +
    annotate('label', x = 23, y = 17.5, label = 'Allosteric functional\npositions', size=3) +
    geom_point(shape=21, size = 3, alpha=0.6, stroke=0.6, color='gray', fill='#E5E5E5', data=filter(df, cat=='None')) +
    geom_point(aes(fill=cat, color=cat), pch=19, alpha=1, size = 3, data=filter(df, cat=='Active site regions')) +
    geom_point(aes(fill=cat, color=cat), pch=19, alpha=1, size = 3, data=filter(df, ! cat %in% c('Active site regions', 'None'))) +
    geom_text_repel(
        aes(label = label, color=cat), show.legend=F, size=2,
        nudge_x = 0.2, nudge_y = 0.1,
        box.padding = 0.3, min.segment.length = 0.3, segment.size = 0.2
        ) + 
    scale_color_manual(
        name='Functional annotation',
        values=c('gray70','#EE3530','#BF40BF','#1ABCBD', '#e8e8e8'),
        breaks=c(
        'Active site regions',
        'Distal sites affecting switching',
        'PTM sites',
        'Regulator interfaces',
        'None')) +
    scale_x_continuous(limits = c(0,31), expand = c(0,0)) +
    scale_y_continuous(limits = c(-1,20), expand = c(0,0)) +
    xlab('Closest sidechain atom to nucleotide (Å)') +
    ylab('Number of toxic mutations') +
    theme_custom +
    theme(legend.position='none')

ggsave('figures/Fig4/Fig4C_dist_to_nuc.pdf', width=3.5, height=3.5)



