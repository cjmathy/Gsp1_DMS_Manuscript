# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

df_sets <- read_csv('data/Ras_data/Ras_sca_gof_table.csv', show_col_types=F)
df_annot <- read_csv('scripts/Fig4/toxic_annotations.txt', show_col_types=F)
df <- 
    left_join(df_sets, df_annot, by=c('pos'='position')) %>%
    select(-`C-terminal switch`) %>%  
    mutate(
        `GTPase regions`=ifelse(`GTPase regions`==0, F, T),
        `Contacts Nucleotide`=ifelse(`Contacts Nucleotide`==0, F, T),
        `Distal sites affecting switching`=ifelse(`Distal sites affecting switching`==0, F, T),
        PTMs=ifelse(PTMs==0, F, T),
        `Regulator Interface`=ifelse(`Regulator Interface`==0, F, T),
    )


df_to_plot <-
    df %>% 
    select(-sca, -gof_all) %>% 
    filter(tox | gof_2plus) %>% 
    mutate(
        set = case_when(
            tox & gof_2plus ~ 'Gsp1\ntoxic/GOF\nand HRas\nGOF',
            tox & !gof_2plus ~ 'Gsp1\ntoxic/GOF\nonly',
            !tox & gof_2plus ~ 'HRas\nGOF\nonly'
        )
    ) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nGOF',
        'HRas\nGOF\nonly'
        )))  %>% 
    select(-tox, -gof_2plus) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(
        names_to='category',
        values_to='count',
        cols=c(`GTPase regions`,
                `Contacts Nucleotide`,
                `Distal sites affecting switching`,
                PTMs,
                `Regulator Interface`)
    ) %>% 
    mutate(category=gsub(x=category, pattern=' ', replacement='\n')) %>% 
    mutate(category=ifelse(
        category=='Distal\nsites\naffecting\nswitching',
        'Distal sites\naffecting switching',
        category)) %>% 
    mutate(category = factor(category, levels=c(
        'GTPase\nregions',
        'Contacts\nNucleotide',
        'Distal sites\naffecting switching',
        'PTMs',
        'Regulator\nInterface'
    )))

# make a dataframe with just the GTPase regions that are aligned with HRas
# i.e. without the Ran/Gsp1 specific extensions defined in Vetter 1999.

gtpase_sites_conserved <- c(seq(19,26), seq(41,47), seq(69,77), seq(124,127), seq(152,154))

df_HRas_canonical <-
    df %>% 
    mutate(`GTPase regions` = ifelse(pos %in% gtpase_sites_conserved, T, F)) %>% 
    select(-sca, -gof_all) %>% 
    filter(tox | gof_2plus) %>% 
    mutate(
        set = case_when(
            tox & gof_2plus ~ 'Gsp1\ntoxic/GOF\nand HRas\nGOF',
            tox & !gof_2plus ~ 'Gsp1\ntoxic/GOF\nonly',
            !tox & gof_2plus ~ 'HRas\nGOF\nonly'
        )
    ) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nGOF',
        'HRas\nGOF\nonly'
        )))  %>% 
    select(-tox, -gof_2plus, -`Contacts Nucleotide`, -`Distal sites affecting switching`, -PTMs, -`Regulator Interface`) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(
        names_to='category',
        values_to='count',
        cols=c(`GTPase regions`)
    ) %>% 
    mutate(category=gsub(x=category, pattern=' ', replacement='\n')) %>% 
    mutate(category = factor(category, levels=c(
        'GTPase\nregions',
        'Contacts\nNucleotide',
        'Distal sites\naffecting switching',
        'PTMs',
        'Regulator\nInterface'
    )))

df_to_plot %>% 
    ggplot(aes(x=category, y=count, fill=category)) +    
    geom_bar(stat='identity', color='black', lwd=0.2, show.legend=F) +
    scale_fill_manual(values=c('lightgray', 'lightgray', '#CD2027', '#CD2027', '#CD2027')) +
    geom_bar(stat='identity', show.legend=F, color='black', lwd=0.2, fill='white', data=df_HRas_canonical) +
    facet_grid(rows=vars(set)) +
    theme_custom +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_text(angle = 0, margin = margin(r = 3, l = 3, unit = "pt"))
        )

ggsave('figures/Fig4/Fig4D_vennbar_GOF.pdf', width=1.8, height=3)




#### SCA

df_to_plot <-
    df %>% 
    select(-gof_2plus, -gof_all) %>% 
    filter(tox | sca) %>% 
    mutate(
        set = case_when(
            tox & sca ~ 'Gsp1\ntoxic/GOF\nand HRas\nSCA',
            tox & !sca ~ 'Gsp1\ntoxic/GOF\nonly',
            !tox & sca ~ 'HRas\nSCA\nonly'
        )
    ) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nSCA',
        'HRas\nSCA\nonly'
        )))  %>% 
    select(-tox, -sca) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(
        names_to='category',
        values_to='count',
        cols=c(`GTPase regions`,
                `Contacts Nucleotide`,
                `Distal sites affecting switching`,
                PTMs,
                `Regulator Interface`)
    ) %>% 
    mutate(category=gsub(x=category, pattern=' ', replacement='\n')) %>% 
    mutate(category=ifelse(
        category=='Distal\nsites\naffecting\nswitching',
        'Distal sites\naffecting switching',
        category)) %>% 
    mutate(category = factor(category, levels=c(
        'GTPase\nregions',
        'Contacts\nNucleotide',
        'Distal sites\naffecting switching',
        'PTMs',
        'Regulator\nInterface'
    ))) 

df_HRas_canonical <-
    df %>% 
    mutate(`GTPase regions` = ifelse(pos %in% gtpase_sites_conserved, T, F)) %>% 
    select(-gof_2plus, -gof_all) %>% 
    filter(tox | sca) %>% 
    mutate(
        set = case_when(
            tox & sca ~ 'Gsp1\ntoxic/GOF\nand HRas\nSCA',
            tox & !sca ~ 'Gsp1\ntoxic/GOF\nonly',
            !tox & sca ~ 'HRas\nSCA\nonly'
        )
    ) %>%
    mutate(set=factor(set, levels = c(
        'Gsp1\ntoxic/GOF\nonly',
        'Gsp1\ntoxic/GOF\nand HRas\nSCA',
        'HRas\nSCA\nonly'
        )))  %>% 
    select(-tox, -sca, -`Contacts Nucleotide`, -`Distal sites affecting switching`, -PTMs, -`Regulator Interface`) %>% 
    group_by(set) %>% 
    summarise(across(-pos, ~ sum(., is.na(.), 0))) %>% 
    ungroup() %>% 
    pivot_longer(
        names_to='category',
        values_to='count',
        cols=c(`GTPase regions`)
    ) %>% 
    mutate(category=gsub(x=category, pattern=' ', replacement='\n')) %>% 
    mutate(category = factor(category, levels=c(
        'GTPase\nregions',
        'Contacts\nNucleotide',
        'Distal sites\naffecting switching',
        'PTMs',
        'Regulator\nInterface'
    )))


df_to_plot %>% 
    ggplot(aes(x=category, y=count, fill=category)) +    
    geom_bar(stat='identity', color='black', lwd=0.2, show.legend=F) +
    scale_fill_manual(values=c('lightgray', 'lightgray', '#CD2027', '#CD2027', '#CD2027')) +
    geom_bar(stat='identity', show.legend=F, color='black', lwd=0.2, fill='white', data=df_HRas_canonical) +
    facet_grid(rows=vars(set)) +
    theme_custom +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_text(angle = 0, margin = margin(r = 3, l = 3, unit = "pt"))
        )

ggsave('figures/Fig4/Fig4E_vennbar_SCA.pdf', width=1.8, height=3)



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

