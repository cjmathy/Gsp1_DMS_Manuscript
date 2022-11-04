# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')
    
df_activating <- read_delim('data/Ras_data/n_activating_for_HRas_pos.txt',
                            delim='\t', col_types=cols(), col_names=c('pos_HRas','n_toxics'))

df_toxics <- read_delim('data/Ras_data/n_toxic_for_HRas_pos.txt',
                        delim='\t', col_types=cols(), col_names=c('pos_HRas','n_activating'))

ras_regions <- list(
    `Switch I (G2)\n(30-40)` = seq(30,40,1),
    `Switch II (G3)\n(60-76)` = seq(60,76,1),
    `Helix 3, N-term\n(86-97)` = seq(86,97,1),
    `Helix 3, C-term\n(98-103)` = seq(98,103,1),
    `Loop 7 /\nSite I\n(104-110)` = seq(104,110,1), # site 1 is 103-108 in Gorfe2007
    `Loop 8 (G4) /\nSite II\n(117-126)` = seq(117,126,1), # site 2 is 117-126 in Gorfe2007
    `Loop 10 (G5) /\nSite III\n(144-151)` = seq(144,151,1) # site 3 is 144-151 in Gorfe2007
)

df_regions <-
    enframe(ras_regions) %>% 
    unnest(value) %>% 
    rename('region'=name, 'pos_HRas'=value) %>% 
    mutate(region=factor(region, levels = c(
        'Switch I (G2)\n(30-40)',
        'Switch II (G3)\n(60-76)',
        'Helix 3, N-term\n(86-97)',
        'Helix 3, C-term\n(98-103)',
        'Loop 7 /\nSite I\n(104-110)',
        'Loop 8 (G4) /\nSite II\n(117-126)',
        'Loop 10 (G5) /\nSite III\n(144-151)'
    )))

df_toxics %>% 
    left_join(df_activating, by='pos_HRas') %>% 
    left_join(df_regions, by='pos_HRas') %>%     
    filter(!is.na(region)) %>% 
    pivot_longer(cols=c(n_toxics,n_activating), names_to='type', values_to='n') %>% 
    mutate(type=factor(type, levels=c('n_toxics', 'n_activating'))) %>% 
    group_by(region,type) %>% 
    summarise(n_tot = sum(n)) %>% 
    ggplot(aes(x=region, y=n_tot, fill=type)) +
    geom_bar(stat='identity', position='dodge', color='black', lwd=0.4) +
    scale_fill_manual(
        values=c('black','gray'), 
        labels=c('Activating in H-Ras','Toxic/GOF in Gsp1'),
        ) +
    xlab('Ras region') + ylab('Number of mutations') +
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    theme_custom +
    theme(
        legend.position=c(0.6,0.8),
        legend.box.background = element_rect(colour = 'black'),
        legend.title=element_blank()
    )

ggsave('figures/FigS7/FigS7_bar.pdf', width =6, height = 2)



