source('scripts/config_workspace.R')

sample_order <- c('WT','F28V','F28Y','H50N','H50R','F54A','F54W',
                  'N156A','N156W','F159L','F159W','F163L','F163Y')

# file made in scripts/Fig3/Fig3D_GEF.R
df <- read_csv('data/GEF_assay_nucleotide_preferences.csv', show_col_types = FALSE) 



df %>% 
  mutate(sample=factor(sample, levels=sample_order)) %>% 
  ggplot(aes(x=sample, y=r)) + 
  geom_bar(aes(fill=type), stat='identity') +
  geom_errorbar(aes(ymin=r-sd,ymax=r+sd), width=0.3, alpha=0.7) +
  xlab('Gsp1 mutant') + ylab('Inital rate (absolute value of v0),\nrelative to WT') +
  xlab('Fitness score') + ylab('Relative preference for GTP\n(ratio of v0, GTP/GDP)') + 
  scale_fill_manual(values=c('red2','gray','blue2'), name='Mutant Type') +
  theme_custom +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = 'none')

ggsave('figures/FigS7/FigS7C_ratio_bar.pdf', height=2, width=3.5)




