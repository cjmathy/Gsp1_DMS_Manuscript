source('scripts/config_workspace.R')

compute_ratio_sd <- function(a,b,sd_a,sd_b) {
  # function for the formula to propagate standards
  # deviations across the division function
  # f = a/b
  # sd_f = |f| * sqrt( (sd_a/a)^2 + (sd_b/b)^2 )
  abs(a/b) * sqrt( (sd_a/a)^2 + (sd_b/b)^2 )
}

sample_order <- c('WT','F28V','F28Y','F54A','F54W','F159L','F159W')

df <-
  read_delim('labnotebook/analysis/GAP_computed_kcat_Km/GAP_kinetics_MichaelisMenten_parameters.txt', delim='\t', col_types=cols()) %>% 
  mutate(mutant_type = case_when(
    mutant %in% c('WT') ~ 'WT',
    mutant %in% c('F28V','H50N','F54A','N156A','F159L','F163L') ~ 'Toxic/GOF',
    mutant %in% c('F28Y','H50R','F54W','N156W','F159W','F163Y') ~ 'Tolerant'
  )) %>% 
  mutate(mutant = factor(mutant, levels = sample_order)) %>% 
  filter(mutant %in% sample_order) %>% 
  filter(mutant != 'F163L') %>% 
  filter(! ((mutant == 'F28V') & (individual_kcat_Km > 300))) %>% 
  mutate(mean_kcat_Km = ifelse(mutant == 'F28V', individual_kcat_Km, mean_kcat_Km))


df %>% 
  select(mutant, mean_kcat_Km, mutant_type) %>% 
  unique() %>% 
  ggplot(aes(x=mutant, y=mean_kcat_Km, fill=mutant_type)) +
  geom_bar(stat='identity') +
  geom_point(aes(mutant, individual_kcat_Km), data=df, show.legend=F) +
  scale_fill_manual(values = c('blue2','red2','gray'), name='Mutant Type') +
  xlab('Gsp1 mutant') + ylab('kcat/Km of GAP-activated\nGTP hydrolysis, relative to WT') +
  ylim(c(0,80)) +
  theme_custom +
  theme(
    legend.key.size = unit(10,'pt'),
    legend.position = 'bottom'
  )

ggsave('figures/FigS7/FigS7C_GAP.pdf', height=2, width=3)
