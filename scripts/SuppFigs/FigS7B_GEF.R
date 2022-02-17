source('scripts/config_workspace.R')

library(ggrepel)

sample_order <- c('WT','F28V','F28Y','H50N','H50R','F54A','F54W',
                  'N156A','N156W','F159L','F159W','F163L','F163Y')

df <- read_csv('data/GEF_assay_data_and_fits.csv', show_col_types = FALSE) %>% 
  mutate(
      date=factor(date),
      plate=as.numeric(plate),
      sample=factor(sample, levels=sample_order),
      nucleotide=factor(nucleotide, levels=c('mant-dGTP','mant-dGDP')),
      replicate=factor(replicate),
      GEF_conc=factor(GEF_conc)
      ) %>% 
      filter(!grepl('G35', id))

df_slopes <- 
  df %>% 
  filter(replicate!=0) %>%  # remove no enzyme wells
  select(id, date, sample, nucleotide, replicate, slope) %>% 
  unique %>% 
  mutate(name = case_when(
    ((sample == 'WT') & (date=='20211213')) ~ 'WT_1',
    ((sample == 'WT') & (date=='20220118')) ~ 'WT_2',
    TRUE ~ as.character(sample)
  )) %>% 
  mutate(name = factor(name, levels=c('WT_1','WT_2',sample_order))) %>% 
  mutate(type=case_when(
      sample %in% c('F28V','H50N','F54A','N156A','F159L','F163L') ~ 'Toxic/GOF',
      sample %in% c('F28Y','H50R','F54W','N156W','F159W','F163Y') ~ 'Tolerant',
      sample == 'WT' ~ 'WT')
  )


compute_ratio_sd <- function(a,b,sd_a,sd_b) {

  # function for the formula to propagate standards
  # deviations across the division function
  # f = a/b
  # sd_f = |f| * sqrt( (sd_a/a)^2 + (sd_b/b)^2 )

  abs(a/b) * sqrt( (sd_a/a)^2 + (sd_b/b)^2 )

}

# Compute v0 values normalized to WT
df_WT_v0_sd <-
  df_slopes %>% 
  filter(sample=='WT') %>% 
  group_by(name, nucleotide) %>% 
  mutate(
    WT_v0 = -mean(slope),
    WT_sd = sd(slope)
  ) %>% 
  ungroup() %>% 
  select(date, nucleotide, WT_v0, WT_sd) %>% 
  unique

df_to_plot <-
  df_slopes %>% 
  left_join(df_WT_v0_sd, by=c('date','nucleotide')) %>% 
  mutate(v0_norm = -slope/WT_v0) %>% 
  group_by(sample, nucleotide) %>% 
  mutate(
    mean_v0_norm = mean(v0_norm),
    sd_norm = sd(v0_norm)) %>% 
  ungroup() %>% 
  select(sample, nucleotide, v0_norm, mean_v0_norm, sd_norm, WT_v0, type)

df_to_plot %>% 
  select(sample, nucleotide, type, mean_v0_norm, sd_norm) %>% 
  unique() %>% 
  ggplot(aes(x=sample, y=mean_v0_norm)) + 
  geom_bar(aes(fill=type), stat='identity') +
  geom_errorbar(aes(ymin=mean_v0_norm-sd_norm,ymax=mean_v0_norm+sd_norm), width=0.3, alpha=0.7) +
  geom_point(data=df_to_plot, aes(x=sample, y=v0_norm), size=0.3, alpha=0.5) +
  facet_grid(
    rows=vars(nucleotide),
    labeller = as_labeller(c(
      'mant-dGTP' = 'Exchange to mant-dGTP',
      'mant-dGDP' = 'Exchange to mant-dGDP'
      ))
    ) +
  xlab('Gsp1 mutant') + ylab('Inital rate (absolute value of v0),\nrelative to WT') +
  scale_fill_manual(values=c('blue2','red2','gray'), name='Mutant Type') +
  theme_custom +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.key.size = unit(10,'pt'))

ggsave('figures/FigS7/FigS7B_relative_rates_bar.pdf', height=3, width=4.2)
