source('scripts/config_workspace.R')

library(ggrepel)

df_fit <- 
  read_csv('data/Gsp1_fitness_scores.csv', col_types = cols()) %>% 
  select('sample'=mutant, score)

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
  select(sample, nucleotide, type, mean_v0_norm, sd_norm) %>% 
  unique() %>% 
  pivot_wider(
    id_cols=c(sample, type),
    names_from=nucleotide,
    values_from=c(mean_v0_norm, sd_norm)) %>% 
  mutate(
    ratio = `mean_v0_norm_mant-dGTP`/`mean_v0_norm_mant-dGDP`,
    sd_ratio = compute_ratio_sd(
      `mean_v0_norm_mant-dGTP`,
      `mean_v0_norm_mant-dGDP`,
      `sd_norm_mant-dGTP`,
      `sd_norm_mant-dGDP`
    )  
  ) %>% 
  select(sample, type, r=ratio, sd=sd_ratio) %>% 
  left_join(df_fit, by='sample') %>% 
  mutate(score = ifelse(is.na(score), 0, score)) %>% 
  mutate(position = ifelse(
    sample == 'WT',
    NA,
    str_sub(sample, 2, nchar(sample)-1)))

write_csv(df_to_plot, 'data/GEF_assay_nucleotide_preferences.csv')

df_to_plot %>% 
  filter((grepl('F', sample)) | (sample=='WT')) %>% 
  ggplot(aes(x=score, y=r, color = type, group=position)) +
    geom_line(color='black', linetype='dashed', lwd=0.25) +
    geom_errorbar(aes(ymin=r-sd, ymax=r+sd), width=0.25, alpha=1, show.legend = F) +
    geom_point(shape=16, size = 2.5, alpha=1, show.legend = T) +
    geom_text_repel(aes(label=sample), size=2.5, show.legend = F) +
    xlab('Fitness score') + ylab('Relative preference for GTP\n(ratio of v0, GTP/GDP)') + 
    scale_color_manual(values=c('blue','red','black'), name='Mutant Type') +
    theme_custom +
    theme(legend.key.size = unit(5, 'points'))

ggsave('figures/Fig3/Fig3E_pointplot.pdf', height=2.9, width=4.3)




