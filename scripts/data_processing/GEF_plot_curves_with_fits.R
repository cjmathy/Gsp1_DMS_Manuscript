source('scripts/config_workspace.R')

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
      )

# plot all points, by sample and plate
facet_plot <- function(data, sample, plate, date) {
  ggplot(data, aes(x=time, y=fluorescence)) +
  facet_wrap(facets=vars(nucleotide, replicate)) +
  geom_point(aes(color=GEF_conc), size=1) +
  geom_vline(aes(xintercept=t1), linetype='dashed', alpha=0.5) +
  geom_vline(aes(xintercept=t2), linetype='dashed', alpha=0.5) +
  geom_abline(aes(slope=slope,intercept=intercept)) +
  geom_abline(aes(slope=slope_with_bg,intercept=intercept), color='red') +
  ggtitle(paste0('Mutant: ', sample, '     Date: ', date, '     Plate: ', plate)) +
  theme_custom_drafting +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
}

df_plots <-
  df %>%
  mutate(slope_with_bg = ifelse(fittype=='full', slope+slope_no_enzyme, NA)) %>% 
  group_by(sample, date, plate) %>% 
  nest() %>% 
  mutate(plots = pmap(.l=list(data,sample,plate,date), facet_plot)) %>% 
  arrange(sample)

pdf('figures/GEF/GEF_individual_plots_with_fits.pdf',width=10,height=8)
invisible(print(df_plots$plots))
dev.off()
