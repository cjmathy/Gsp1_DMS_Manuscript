# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

hydrophobics <- c('W','Y','F','M','L','I','V','A')

df <-
    read_csv('data/ddg/ddg_benchmark_dataset.csv', col_types=cols()) %>% 
    mutate(position = as.integer(substring(mutation_short, 2, nchar(mutation_short)-1)),
           aa_from = substring(mutation_short, 1, 1),
           aa_to = substring(mutation_short, nchar(mutation_short), nchar(mutation_short))) %>% 
    mutate(ddg_expt = case_when(
              is.na(ddg_expt_kellogg) ~ ddg_expt_protherm,
              is.na(ddg_expt_protherm) ~ ddg_expt_kellogg,
              T ~ ddg_expt_protherm)) %>% 
    mutate(ddg_calc_adj = ddg_calc*0.298) %>% 
    mutate(is_ILV = ifelse(aa_from %in% c('I','L','V'), T, F)) %>% 
    mutate(is_hydrophobic = ifelse(aa_from %in% hydrophobics, 'Mutations at hydrophobics', 'Mutations at other amino acids'))

label <-
    df %>% 
    select(is_hydrophobic) %>% 
    group_by(is_hydrophobic) %>% 
    summarise(n=n()) %>% 
    mutate(label = paste0('n = ', n))

df %>% 
    ggplot(aes(x=ddg_expt, y=ddg_calc_adj)) + 
    geom_point(alpha=0.3, shape=16, size=1) +
    # geom_point(alpha=0.5, shape=16, size=1, color='black', data=filter(df, !is_ILV)) +
    # geom_point(alpha=0.5, shape=16, size=1, color='red', data=filter(df, is_ILV)) +
    geom_text(data=label, aes(label=label, x=10, y=18), size = 2, hjust = 1) +
    facet_grid(~is_hydrophobic) +
    geom_abline(slope=1,intercept=0) +
    xlab('ddG (expt)') + ylab('ddG (calc)') +
    ggtitle('Accuracy of Rosetta ddG modeling in benchmark set') +
    theme_custom

ggsave('figures/FigS4/FigS4C_scatter.pdf', width=4.5, height = 2.6)

df %>% 
    mutate('abs_error'=ddg_expt-ddg_calc_adj) %>% 
    group_by(is_hydrophobic) %>% 
    summarise(mean(abs_error), sd(abs_error))
