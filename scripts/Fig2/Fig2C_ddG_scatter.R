# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# load fitness data
df_fit <- read_csv('data/Gsp1_fitness_scores.csv', col_types = cols())

# load ddg data
df_seq <-
    read_csv('data/Gsp1_Ran_sequence_alignment.csv', col_types=cols()) %>% 
    mutate(pos_Sc = as.integer(pos_Sc), pos_Hs = as.integer(pos_Hs))

df_ddg_3m1i <-
    read_csv('data/ddg/20200806_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Sc', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(struct='3m1i', species='Sc') %>% 
    right_join(df_seq, by='pos_Sc') %>% 
    rename('position'='pos_Sc') %>% 
    filter(position != 71)  # don't use this position, as the structure is mutated Q71L


# burial, based on SASA across all structures
df_bur <- read_csv('data/burial.csv', col_types=cols())  %>% 
    select(position, sasa_group_category) %>% 
    mutate(sasa_group_category = ifelse(sasa_group_category=='mixed', 'surface',sasa_group_category))


# Scatterplots by burial
df_tmp <-
    left_join(df_ddg_3m1i, df_fit, by=c('mutation'='mutant','aa_from','position','aa_to')) %>%
    left_join(df_bur, by='position') %>% 
    filter(!low_reads_flag) %>% 
    select(mutation, aa_from, position, aa_to, ddg, score, bin, sasa_group_category)

df <-
    bind_rows(df_tmp, mutate(df_tmp, sasa_group_category='all positions')) %>% 
    mutate(sasa_group_category = factor(sasa_group_category,
              levels=c('all positions','surface','interface core','structure core')))


# compute correlations with and without toxic mutations
fit_trend_fitness_ddg <- function(df) {
    df %>% 
    select(mutation, ddg, score, sasa_group_category) %>% 
    group_by(sasa_group_category) %>% 
    mutate(corr = cor(score, ddg)) %>% 
    nest(data = c(mutation, ddg, score)) %>% 
    mutate(fit = map(data, ~ lm(ddg ~ score, data = .)),
           tidied = map(fit, broom::tidy)) %>% 
    unnest(cols=c(tidied)) %>% 
    select(sasa_group_category, corr, term, estimate) %>% 
    pivot_wider(names_from=term, values_from=estimate) %>% 
    rename('slope'=score, 'intercept'=`(Intercept)`) %>% 
    mutate(corr=round(corr, 3), slope=round(slope, 3), intercept = round(intercept, 2))
}

linear_fits <- left_join(
    df %>% 
        fit_trend_fitness_ddg() %>% 
        rename('corr_with'='corr','intercept_with'='intercept','slope_with'='slope'),
    df %>% 
        filter(bin != 'toxic') %>% 
        fit_trend_fitness_ddg() %>% 
        rename('corr_wout'='corr','intercept_wout'='intercept','slope_wout'='slope'),
    by = 'sasa_group_category'
) 

linear_fits

df %>% 
    ggplot(aes(x=score, y=ddg, fill=bin)) +
    geom_point(color='black', shape=21, stroke=0.1, size=0.7) +
    facet_grid(~sasa_group_category, drop = FALSE) +
    geom_smooth(data=filter(df, bin!='toxic'),
                aes(group=1), method="lm", se=FALSE,
                formula='y~x', color='chartreuse2', size=0.7,
                show.legend = FALSE) +
    geom_smooth(data=df,
                aes(group=1), method="lm", se=FALSE,
                formula='y~x', color='dodgerblue1', size=0.5,
                show.legend = FALSE) +    
    scale_fill_manual(values = BIN_COLORS, breaks=BIN_ORDERING) +
    xlab('Fitness') + ylab('calculated ddG of mutation (REU)') +
    theme_custom +
    theme(legend.position='none')
    
ggsave('figures/Fig2/Fig2D_ddG_scatter.pdf', width=6.5, height = 2)

