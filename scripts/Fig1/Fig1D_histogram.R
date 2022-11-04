# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

df <-
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>% 
    mutate(score = as.double(score)) %>% 
    filter(!low_reads_flag) %>%  
    mutate(bin=case_when(aa_to==aa_from & !is.na(score) ~ 'WT',
                         aa_to=='*'  & !is.na(score) ~ 'STOP',
                         T ~ bin)) %>% 
    mutate(table_cat = bin) %>%
    mutate(table_cat = case_when(
        table_cat != 'STOP' ~ table_cat,
        (table_cat == 'STOP') & (position < 175) ~ 'STOP1',
        (table_cat == 'STOP') & (position >= 175) ~ 'STOP2'
    )) %>% 
    mutate(
        bin=factor(bin, levels=rev(BIN_ORDERING)),
        table_cat = factor(table_cat, levels=c(
            'WT','STOP1', 'STOP2', 'beneficial',
            'WT-like', 'intermediate', 'STOP-like', 'toxic'
        ))) %>% 
    group_by(table_cat) %>% 
    mutate(n=n(),
           mean_score=mean(score, na.rm=T),
           sd = sd(score, na.rm=T)) %>% 
    ungroup()


# get number of low-reads
read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    mutate(score = as.double(score)) %>% 
    filter(low_reads_flag) %>% 
    nrow

# use the counts to edit the legend in the final figure, adding the n for each category
df_means <-
    df %>% 
    select(table_cat, n, mean_score, sd) %>% 
    unique %>% 
    arrange(table_cat)

df_means


bw <- 0.1  # bin width
lt <- 0.05 # line thickness

df %>% 
    ggplot(aes(x=score, fill=bin)) + 
    geom_bar(data=df %>% filter(!bin %in% c('WT','STOP')), 
                   aes(y = ..count..), stat='bin',
                   binwidth=bw,color='black',size=lt) +
    geom_bar(data=df %>% filter(bin=='STOP'), 
                   aes(y = ..count..), stat='bin',
                   binwidth=bw, color='black', size=lt, alpha=0.8) +
    geom_bar(data=df %>% filter(bin=='WT'), 
                   aes(y = ..count..), stat='bin',
                   binwidth=bw, color='black', size=lt, alpha=0.8) +
    scale_fill_manual(values = BIN_COLORS, breaks=BIN_ORDERING) +
    xlab('Fitness Score') + ylab('Count') +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    theme_custom +
    theme(legend.key.size = unit(0.5,"line"),
          legend.position="bottom",
          panel.border = element_rect(colour = "black", fill=NA, size=0.25))

ggsave('figures/Fig1/Fig1D_histogram.pdf', height=2.9, width=3.3)

