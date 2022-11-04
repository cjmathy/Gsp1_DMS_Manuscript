# load tidyverse and circlize, and set theme for plotting
source('scripts/config_workspace.R') 

# load additional package(s)
library(cowplot)

df <-
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    mutate(Category=case_when(aa_to==aa_from & !is.na(score) ~ 'WT',
                              aa_to=='*'  & !is.na(score) ~ 'STOP',
                              T ~ 'Mutant')) %>% 
    mutate(Category = factor(Category, levels = c('Mutant', 'STOP', 'WT'))) %>% 
    filter(!low_reads_flag)

# max/min of STOPs before position 175
filter(df, position < 175, aa_to=='*') %>% summarise(min(score))  # -2.90
filter(df, position < 175, aa_to=='*') %>% summarise(max(score))  # -1.10
filter(df, position < 175, aa_to=='*') %>% summarise(mean(score)) # -2.06

# max/min of STOPs at or after position 175 (inclusive)
filter(df, position >= 175, aa_to=='*') %>% summarise(min(score))  # -10.5
filter(df, position >= 175, aa_to=='*') %>% summarise(max(score))  # 0.452
filter(df, position >= 175, aa_to=='*') %>% summarise(mean(score)) # -6.66


df

p1 <- df %>% 
    ggplot() +
    geom_point(data=df %>% filter(Category == 'Mutant'),
               mapping=aes(x=position, y=score),
               fill='gray', color='transparent',
               size=0.5, shape=21) +
    geom_point(data=df %>% filter(Category == 'WT'),
               mapping=aes(x=position, y=score),
               fill='forestgreen', color='transparent',
               size=0.75, shape=21) +
    geom_point(data=df %>% filter(Category == 'STOP'),
               mapping=aes(x=position, y=score),
               fill='black', color='transparent',
               size=0.75, shape=21) +
    ylim(c(-12,2)) + xlab('Gsp1 sequence position') + ylab('Fitness Score') +
    theme_custom + 
    theme(legend.position = 'none',
          plot.margin = unit(c(5.5, 1.5, 5.5, 5.5), 'points'),
          panel.border = element_rect(colour = "black", fill=NA, size=0.25)
    )

p2 <- df %>% 
    ggplot(aes(x=score, fill=Category)) + 
    geom_density(aes(y=..count..), size=0.2) +
    scale_fill_manual(values = list( 'WT'='forestgreen','STOP'='black','Mutant'='gray')) +
    coord_flip() +
    xlim(c(-12,2)) + ylab('Counts') + xlab('') +
    theme_custom +
    theme(legend.position=c(0.62,0.25),
          legend.key.size = unit(0.4, 'lines'),
          legend.box.background = element_rect(colour = "black"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = unit(c(5.5, 5.5, 5.5, 1.5), 'points'),
          panel.border = element_rect(colour = "black", fill=NA, size=0.25)
          )

plot_grid(p1, p2, rel_widths = c(2.75,1))

ggsave('figures/Fig1/Fig1E_STOPs.pdf', width = 3.3, height = 1.3)
