# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

df <- 
    read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    filter(!low_reads_flag) %>% 
    mutate(is_toxic = ifelse((aa_to != '*') & (bin == 'toxic'), T, F))


# STOP codon scores for main text
filter(df, aa_to=='*' & position < 175) %>% summarise(min(score), max(score))
filter(df, aa_to=='*' & position >= 175) %>% summarise(min(score), max(score))

# count the frequency of positions with a given number of toxics
counts_true <- df %>% 
    group_by(position, is_toxic) %>%
    summarise(n = n(), .groups='drop') %>% 
    complete(position, is_toxic, fill=list(n=0)) %>% 
    filter(is_toxic) %>% 
    select(n) %>% 
    table() %>% 
    as.data.frame() %>% 
    rename('n_tox' = '.')

counts_true

# compute an independent distribution for count frequencies
# by randomly assigning the toxics to other positions

# n. of toxic mutations is 1188
m <- df %>% 
    summarise(n = sum(is_toxic)) %>% 
    pull(n)

# total number of mutations is 219*21 = 4599
# but we have some positions with low counts before selection removed
# so we have 4519
t <- nrow(df)

# n. of tolerant = 4519-1188 = 3331
n <- t-m

# we're computing the expectation of pulling x toxics when pulling 21 values
# for a given position, with x from 0-21
k <- 21
x <- 0:k

# expected number of positions when independently drawing
counts_null <- data.frame(list(n_tox=x, Freq=219*dhyper(x, m, n, k))) %>% 
    mutate(n_tox = factor(n_tox))


counts_true %>% 
    left_join(counts_null, by='n_tox') %>% 
    rename('obs'='Freq.x', 'ind'='Freq.y') %>%
    mutate(n_tox = as.integer(n_tox)) %>% 
    pivot_longer(cols = c('obs', 'ind'), names_to = 'dist', values_to='freq') %>% 
    ggplot(aes(x=n_tox, y=freq, color=dist)) +
    geom_point(stroke = 0.5, size = 0.5) + 
    scale_color_manual(name = 'Distribution',
                       values = c('ind' = 'black', 'obs' = 'darkgray'),
                       labels=c('Independent Model','Experimental data')) +
    xlab('Number of toxic mutations at a position') + ylab('Observations') +
    guides(color = guide_legend(nrow=2)) +
    theme_custom +
    theme(legend.position = "bottom")

ggsave('figures/FigS1/FigS1_independent_dist.pdf', width = 3.6, height = 2.5)


# Select toxic positions as 10+ toxic mutations, write out to file
df %>%
    filter(is_toxic) %>%
    group_by(position) %>%
    summarise('n' = n(), .groups='keep') %>%
    filter(n >= 10) %>%
    arrange(position) %>%
    write_csv('data/toxic_positions.csv')
