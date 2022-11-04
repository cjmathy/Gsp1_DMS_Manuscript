source('scripts/config_workspace.R')

# this code is identical to raw_counts_to_scores.R up through the making of df_binned
# then it plots the counts, stratified by bin

# read in counts table
df_raw <- read_csv('data/raw_counts.csv', col_types = cols()) %>% 
    filter(position!=220)

read_avg <- mean(df_raw$counts_0gen, na.rm=T)
read_threshold <- floor(0.02 * read_avg)

df_raw_labeled <-
    df_raw %>% 
    mutate(low_reads_flag = ifelse((counts_0gen <= read_threshold) | is.na(counts_0gen), T, F))

# get the total number of reads in the library
sum_0gen <- sum(df_raw$counts_0gen, na.rm=T)
sum_6gen <- sum(df_raw$counts_6gen, na.rm=T)

# add 0.5 to the allele counts so that there are no log-transformed zeros
df_uncentered <- 
    df_raw_labeled %>% 
    mutate('counts_0gen_adj' = counts_0gen + 0.5,
           'counts_6gen_adj' = counts_6gen + 0.5) %>% 
    mutate(freq_0gen = counts_0gen_adj/sum_0gen,
           freq_6gen = counts_6gen_adj/sum_6gen) %>% 
    mutate(score = log(freq_6gen, base=2)-log(freq_0gen, base=2))

# compute the average score of a WT synonomous codon allele
WT_syn_mean_score <- 
    df_uncentered %>% 
    filter(aa_from==aa_to) %>% 
    filter(!low_reads_flag) %>% 
    pull(score) %>% 
    mean(na.rm=T)

# adjust the scores so that the WT synonomous scores are centered at 0
# since these are log transformed values, this is akin to normalization
# by dividing by the mean of the WT distribution
df <- mutate(df_uncentered, score = score-WT_syn_mean_score)

WT_mean <- filter(df, aa_to==aa_from, !low_reads_flag) %>% pull(score) %>% mean(na.rm=T)
WT_sd <- filter(df, aa_to==aa_from, !low_reads_flag) %>% pull(score) %>% sd(na.rm=T) 
STOP_mean <- filter(df, aa_to=='*', position < 175, !low_reads_flag) %>% pull(score) %>% mean(na.rm=T)
STOP_sd <- filter(df, aa_to=='*', position < 175, !low_reads_flag) %>% pull(score) %>% sd(na.rm=T) 

Z <- 2

df_binned <-
    df %>% 
    mutate('Z_WT' = (score - WT_mean) / WT_sd,
           'Z_STOP' = (score - STOP_mean) / STOP_sd) %>% 
    mutate(bin=case_when(
                low_reads_flag ~ 'low-reads',
                is.na(score) ~ 'drop-outs',
                (Z_WT >= Z) ~ 'beneficial',
                ((Z_WT < Z) & (Z_WT > -Z)) ~ 'WT-like',
                ((Z_WT <= -Z) & (Z_STOP >= Z)) ~ 'intermediate',
                ((Z_STOP < Z) & (Z_STOP > -3)) ~ 'STOP-like',
                (Z_STOP <= -3) ~ 'toxic'
          )
    ) 



df_binned <- df_binned %>%
    mutate(bin=factor(bin, levels=c(
        'toxic',
        'STOP-like',
        'intermediate',
        'WT-like',
        'beneficial'
    )))


BIN_ORDERING <- c('toxic','STOP-like','intermediate','WT-like','beneficial')
BIN_COLORS <- setNames(as.list(RColorBrewer::brewer.pal(7, 'RdBu')[1:5]),BIN_ORDERING)

df_binned %>% 
    filter(!low_reads_flag) %>% 
    filter(counts_0gen < 20000) %>% 
    ggplot(aes(counts_0gen, fill=bin)) +
    geom_bar(stat='bin',color='black',size=0.2, bins=100) +
    facet_grid(rows='bin', scales='free') +
    scale_fill_manual(values = BIN_COLORS, breaks=BIN_ORDERING, name='Mutant Category') +
    xlab('Read counts in initial library') + ylab('Allele count') +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_custom +
    theme(legend.key.size = unit(0.5,"line"),
          legend.position="bottom",
          panel.border = element_rect(colour = "black", fill=NA, size=0.25))

df_binned %>% filter(!low_reads_flag) %>% nrow

# 22 alleles had 20000+ reads, none were beneficial
df_binned %>% 
    filter(!low_reads_flag) %>% 
    filter(counts_0gen >= 20000) %>% 
    print(n=22)

ggsave('figures/FigS2/FigS2_counts_hist.pdf', height=4, width=5)

