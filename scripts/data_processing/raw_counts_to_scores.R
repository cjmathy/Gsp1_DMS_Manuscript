# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R') 

# read in counts table
df_raw <- read_csv('data/raw_counts.csv', col_types = cols()) %>% 
    filter(position!=220)

# Set a threshold for low read counts at 2% of the average variant's reads.
# Alleles below this are considered not well represented in the library, and
# are ignored in downstream analyses
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

### binning

WT_mean <- filter(df, aa_to==aa_from, !low_reads_flag) %>% pull(score) %>% mean(na.rm=T)
WT_sd <- filter(df, aa_to==aa_from, !low_reads_flag) %>% pull(score) %>% sd(na.rm=T) 
STOP_mean <- filter(df, aa_to=='*', position < 175, !low_reads_flag) %>% pull(score) %>% mean(na.rm=T)
STOP_sd <- filter(df, aa_to=='*', position < 175, !low_reads_flag) %>% pull(score) %>% sd(na.rm=T) 


# bin scores:
#   - beneficial: better than two SD above WT
#   - WT-like: within two SD (of WT dist.) of WT mean
#   - STOP-like: within two SD (of STOP dist.) above or 3 below STOP mean
#   - intermediate: in between WT-like and STOP-like
#   - toxic: worse than 3 SD from STOP mean
#   - drop-outs: NAs in 6gen
#   - low-reads: present in 6gen, but less than or equal to 28 reads in 0gen
Z <- 2

# could alternately use Z-score for a certain percent of distribution 1-p
# p <- 0.01
# Z <- qnorm(1-p/2)


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


# note that for this dataset, the 6gen dataset, all drop-outs are
# due to low reads, so nothing gets binned as "drop-outs"

# save the dataframe
df_binned %>% 
    select(mutant, aa_from, position, aa_to, score, bin, low_reads_flag) %>% 
    write_csv('/Users/cjmathy/gdrive/gsp1_dms/data/Gsp1_fitness_scores.csv')

