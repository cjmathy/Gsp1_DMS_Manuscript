source('scripts/config_workspace.R')

library(glue)
library(ROCR)

df_aln <- read_delim('data/Ras_data/HRas_Gsp1_index.txt', delim='\t', show_col_types=F)

df_Gsp1 <-
    read_csv('data/Gsp1_fitness_scores.csv', show_col_types=F) %>% 
    left_join(select(df_aln, aln_num, 'position'=Gsp1_res_num),by='position') %>% 
    select(aln_num,'aa_from_Gsp1'=aa_from, 'position_Gsp1'='position', aa_to,
          'score_Gsp1'=score, 'bin_Gsp1'=bin,'low_reads_flag_Gsp1'=low_reads_flag)

df_HRas <-
    read_csv('data/Ras_data/HRas188_BaF3_processed.csv', show_col_types=F) %>% 
    left_join(select(df_aln, aln_num, 'position'=HRas_res_num),by='position') %>% 
    select(aln_num, 'aa_from_HRas'=aa_from, 'position_HRas'=position, aa_to,
           'score_HRas'=score)

activating_threshold <- 1.5*sd(pull(df_HRas, score_HRas))
activating_threshold

# Positions unique to one protein or the other do not appear in df
# i.e aln_num = 105,106,133 in HRas are insertions compared to Gsp1
# and aln_num = 38 in Gsp1 is not in HRas
# none of the HRas positions 105, 106, 133 have activating mutations
# also the termini that are not included in both experiments, i.e. aln_num = 1
# gets removed because it is present in the Gsp1 dataset but not in the HRas one

df <- 
    inner_join(df_Gsp1, df_HRas, by=c('aln_num', 'aa_to')) %>% 
    filter(!low_reads_flag_Gsp1) %>% 
    select(-low_reads_flag_Gsp1)


## ROC analysis

# From the preprint (Hidalgo et al 2022 bioRxiv) results:

# To generate a ROC curve, the fitness data for individual mutations in a 
# particular bacterial dataset (e.g., H-Ras2-166+GAP) are used to predict 
# the fitness of mutations in the Ba/F3 experiment. A variable threshold value 
# of fitness is used, and for each threshold value, mutations in the bacterial 
# dataset with a fitness value greater than that threshold are considered to 
# predict activation in Ba/F3 data. Mutations with a fitness score greater than 
# 1.5 times the standard deviation in the Ba/F3 dataset are considered activating 
# (i.e., true positives). For each of the bacterial datasets, an ROC curve is generated
# by graphing the fraction of true positives versus the fraction of false positives at
# various threshold settings, and the estimated area under the curve (AUC) gives a
# measure of the overall prediction accuracy.

# From the preprint (Hidalgo et al 2022 bioRxiv) methods:
# When using a second saturation-mutagenesis dataset as the reference, 
# mutations with a score of >1.5 times the standard deviation are considered
# true positives, and  the rest as true 6negatives. The ROC curve is generated
# by graphing the fraction of true positives versus the fraction of false 
# positives at various threshold settings, and the estimated area under the 
# curve (AUC) gives a measure of the overall accuracy of the predictions.

df_ROC <- mutate(df, 'TP'=ifelse(score_HRas > activating_threshold, T, F))    
pred <- prediction(-1*df_ROC$score_Gsp1, df_ROC$TP)
perf <- performance(pred,"tpr","fpr")

auc.tmp <- performance(pred,"auc")
auc <- as.numeric(auc.tmp@y.values)
auc

# randomize Gsp1 scores
set.seed(2022)
df_ROC_rand <-
    df_ROC %>% 
    transform('score_Gsp1'=sample(score_Gsp1)) %>% 
    tibble()

pred_rand <- prediction(-1*df_ROC_rand$score_Gsp1, df_ROC_rand$TP)
perf_rand <- performance(pred_rand,"tpr","fpr")

auc_rand.tmp <- performance(pred_rand,"auc")
auc_rand <- as.numeric(auc_rand.tmp@y.values)
auc_rand

pdf('figures/Fig4/Fig4C_HRas_ROC_new.pdf', family = 'Helvetica', width=2.1,height=2.3,pointsize=6)
plot(perf, col='black', lwd=2,
     main='ROC curve: predictions of\nHRas GOF mutations\nusing Gsp1 fitness scores', 
     xlab='False Positive Rate',ylab='True Positive Rate'
)
auc_round <- round(auc, 3)
plot(perf_rand, add=TRUE, col='blue', lwd=2)
auc_rand_round <- round(auc_rand, 3)
abline(a=0,b=1, lty='dotted')
legend(x=0.15, y=0.2,
      legend=c(glue('AUC: {auc_round}'),
               glue('AUC (randomized): {auc_rand_round}')
      ),
      fill=c('black','blue'),
      bty='n'
)
dev.off()


# easy visualization of the predictive power of Gsp1 score in
# identifying HRas activating mutations (TPs = True Positives)
df_ROC %>% 
    ggplot(aes(x=TP, y=score_Gsp1)) +
    geom_violin() +
    theme_custom





#### Import SCA data to make a table of annotated positions

# get positions of HRas sector from https://github.com/ranganathanlab/pySCA/blob/master/notebooks/SCA_G.ipynb

sector1 <- c(
    5, 10, 11, 14, 15, 16, 22, 28, 32, 34, 35, 36, 39, 42, 54, 56, 57, 58,
    59, 60, 61, 62, 63, 64, 68, 71, 72, 73, 75, 81, 83, 85, 96, 99, 110,
    116, 117, 119, 144, 145, 146, 147, 156
)
sector2 <- c(17, 23, 82, 84, 90, 115, 123, 125, 129, 130, 134, 141, 143)

sca_HRas <- c(sector1, sector2)

# Prepare a table for Fig 4D
toxics <- read_csv('data/toxic_positions.csv', col_types=cols()) %>% 
    filter(position != 220) %>% 
    pull(position)

sca_gsp1 <-
    df_aln %>% 
    filter(HRas_res_num %in% sca_HRas) %>% 
    pull(Gsp1_res_num)

gofs_all <-
    df_ROC %>% 
    group_by(position_Gsp1, TP) %>% 
    summarise(n=n()) %>% 
    filter(TP) %>% 
    pull(position_Gsp1)

gofs_2plus <-
    df_ROC %>% 
    group_by(position_Gsp1, TP) %>% 
    summarise(n=n()) %>% 
    filter(TP) %>% 
    filter(n>1) %>% 
    pull(position_Gsp1)

# Chi-square shows that we should use the 2+ GOF mutations at a position to label it as GOF in HRas
df %>% 
    select('pos'=position_Gsp1) %>%  # only use positions aligned with HRas
    unique() %>% 
    mutate(
        tox = ifelse(pos %in% toxics, 'Toxic', 'Nontoxic'),
        sca = ifelse(pos %in% sca_gsp1, 'T', 'F'),
        gof_all = ifelse(pos %in% gofs_all, 'T', 'F'),
        gof_2plus = ifelse(pos %in% gofs_2plus, 'T', 'F')) %>% 
    pivot_longer(cols=c(sca, gof_all, gof_2plus), names_to='category', values_to='in_category') %>% 
    group_by(tox, category, in_category) %>% 
    count() %>% 
    pivot_wider(names_from=tox, values_from=n) %>% 
    nest(data = c(in_category, Toxic, Nontoxic)) %>% 
    mutate(
        chisquare = map(data,
            function(df) {
                df %>% 
                column_to_rownames('in_category') %>% 
                chisq.test()
            }),
        tidied = map(chisquare, broom::tidy)
    ) %>% 
    unnest(cols=c(tidied)) %>% 
    select(category, p.value)

# Results:
#    - Sectors have a significant association with Toxics
#    - Gain-of-functions have a significant association with Toxics if
#      positions with more than one GOF mutation in the HRas data are
#      considered. If even on GOF mutation labels a position as GOF,
#      there is just barely a significant association (P = 0.0411)

# # A tibble: 3 x 2
# # Groups:   category [3]
#   category    p.value
#   <chr>         <dbl>
# 1 gof_2plus 0.0000711
# 2 gof_all   0.0411   
# 3 sec       0.00445  

df %>% 
    select('pos'=position_Gsp1) %>%  # only use positions aligned with HRas
    unique() %>% 
    mutate(
        tox = ifelse(pos %in% toxics, 'true', 'false'),
        sca = ifelse(pos %in% sca_gsp1, 'true', 'false'),
        gof_all = ifelse(pos %in% gofs_all, 'true', 'false'),
        gof_2plus = ifelse(pos %in% gofs_2plus, 'true', 'false')) %>% 
    write_csv('data/Ras_data/Ras_sca_gof_table.csv')
