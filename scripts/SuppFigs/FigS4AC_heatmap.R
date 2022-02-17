# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# load additional package
library(ComplexHeatmap)
library(ggrepel)
library(scales)

AA_ORDERING <- c('W','Y','F','M','L','I','V','A','C','G',
                'P','T','S','Q','N','E','D','H','R','K')

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
    rename('position'='pos_Sc')

df <- 
    left_join(df_ddg_3m1i, df_fit, by=c('mutation'='mutant','aa_from','position','aa_to')) %>% 
    mutate(bin = factor(bin,  levels=BIN_ORDERING)) %>% 
    filter(!low_reads_flag) %>% 
    filter(! (aa_from=='L' & position == 71)) %>%  # mutant in the structure
    filter(! (('*' %in% aa_from) | ('*' %in% aa_to)) ) # no truncations in ddG dataset

## What positions does Rosetta think have ddG > 3, but actually are WT-like

mat <-
    df %>% 
    mutate(aa_to = factor(aa_to, levels=AA_ORDERING),
           aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
    filter(bin=='WT-like', ddg > 3) %>% 
    group_by(aa_from, aa_to) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    select(aa_from, aa_to, n) %>% 
    unique() %>% 
    complete(aa_from, aa_to, fill=list('n'=0)) %>% 
    pivot_wider(names_from='aa_from', values_from='n') %>%
    column_to_rownames('aa_to') %>% 
    as.matrix

# get total number of counts for each aa_from

counts <- 
    df_fit %>% 
    mutate(aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
    group_by(aa_from) %>% 
    summarise(n=n()) %>% 
    filter(aa_from != '*') %>% 
    pull(n)


# WT-like mutations predicted as\ndisruptive by Rosetta (ddG > 3)

ha = HeatmapAnnotation(
	barplot = anno_barplot(1:21, axis = TRUE, 
		
        height = unit(1.5, "cm")
))

hm <- 
    Heatmap(
        mat[rev(AA_ORDERING),AA_ORDERING],
        
        rect_gp = gpar(col = "black", lwd = 0.2),
        name='Number of outliers', 
        show_heatmap_legend = T,
        col = colorRamp2(seq(0, 8, by=1), rev(RColorBrewer::brewer.pal(9, 'Spectral'))),
        width = unit(4, 'cm'), height = unit(4, 'cm'),
         
        heatmap_legend_param = list(
            title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
            labels_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
            legend_height = unit(3, 'cm'),
            grid_width = unit(3, 'mm')),

        top_annotation = HeatmapAnnotation(
            'number of mutants\nin the dataset' = anno_barplot(
                counts,
                ylim = c(0, 400),
                axis_param = list(
                    gp=gpar(fontsize=6),
                    at = c(0, 200, 400), 
                    labels = c(0, 200, 400)
            )
            ),
            height = unit(5, 'mm'),
            annotation_name_gp= gpar(fontsize = 7)
        ),

        row_title = 'Substituted Amino Acid', cluster_rows = F, row_names_side = 'left',
        row_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
        row_names_gp = gpar(fontsize = 7, fontfamily='Courier'),

        column_title = 'Original Amino Acid', column_title_side = 'bottom',
        cluster_columns = F, column_names_side = 'bottom',
        column_names_rot = 0, column_names_centered = TRUE,
        column_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
        column_names_gp = gpar(fontsize = 7, fontfamily='Courier')
  )

pdf('figures/FigS4/FigS4A_ddG_outliers_heatmap.pdf', height = 3, width = 3)
draw(hm)
dev.off()

hydrophobics <- c('W','Y','F','M','L','I','V','A')

df_to_plot <-
    df %>% 
    mutate(aa_to = factor(aa_to, levels=AA_ORDERING),
           aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
    filter(bin=='WT-like', ddg > 3) %>% 
    group_by(position) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    select(aa_from, position, n) %>% 
    unique() %>% 
    right_join(unique(select(df,aa_from, position)), by=c('aa_from','position')) %>% 
    filter(!is.na(aa_from)) %>% 
    arrange(position)  %>% 
    mutate(n = ifelse(is.na(n), 0, n)) %>% 
    mutate(label=ifelse(n>1, paste0(aa_from,position), NA)) %>% 
    mutate(hydrophobic=ifelse(aa_from %in% hydrophobics,T,F))

df_to_plot %>% 
    ggplot(aes(x=position, y=n)) +
    geom_vline(xintercept = 95.5, lty='dotted') +
    geom_point(aes(color=hydrophobic), show.legend = F) +
    geom_text_repel(aes(label=label), size=2) +
    scale_color_manual(values=c('black','red')) +
    scale_y_continuous(breaks= pretty_breaks()) +
    xlab('Gsp1 sequence position') + ylab('Number of WT-like mutations predicted as\ndisruptive by Rosetta (ddG > 3)') +
    theme_custom

ggsave('figures/FigS4/FigS4C_ddG_outliers_position.pdf', width=4, height = 3)

# Show mis-represented hydrophobics in c-terminal lobe in structure
df_to_plot %>% 
    filter(position > 87 & !is.na(label) & hydrophobic) %>% 
    pull(position) %>% 
    paste(collapse='+')


## Are these positions also enriched as wrong in the benchmark dataset?
df_bmrk <-
    read_csv('data/ddg/ddg_benchmark_dataset.csv',col_types=cols()) %>% 
    extract(col=mutation_short, into=c('aa_from','position', 'aa_to'), regex='([A-Z])([0-9]{1,4})([A-Z])', convert=T, remove=F) %>% 
    mutate(ddg_calc = ddg_calc*0.298) %>%  # scaling factor from benchmarking
    select(record_ID, aa_from, aa_to,ddg_expt_protherm, ddg_expt_kellogg, ddg_calc)


df_bmrk %>% 
    filter(!is.na(ddg_expt_kellogg) & !is.na(ddg_expt_protherm))


# count the number of records we consider for the benchmark heatmap: 1106
df_bmrk %>% 
    select(record_ID, aa_from, aa_to, ddg_expt_kellogg, ddg_calc) %>% 
    filter(!is.na(ddg_expt_kellogg)) %>% 
    mutate(difference = abs(ddg_expt_kellogg-ddg_calc)) %>% 
    mutate(aa_to = factor(aa_to, levels=AA_ORDERING),
           aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
    nrow()

# count the number of aa_from=hydrophobic: 537
df_bmrk %>% 
    select(record_ID, aa_from, aa_to, ddg_expt_kellogg, ddg_calc) %>% 
    filter(!is.na(ddg_expt_kellogg)) %>% 
    mutate(difference = abs(ddg_expt_kellogg-ddg_calc)) %>% 
    mutate(aa_to = factor(aa_to, levels=AA_ORDERING),
           aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
    filter(aa_from %in% hydrophobics) %>% 
    nrow()

# count the number of aa_from=hydrophobic with difference > 3, 2, or 1
df_bmrk_dif_hydrophobic <-
    df_bmrk %>% 
    select(record_ID, aa_from, aa_to, ddg_expt_kellogg, ddg_calc) %>% 
    filter(!is.na(ddg_expt_kellogg)) %>% 
    mutate(difference = abs(ddg_expt_kellogg-ddg_calc)) %>% 
    mutate(aa_to = factor(aa_to, levels=AA_ORDERING),
           aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
    filter(aa_from %in% hydrophobics)

filter(df_bmrk_dif_hydrophobic, difference > 3) %>% nrow  # 16
filter(df_bmrk_dif_hydrophobic, difference > 2) %>% nrow  # 40
filter(df_bmrk_dif_hydrophobic, difference > 1) %>% nrow  # 158

# df_bmrk %>% 
#     select(record_ID, aa_from, aa_to, ddg_expt_kellogg, ddg_calc) %>% 
#     filter(!is.na(ddg_expt_kellogg)) %>% 
#     mutate(difference = abs(ddg_expt_kellogg-ddg_calc)) %>% 
#     mutate(aa_to = factor(aa_to, levels=AA_ORDERING),
#            aa_from = factor(aa_from, levels=AA_ORDERING)) %>% 
#     group_by(aa_from, aa_to) %>% 
#     mutate(n_tot = n()) %>% 
#     filter(difference > 3) %>% 
#     group_by(aa_from, aa_to) %>% 
#     mutate(n_diff = n()) %>% 
#     ungroup() %>% 
#     complete(aa_from, aa_to, fill=list('n_diff'=0, 'n_tot'=1)) %>% 
#     mutate(x= n_diff/n_tot) %>% 
#     ggplot(aes(aa_from, aa_to, fill=n_diff/n_tot)) + 
#     geom_tile(color='black') +
#     xlab('Original Amino Acid') + ylab('Substituted Amino Acid') +
#     ggtitle('Fraction of mutations poorly predicted by\nRosetta in Benchmark (ddG off by > 3)') +
#     scale_fill_distiller(name='Fraction of\nmutations',palette = "Spectral") +
#     guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5)) +
#     theme_custom +
#     theme(axis.text.x= element_text(family='Courier'),
#           axis.text.y= element_text(family='Courier')) 

# ggsave('figures/FigS4/FigS4_ddG_outliers_heatmap_bmark.pdf', width=3, height = 2.65)



# How many hydrophobics in each lobe
df_seq %>% 
    select(aa_Sc, pos_Sc) %>% 
    filter(!is.na(pos_Sc)) %>% 
    mutate('lobe'=case_when(
        (pos_Sc >= 1) & (pos_Sc <= 95) ~ 'effector',
        (pos_Sc >= 96) & (pos_Sc <= 171) ~ 'allosteric',
        T ~ 'C-terminal extension')
    ) %>% 
    # filter(aa_Sc %in% c('I','L','V')) %>% 
    filter(aa_Sc %in% hydrophobics) %>% 
    group_by(lobe) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n))

