source('scripts/config_workspace.R')

make_point_plot <- function(df, x, WT_val, remove_y_axis=FALSE) {
    gg <-
        ggplot(df, aes({{x}}, y=mutant_type, color=mutant_type, group=position)) +
        geom_point() +
        geom_line(color='black', linetype='dashed', lwd=0.25) +
        ggrepel::geom_text_repel(
            data=filter(df, mutant_type=='Toxic/GOF'),
            aes(label=mutant), nudge_y=0.15, size=2, show.legend = F, max.overlaps=25) +
        ggrepel::geom_text_repel(
            data=filter(df, mutant_type=='WT-like'),
            aes(label=mutant), nudge_y=-0.15, size=2, show.legend = F, max.overlaps=25) +
        scale_color_manual(name='Type',values=c('blue','red','black')) +
        geom_vline(xintercept=WT_val, linetype='dotted', alpha=0.5) +
        annotate('text', x=WT_val, y=0.5, label='WT', hjust=-0.1, size=2) +
        ylab('Mutant Type') +
        theme_custom +
        theme(legend.position='none')
    if(remove_y_axis) {
        gg <- gg + theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    }
    gg
}

df_fit <- read_csv('data/Gsp1_fitness_scores.csv', col_types = cols()) %>% 
    filter(mutant %in% sample_order) %>% 
    mutate(mutant_type = ifelse(bin=='toxic', 'Toxic/GOF', 'WT-like')) %>% 
    mutate(mutant_type = factor(mutant_type, levels = c('WT-like', 'Toxic/GOF'))) %>% 
    select(mutant, aa_from, position, aa_to, score, mutant_type)

sample_order <- c('WT','F28V','F28Y', 'H50N','H50R','F54A','F54W',
                  'N156A','N156W','F159L','F159W','F163L','F163Y')

df_seq <-
    read_csv('data/Gsp1_Ran_sequence_alignment.csv', col_types=cols()) %>% 
    mutate(pos_Sc = as.integer(pos_Sc), pos_Hs = as.integer(pos_Hs))

# Plot 1: fitness
p1 <-
    make_point_plot(df_fit, x=score, WT_val=0.0) +
    xlim(c(-10,1.5)) +
    xlab('Fitness score')

# Plot 2: melts
df_fit_params <- read_csv('data/CD_melt_fits.csv', col_types=cols())
WT_Tm <- filter(df_fit_params, sample == 'WT GDP')$Tm

p2 <-
    df_fit_params %>% 
    filter(mutant %in% sample_order) %>% 
    filter(mutant!='WT') %>% 
    left_join(select(df_fit, mutant, score), by='mutant') %>% 
    mutate('score' = ifelse(mutant=='WT', 0, score)) %>% 
    make_point_plot(x=Tm, WT_val=WT_Tm, remove_y_axis=TRUE) +
    xlab('Apparent Tm')

# Plot 3: ddG values
p3 <-
    read_csv('data/ddg/3gj0_gsp1_ddg.csv', col_types=cols()) %>% 
    extract(col=mutation, into=c('aa_from','pos_Hs', 'aa_to'), regex='([A-Z])([0-9]{1,3})([A-Z])', convert=T, remove=F) %>% 
    mutate(struct='3gj0', species='Hs') %>% 
    right_join(df_seq, by='pos_Hs') %>% 
    select('position'=pos_Sc, aa_to, ddg) %>% 
    inner_join(df_fit, by=c('position','aa_to')) %>% 
    make_point_plot(x=ddg, WT_val=0.0, remove_y_axis=TRUE) +
    xlab('Calculated ddG of mutation (REU)')


# Plot 4: GEF
p4 <-
    read_csv('data/GEF_assay_nucleotide_preferences.csv', show_col_types=FALSE) %>% 
    filter(sample != 'WT') %>% 
    rename('mutant'=sample, 'mutant_type'=type) %>% 
    make_point_plot(x=r, WT_val=1.0, remove_y_axis=TRUE) +
    xlim(c(0,7)) +
    xlab('Rel. Ratio v0_GTP / v0_GDP')

cowplot::plot_grid(p1, p2, p3, p4, nrow=1, align='v') 

ggsave('figures/FigS7/FigS7A_biophysics.pdf', width = 9, height = 2)

