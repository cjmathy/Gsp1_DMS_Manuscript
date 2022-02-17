# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')
library(glue)

# burial, based on SASA across all structures
df_bur <- read_csv('data/burial.csv', col_types=cols()) %>% 
    select(position, sasa_group_category) %>% 
    mutate(sasa_group_category = ifelse(sasa_group_category=='mixed', 'surface',sasa_group_category))

# get residues within 4 angstroms of nucleotide from PDB 3m1iA
contacts_nucleotide <- scan('data/residues_contacting_nucleotide.txt', comment.char='#')

gtpase_regions <- list(
    'P-loop'=seq(19,26),              # Vetter chapter in Wittinghoffer book 2014
    'Extended Switch I'=seq(34,47),   # Vetter et al. Cell 1999
    'Extended Switch II'=seq(67,81),  # Vetter et al. Cell 1999 + Mg coordinating 67 from Scheffzek 1995
    'N/TKxD motif'=seq(124,127),      # Vetter chapter in Wittinghoffer book 2014
    'SAK motif'=seq(152,154)          # Vetter chapter in Wittinghoffer book 2014
)

# define active site as canonical GTPase regions + contacting residues
activesite <- union(contacts_nucleotide, unlist(gtpase_regions)) %>% sort

df <-
    read_csv('~/gdrive/gsp1_dms/data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    filter(!low_reads_flag) %>% 
    mutate(bin=factor(bin, levels=rev(BIN_ORDERING))) %>% 
    left_join(df_bur, by='position') %>% 
    mutate(sasa_group_category = ifelse(position %in% activesite, 'active site', sasa_group_category)) %>% 
    mutate(sasa_group_category = factor(sasa_group_category, levels=c('active site','surface','interface core','structure core')))
    

# Total residues in each group, for figure
totals <-
    df %>% 
    filter(position!=220) %>% 
    filter(aa_to!='*') %>% 
    filter(aa_to!=aa_from) %>% 
    group_by(sasa_group_category) %>% 
    summarise(n=n()) %>% 
    deframe()

df_names <- data.frame(
    'sasa_group_category' = c('active site','surface','interface core','structure core'),
    'name' = c(
        glue('GTPase\nregions\nand\nnucleotide\ncontacts\nn = {n}', n=totals['active site']),
        glue('Surface\nn = {n}', n=totals['surface']),
        glue('Interface\ncore\nn = {n}', n=totals['interface core']),
        glue('Structure\ncore\nn = {n}', n=totals['structure core'])
    )
) %>% mutate(name=factor(name, levels=.$name))

bw <- 0.1  # bin width
lt <- 0.05 # line thickness

df %>% 
    left_join(df_names, by='sasa_group_category') %>% 
    filter(position!=220) %>% 
    filter(aa_to!='*') %>% 
    filter(aa_to!=aa_from) %>% 
    ggplot(aes(x=score, fill=bin)) + 
    geom_bar(aes(y = ..count..), stat='bin',
             binwidth=bw,color='black',size=lt) +
    facet_grid(rows=vars(name)) +
    scale_fill_manual(values = BIN_COLORS, breaks=BIN_ORDERING, name='Mutant\nCategory') +
    xlab('Fitness Score') + ylab('Count') +
    guides(fill = guide_legend(nrow = 3, bycol = TRUE)) +
    theme_custom +
    theme(
        legend.position="none",
        strip.text.y = element_text(angle=0, margin = margin(r = 3, l = 3, unit = "pt"))
    )
    
ggsave('figures/Fig2/Fig2A_histogram_region.pdf', height=3.3, width=2)

df %>% 
    filter(position!=220) %>% 
    filter(aa_to!='*') %>% 
    filter(aa_to!=aa_from) %>% 
    group_by(sasa_group_category, bin) %>% 
    mutate(n=n()) %>% 
    ungroup() %>% 
    group_by(sasa_group_category) %>% 
    mutate(n_region=n()) %>% 
    ungroup() %>% 
    select(sasa_group_category,bin,n,n_region) %>% 
    unique() %>% 
    mutate(percent = n/n_region*100) %>% 
    select(sasa_group_category, bin, percent) %>% 
    pivot_wider(names_from='bin',values_from='percent')  %>% 
    select(sasa_group_category, `WT-like`, `STOP-like`, toxic)


