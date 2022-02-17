# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# load additional packages
library(ComplexHeatmap)

# Write a tab-delimited file for the pymol script that makes Fig4B
read_csv('data/toxic_positions.csv', col_types=cols()) %>% 
write_tsv('scripts/Fig4/toxics.txt', col_names=F)

# get toxics, ordered by severity
toxics <-
    read_csv('data/toxic_positions.csv', col_types=cols()) %>% 
    filter(position != 220) %>% 
    arrange(desc(n)) %>% 
    pull(position)

df_fit <- read_csv('data/Gsp1_fitness_scores.csv', col_types=cols())

df <- 
    df_fit %>% 
    mutate('bin' = factor(bin, levels=BIN_ORDERING),
           'aa_to'=factor(aa_to, levels=AA_ORDERING),
           'aa_from'=factor(aa_from, levels=AA_ORDERING)) %>%
    filter(position %in% toxics)

# prepare a sequence vector for WT Gsp1 so we can box cells in heatmaps that are WT in sequence (i.e. T34T)
gsp1_seq <- select(df_fit, aa_from, position) %>% unique() %>% pull(aa_from) %>% as.character()

# convert dataframes to matrices for each dataset, for plotting
mat <- 
  df %>%
  select(position, aa_to, score) %>%
  pivot_wider(id_cols = aa_to, names_from = position, values_from = score) %>%
  arrange(desc(aa_to)) %>%
  column_to_rownames('aa_to') %>%
  as.matrix()


low_read_mat <-
  df %>%
  select(position, aa_to, low_reads_flag) %>%
  pivot_wider(id_cols = aa_to, names_from = position, values_from = low_reads_flag) %>%
  arrange(desc(aa_to)) %>%
  column_to_rownames('aa_to') %>%
  as.matrix()

mat <- mat[1:21, as.character(toxics)]
low_read_mat <- low_read_mat[1:21, as.character(toxics)]

# function for plotting a heatmap of empiric data
# input dataset should be a matrix, with 21 aa labels as rows
# and each residue position as columns
draw_fitness_heatmap <- function(mat, name, col, label_WT=F, sequence=NULL, gray_out=low_read_mat) {
  Heatmap(
    mat,
    name=name,
    col=col,
    show_heatmap_legend = T,
    na_col='black',

    width = ncol(mat)*unit(1.8, "mm"), 
    height = nrow(mat)*unit(1.8, "mm"),

    row_title = 'Amino Acid Substituted',
    cluster_rows = F,
    row_names_side = 'right',
    row_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
    row_names_gp = gpar(fontsize = 6, fontfamily='Courier'),

    column_title = 'Gsp1 sequence position',
    cluster_columns = F,
    column_names_side = 'top',
    column_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
    column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
    # column_labels = replace(colnames(mat), seq(from = 2, to = length(colnames(mat)), by = 2), ''),

    # borders around all cells, green if synonymous mutation
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (gray_out[i,j]) {
            grid.rect(x, y, width, height, gp = gpar(fill = 'white', lwd=0.5))
            grid.lines(x=c(x-width/2, x+width/2), y=c(y+height/2, y-height/2),
                       default.units = 'npc', gp=gpar(lwd=0.2))
            grid.lines(x=c(x-width/2, x+width/2), y=c(y-height/2, y+height/2),
                       default.units = 'npc', gp=gpar(lwd=0.2))
        }
        if ((rownames(mat)[i] == sequence[j]) & label_WT) {
            grid.circle(x, y, r = min(unit.c(width, height)/2), gp = gpar(fill = 'black'))
        }
        grid.rect(x, y, width, height, gp = gpar(col = 'black', lwd = 0.25, fill = "transparent")) 
      },
          
    # add sequence
    bottom_annotation = HeatmapAnnotation(
        which = 'column',
        'sequence' = anno_text(
            sequence, gp=gpar(fontsize = 6, fontfamily='Courier'),
            rot=0, just = "center", location = 0.5
            )
        ),

    heatmap_legend_param = list(
        title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
        labels_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
        legend_height = unit(1, 'in'),
        grid_width = unit(0.125, 'in')
    )
  )
}

# heatmap coloring from figure 1

ramp <- scales::colour_ramp(c("black", "red", 'white', "blue"))
colorramp <- ramp(seq(0, 1, length = 13))
# show_col(colorramp)
colors <- c(colorramp[c(1,3,7)], rep(colorramp[c(9)],3), colorramp[10])
col_fn <- colorRamp2(c(-6.8, -2.15, -1.2, -0.42, 0, 0.42, 1.2), colors)


# make heatmaps
hm1 <- draw_fitness_heatmap(
            mat=mat, name='Fitness', 
            col=col_fn, label_WT=T, sequence=gsp1_seq[toxics],
            gray_out=low_read_mat
        )

## make binary matrix for GTPase regions

gtpase_regions <- list(
    'P-loop'=seq(19,26),                       # Vetter chapter in Wittinghoffer book 2014
    'Extended Switch I (Vetter)'=seq(34,47),   # Vetter et al. Cell 1999
    'Extended Switch II (Vetter)'=seq(67,81),  # Vetter et al. Cell 1999 + Mg contact 67
    'N/TKxD motif'=seq(124,127),               # Vetter chapter in Wittinghoffer book 2014
    'SAK motif'=seq(152,154)                   # Vetter chapter in Wittinghoffer book 2014
)
gtpase_regions <- list('GTPase regions'=unlist(gtpase_regions, use.names=F))

# nucleotide_recognition <- list(
#     'NK.D to SAK interactions'=c(92,100,127,131,138,141,148,150,152)
# )

cterm <- list('C-terminal switch'=seq(181,208))

distal_sites <- list(
    'Distal sites in H32 mechanism'=c(28,30,31,32,33,35,48,50,54,156,159,163,183,186),
    'Allosteric positions (31P NMR)'=c(34,141,147,157))
distal_sites <- list('Distal sites affecting switching'=unlist(distal_sites, use.names=F))

contacts_nucleotide <-   # within 4 angstroms of nucleotide from PDB 3m1iA 
    list('Contacts Nucleotide'=unlist(scan('data/residues_contacting_nucleotide.txt', comment.char='#'), use.names=F))

# from BIOGRID: https://thebiogrid.org/31559/protein
PTMs_detailed <- list(
    'Phosphorylation site'=c(2,155,181),  # https://www.yeastgenome.org/locus/S000004284/protein
    'Succinylation site'=c(25),
    'Ubiquitinylation site'=c(101,125))

PTMs <- list('PTMs'=unlist(PTMs_detailed, use.names=F))

binary_df_gtpase <-
    enframe(c(
        gtpase_regions,
        contacts_nucleotide, 
        cterm,
        distal_sites,
        PTMs)) %>% 
    unnest(cols=c(value)) %>% 
    rename('position'=value) %>% 
    mutate(bool=1) %>% 
    unique() %>% 
    pivot_wider(names_from=name, values_from = bool) %>% 
    replace(is.na(.), 0)


## make binary matrix for regulator interfaces 

interfaces <-
    read_tsv('data/Gsp1_interfaces_SASA_and_conservation.txt', col_types=cols()) %>% 
    filter(protein == 'GSP1', file_path != 'PyMOL_figures/pdbs/clean/3m1i2.pdb') %>%  # see below
    select(partner, yeast_num, interface)    

#########
# Annotate as interface: T56 (Yrb1), D93 (GEF), T99 (GEF), V133 (GEF)
# Though these are not computed as in the interface, because 
# the side chains face inward, the residues around them form 
# the interface, and mutations at these positions would likely
# disrupt the interface geometry
# 
# Also for D93, Rojas et al 2012 (from Alfonso Valencia's group)
# identifies Ran D91 as a "specificity determining residue"
# because it is in the GEF interface and they also label it as
# "coordinating communication between different lobes" along with
# Ran M89, both of which are in the LV.D region (M insted of V)
# defined in Neuwald et al 2003.
# 
# there is no experimental data suggesting they are PTMs

df_interface <- 
    data.frame('yeast_num'=seq(1,220)) %>% 
    left_join(interfaces, by='yeast_num') %>% 
    complete(yeast_num, partner) %>% 
    filter(!is.na(partner)) %>% 
    mutate(interface = ifelse(is.na(interface), 'not_interface', interface))


binary_df_interface <-
    df_interface %>% 
    filter(interface=='core', partner %in% c('SRM1','RNA1','YRB1')) %>% 
    rename('position'='yeast_num') %>% 
    mutate('Regulator Interface'=1) %>% 
    select(position, 'Regulator Interface') %>% 
    unique() %>% 
    add_row('position' = 56, 'Regulator Interface' = 1) %>% 
    add_row('position' = 93, 'Regulator Interface' = 1) %>% 
    add_row('position' = 99, 'Regulator Interface' = 1) %>% 
    add_row('position' = 133, 'Regulator Interface' = 1)


binary_df <-
    data.frame('position'=seq(1,220,1)) %>% 
    left_join(binary_df_gtpase, by = 'position') %>% 
    left_join(binary_df_interface, by = 'position') %>% 
    replace(is.na(.), 0)
    
# write out for venn barplots
write_csv(binary_df, 'scripts/Fig4/toxic_annotations.txt')

# sets for PyMOL script for 4B
x <- sort(toxics)
y <- sort(c(unlist(gtpase_regions, use.names=F), unlist(contacts_nucleotide, use.names=F)) )
paste(setdiff(x,y), collapse='+')

binary_mat <- 
    binary_df %>% 
    filter(position %in% toxics) %>% 
    column_to_rownames('position') %>% 
    as.matrix() %>% 
    t()


# order by severity
binary_mat <- binary_mat[,as.character(toxics)]

# save binary mat for venn bar graphs



novel_positions <- which(! as.integer(colnames(binary_mat)) %in% c(unlist(gtpase_regions), unlist(cterm)))

hm2 <- Heatmap(
    mat=binary_mat, name='Categories',
    col=c('white','white'),
    show_heatmap_legend = F,
    na_col='white',

    width = ncol(binary_mat)*unit(1.8, "mm"), 
    height = nrow(binary_mat)*unit(1.8, "mm"),

    row_title = 'Category',
    cluster_rows = F,
    row_names_side = 'right',
    row_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
    row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),

    column_title = 'Gsp1 sequence position',
    cluster_columns = F,
    column_names_side = 'top',
    column_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
    column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),

    # borders around all cells, green if synonymous mutation
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(j %in% novel_positions) {
            grid.rect(x, y, width, height, gp = gpar(col = 'black', lwd = 0.5, fill = "red3")) 
        } else {
            grid.rect(x, y, width, height, gp = gpar(col = 'black', lwd = 0.5, fill = "transparent")) 
        }
        if(binary_mat[i,j]==1) {
            grid.circle(x, y, r = unit(0.25, "mm"), gp = gpar(fill = 'black'))
        }

      }
)

pdf('figures/Fig4/Fig4_toxics_annotation.pdf', width=6, height=2.5)
hm1 %v% hm2
dev.off()
