# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# load additional package(s)
library(ComplexHeatmap)
library(scales)

df <- 
    read_csv('data/Gsp1_fitness_scores.csv', col_types = cols()) %>% 
    mutate('bin' = factor(bin, levels=BIN_ORDERING),
           'aa_to'=factor(aa_to, levels=AA_ORDERING),
           'aa_from'=factor(aa_from, levels=AA_ORDERING))

# prepare a sequence vector for WT Gsp1 so we can box cells in heatmaps that are WT in sequence (i.e. T34T)
gsp1_seq <- select(df, aa_from, position) %>% unique() %>% pull(aa_from) %>% as.character()

# convert dataframes to matrices for each dataset, for plotting
mat <- 
  df %>%
  select(position, aa_to, score) %>%
  pivot_wider(id_cols = aa_to, names_from = position, values_from = score) %>%
  arrange(desc(aa_to)) %>%
  column_to_rownames('aa_to') %>%
  as.matrix()

# prepare a matrix for the positions with low read counts, which should be gray
mat_low_read <-
  df %>% 
  select(position, aa_to, low_reads_flag) %>%
  pivot_wider(id_cols = aa_to, names_from = position, values_from = low_reads_flag) %>%
  arrange(desc(aa_to)) %>%
  column_to_rownames('aa_to') %>%
  as.matrix()

##### HEATMAPS #####

# function for plotting a heatmap of empiric data
# input dataset should be a matrix, with 21 aa labels as rows
# and each residue position as columns
draw_fitness_heatmap <- function(mat, name, col, label_WT=F, sequence=NULL, low_read=mat_low_read) {
  Heatmap(
    mat,
    name=name,
    col=col,
    show_heatmap_legend = T,
    na_col='black',
          
    row_title = 'Amino Acid Substituted',
    cluster_rows = F,
    row_names_side = 'left',
    row_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
    row_names_gp = gpar(fontsize = 6, fontfamily='Courier'),

    column_title = 'Gsp1 sequence position',
    cluster_columns = F,
    column_names_side = 'top',
    column_title_gp = gpar(fontsize = 7, fontfamily='Helvetica'),
    column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
    column_labels = replace(colnames(mat), seq(from = 2, to = length(colnames(mat)), by = 2), ''),

    cell_fun = function(j, i, x, y, width, height, fill) {

        # Cells with low_reads_flag = T: white, with x through them
        if (low_read[i,j]) {
            grid.rect(x, y, width, height, gp = gpar(fill = 'white', lwd=0.5))
            grid.lines(x=c(x-width/2, x+width/2), y=c(y+height/2, y-height/2),
                       default.units = 'npc', gp=gpar(lwd=0.2))
            grid.lines(x=c(x-width/2, x+width/2), y=c(y-height/2, y+height/2),
                       default.units = 'npc', gp=gpar(lwd=0.2))
            # grid.rect(x, y, width, height, gp = gpar(fill = 'gray', lwd=0.5))
        }
        
        # WT synonymous cells: black dot inside of them
        if ((rownames(mat)[i] == sequence[j]) & label_WT) {
            grid.circle(x, y, r = min(unit.c(width, height)/2), gp = gpar(fill = 'black'))
        }

        # borders around all cells
        grid.rect(x, y, width, height, gp = gpar(col = 'black', lwd = 0.5, fill = "transparent")) 
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

### Get statistics of dataset to set color scheme

# max of heatmap when ignoring the low reads 
df %>% 
  filter(! low_reads_flag) %>% 
  summarise(
    max(score, na.rm=T),  # 1.61
    min(score, na.rm=T)   # -12.0
  )

# 

# this includes WTsyn and STOP
# just used for setting color gradient
read_csv('data/Gsp1_fitness_scores.csv', col_types=cols()) %>%
    filter(!low_reads_flag) %>% 
    mutate(bin=factor(bin, levels=rev(BIN_ORDERING))) %>% 
    filter(! bin %in% c('WT','STOP')) %>% 
    group_by(bin) %>% 
    mutate(n=n(),
           mean_score=mean(score, na.rm=T),
           sd = sd(score, na.rm=T)) %>% 
    select(bin, n, mean_score, sd) %>% 
    unique


# # A tibble: 5 x 4                  
# # Groups:   bin [5]
#   bin              n mean_score    sd
#   <fct>        <int>      <dbl> <dbl>
# 1 STOP-like      677    -2.15   0.431
# 2 WT-like       2396    -0.0289 0.421
# 3 intermediate   203    -1.21   0.123
# 4 beneficial      19     1.20   0.160
# 5 toxic         1243    -6.78   2.43 

# based on this, -2 should already be dark red (because it is deleterious at the level of STOP mutations)

ramp <- colour_ramp(c("black", "red", 'white', "blue"))
colorramp <- ramp(seq(0, 1, length = 13))
# show_col(colorramp)
colors <- c(colorramp[c(1,3,7)], rep(colorramp[c(9)],3), colorramp[10])
col_fn <- colorRamp2(c(-6.8, -2.15, -1.2, -0.42, 0, 0.42, 1.2), colors)

# make heatmaps 
hm1 <- draw_fitness_heatmap(mat=mat[1:21,1:110], name='Fitness', col=col_fn, label_WT=T,
                            sequence=gsp1_seq[1:110], low_read=mat_low_read[1:21,1:110])
hm2 <- draw_fitness_heatmap(mat=mat[1:21,111:220], name='Fitness', col=col_fn, label_WT=T,
                            sequence=gsp1_seq[111:220], low_read=mat_low_read[1:21,111:220])

pdf('figures/Fig1/Fig1C_top.pdf', height = 1.85, width = 7.3)
draw(hm1)
dev.off()

pdf('figures/Fig1/Fig1C_bot.pdf', height = 1.85, width = 7.3)
draw(hm2)
dev.off()

# annotation of secondary structure done in Illustrator, using the annotations from
# crystal structures at:
#   Scheffzek, K., Klebe, C., Fritz-Wolf, K., Kabsch, W., & Wittinghofer, A. 
#     (1995). Crystal structure of the nuclear Ras-related protein Ran in its 
#     GDP-bound form. Nature, 374(6520), 378-381.
#   Vetter, I. R., Nowak, C., Nishimoto, T., Kuhlmann, J., & Wittinghofer, A. 
#     (1999). Structure of a Ran-binding domain complexed with Ran bound to a
#     GTP analogue: implications for nuclear transport. Nature, 398(6722), 39-46.