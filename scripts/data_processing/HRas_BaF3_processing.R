source('scripts/config_workspace.R')

# load additional package(s)
library(ComplexHeatmap)


# On 01/31/2022 downloaded Ras data from https://github.com/fhidalgor/ras_cancer_hidalgoetal/, commit 0dcb01b, using:
# wget https://github.com/fhidalgor/ras_cancer_hidalgoetal/blob/main/Enrichments_DMS/HRas188_BaF3_r0.csv
# wget https://github.com/fhidalgor/ras_cancer_hidalgoetal/blob/main/Enrichments_DMS/HRas188_BaF3_r1.csv

aa_order_hidalgo <- unlist(strsplit('DEKHRGNQASTPCVYMILFW*', split = ""))

df_seq <- read_delim('data/Ras_data/HRas_Gsp1_index.txt', delim='\t', show_col_types=F)

HRas_BaF3_seq <-
    df_seq %>% 
    filter(HRas_res_num > 1, HRas_res_num < 161) %>% 
    pull(HRas)

read_Ras_data <- function(file_in) {
    read_csv(file_in, col_names=as.character(seq(2,160)), show_col_types=F) %>% 
    mutate(aa=factor(
        unlist(strsplit('ACDEFGHIKLMNPQRSTVWY*', split = "")),
        levels=aa_order_hidalgo)) %>% 
    select(aa, everything()) %>% 
    arrange(aa) %>% 
    column_to_rownames('aa') %>% 
    as.matrix()
}


data_r0 <- read_Ras_data('data/Ras_data/HRas188_BaF3_r0.csv')
data_r1 <- read_Ras_data('data/Ras_data/HRas188_BaF3_r1.csv')

# one NA in r0, at 582. Position V29Y
which(is.na(data_r0))
which(is.na(data_r0), arr.ind=TRUE)

# 59 NAs in r1
which(is.na(data_r1))
length(which(is.na(data_r1)))
which(is.na(data_r1), arr.ind=TRUE)

# I can see values for many of these points in the Figure 2A heatmap, so I believe
# the authors took the value of one replicate if there was an NA
data_r0[is.na(data_r0)] = data_r1[is.na(data_r0)]
data_r1[is.na(data_r1)] = data_r0[is.na(data_r1)]


data <- (data_r0 + data_r1)/2

hm <- Heatmap(
    data,
    name='deltaE_i',
    col=colorRamp2(
        seq(-1, 1, by=0.5),
        rev(colorRampPalette(colors = c("red", "white", "blue"))(5))
    ),
    # col=colorRamp2(seq(-0.6, 0.6, by=0.2), RColorBrewer::brewer.pal(9, 'RdBu')[8:2]),
    show_heatmap_legend = T,
    na_col='black',
        
    
    cluster_rows = F,
    row_names_side = 'left',
    row_title_gp = gpar(fontsize = 10, fontfamily='Helvetica'),
    row_names_gp = gpar(fontsize = 10, fontfamily='Courier'),

    column_title = 'HRas sequence position',
    cluster_columns = F,
    column_names_side = 'top',
    column_title_gp = gpar(fontsize = 10, fontfamily='Helvetica'),
    column_names_gp = gpar(fontsize = 10, fontfamily='Helvetica'),
    column_labels = replace(colnames(data), seq(from = 2, to = length(colnames(data)), by = 2), ''),

    bottom_annotation = HeatmapAnnotation(
    which = 'column',
    'sequence' = anno_text(
        HRas_BaF3_seq, gp=gpar(fontsize = 10, fontfamily='Courier'),
        rot=0, just = "center", location = 0.5
        )
    ),

    cell_fun = function(j, i, x, y, width, height, fill) {
        # WT synonymous cells: black dot inside of them
        if ((rownames(data)[i] == HRas_BaF3_seq[j])) {
            grid.rect(x, y, width, height, gp = gpar(col = 'black', lwd = 0.5, fill = "transparent")) 
        }            
    },

)


# pdf('figures/Fig4/Ras_BaF3_hm.pdf', height = 3, width = 18)
draw(hm)
# dev.off()


data %>% 
    as.data.frame() %>% 
    rownames_to_column('aa') %>% 
    tibble() %>% 
    pivot_longer(cols=c(-aa), names_to='position', values_to='score') %>% 
    mutate(position=as.integer(position)) %>% 
    rename('aa_to'=aa) %>% 
    arrange(position) %>% 
    left_join(select(df_seq, 'position'=HRas_res_num, 'aa_from'=HRas), by='position') %>% 
    select(aa_from, position, aa_to, score) %>% 
    write_csv('data/Ras_data/HRas188_BaF3_processed.csv')


