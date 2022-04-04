# load tidyverse and circlize, and set theme for plotting
# also loads global variables for plotting
source('scripts/config_workspace.R')

# load CD functions
source('scripts/data_processing/CD_functions.R')

# load additional packages
library(scales)
library(ggrepel)
library(broom)

df_tox <- 
    list(`Toxic/GOF` = c('F28V','H50N','F54A','N156A','F159L','F163L'), 
         `WT-like` = c('F28Y','H50R','F54W','N156W','F159W','F163Y'),
         WT = c('WT')
         ) %>% 
    stack() %>% 
    select('mutant' = values, 'mutant_type' = ind)

sample_order <- c('F28V GDP','F28Y GDP',
                  'H50N GDP','H50R GDP', 
                  'F54A GDP','F54W GDP',
                  'N156A GDP','N156W GDP',
                  'F159L GDP','F159W GDP',
                  'F163L GDP','F163Y GDP',
                  'WT GDP')

scan_filepaths1 <- dir(path = 'data/CD', pattern = '*Scan-BufferSubtracted.txt', full.names = T, recursive = T)
scan_filepaths2 <- dir(path = 'data/CD', pattern = '*Scan_BufferSubtracted.txt', full.names = T, recursive = T)
all_scan_filepaths <- c(scan_filepaths1, scan_filepaths2)
index <- read_delim(file.path('data/CD', 'index.txt'), delim='\t', col_types=cols()) %>% 
    mutate(filepath = map(filename, grep, all_scan_filepaths, value=TRUE)) %>% 
    unnest(filepath)

scans_df <- 
    index %>% 
    pull(filepath) %>% 
    lapply(function (x) parse_CD_file(x)) %>% 
    bind_rows() %>%
    rename('wavelength' = NANOMETERS) %>% 
    left_join(index, by='filename') %>% 
    mutate('molar_ellipticity' = pmap_dbl(
        ., function(CD_mdeg, path_length_mm, conc_uM, ...) mdeg_to_mol_ellip(
            vec = CD_mdeg, p = path_length_mm, c = conc_uM)
        )) %>% 
    left_join(df_tox, by = 'mutant') %>% 
    mutate(sample = factor(sample, levels=sample_order)) %>% 
    mutate(seqpos = ifelse(position == 'WT', 0, strtoi(substring(position, 2)))) %>% 
    mutate(position = fct_reorder(position, seqpos)) %>% 
    as_tibble() %>% 
    filter(sample %in% sample_order)


ggplot(scans_df) +
    geom_point(aes(x=wavelength, y= molar_ellipticity, color=mutant_type), size=0.5, stroke=0.5) +
    scale_color_manual(values=c('red','blue','gray'), name='Mutant Type') +
    xlab('Wavelength') + ylab('Molar ellipticity [deg * cm * dmol-1]') +
    facet_wrap(vars(sample), ncol=4, scales='free') +
    ylim(c(-3e5, 3e5)) + 
    theme_custom

ggsave('figures/FigS5/FigS5_scans.pdf', width = 7, height = 4)
