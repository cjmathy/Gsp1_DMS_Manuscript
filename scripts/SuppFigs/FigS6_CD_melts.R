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
         Tolerant = c('F28Y','H50R','F54W','N156W','F159W','F163Y'),
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


all_melt_filepaths <- dir(path = 'data/CD', pattern = '*Melt.txt', full.names = T, recursive = T)
index <- read_delim(file.path('data/CD', 'index.txt'), delim='\t', col_types=cols()) %>% 
    mutate(filepath = map(filename, grep, all_melt_filepaths, value=TRUE)) %>% 
    unnest(filepath)
    

melts_df <- 
    index %>% 
    pull(filepath) %>% 
    lapply(function (x) parse_CD_file(x)) %>% 
    bind_rows() %>%
    rename('temperature' = `Temperature [C]`) %>% 
    left_join(index, by='filename') %>% 
    mutate('molar_ellipticity' = pmap_dbl(
        ., function(CD_mdeg, path_length_mm, conc_uM, ...) mdeg_to_mol_ellip(
            vec = CD_mdeg, p = path_length_mm, c = conc_uM)
        )) %>% 
    left_join(df_tox, by = 'mutant') %>% 
    mutate(sample = factor(sample, levels=sample_order)) %>% 
    mutate(seqpos = ifelse(position == 'WT', 0, strtoi(substring(position, 2)))) %>% 
    mutate(position = fct_reorder(position, seqpos))  %>%                
    group_by(filename) %>% 
    mutate(starting_mol_ell = if_else(temperature == 25, molar_ellipticity, NA_real_)) %>% 
    mutate(starting_mol_ell = min(starting_mol_ell, na.rm = TRUE)) %>% 
    mutate(max_mol_ell = max(molar_ellipticity)) %>% 
    ungroup() %>% 
    mutate(mol_ell_clean = (molar_ellipticity - starting_mol_ell)/(max_mol_ell - starting_mol_ell)) %>% 
    filter(sample %in% sample_order)

fits_df <- 
    melts_df %>%
    select(sample, mutant, date, position, filename, mutant_type, temperature, mol_ell_clean) %>% 
    group_by(filename) %>% 
    nest(data=c(temperature, mol_ell_clean)) %>%
    mutate(fit = map2(.x = data, .y = filename, ~ fit_sigmoid_melt(.x, x ='temperature', y='mol_ell_clean', id = .y))) %>% 
    mutate(tidied = map(fit, tidy))

df_to_plot <-
    fits_df %>% 
    unnest(tidied) %>% 
    mutate(augmented = map(fit, augment)) %>% 
    unnest(augmented) %>% 
    select(sample, mutant, date, position, filename, mutant_type, term, estimate, x, y, .fitted) %>% 
    pivot_wider(names_from = 'term', values_from = 'estimate')


ggplot(df_to_plot, aes(x)) +
    geom_point(aes(y=y, color=mutant_type), size=0.5, stroke=0.5) +
    scale_color_manual(values=c('red','blue','gray'), name='Mutant Type') +
    geom_line(aes(y=.fitted)) +
    geom_text(data=filter(df_to_plot, x==25),
              aes(x=45, y = 0.9, label = paste0('Tm = ', round(Tm, 1))),
              size = 2, family='Helvetica'
              ) +
    xlab('Temperature (C)') + ylab('Fraction unbound') +
    facet_wrap(vars(sample), ncol=4, scales='free') +
    theme_custom

ggsave('figures/FigS6/FigS6_melts.pdf', width = 7, height = 4)
 

df_fit_params <- fits_df %>% 
    unnest(tidied) %>% 
    mutate(augmented = map(fit, augment)) %>% 
    unnest(augmented) %>% 
    select(sample, mutant, date, position, filename, mutant_type, term, estimate) %>% 
    unique() %>% 
    pivot_wider(names_from = 'term', values_from = 'estimate') %>% 
    arrange(sample) %>% 
    ungroup(filename)

write_csv(df_fit_params, 'data/CD_melt_fits.csv')