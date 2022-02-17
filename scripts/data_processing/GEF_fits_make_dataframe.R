library(tidyverse)
library(broom)

source('scripts/config_workspace.R')

sample_order <- c('WT','F28V','F28Y','H50N','H50R','F54A','F54W',
                  'N156A','N156W','F159L','F159W','F163L','F163Y')

read_assay_file <- function(file) {
  name <- str_split(file, '/')[[1]][3]
  date <- str_split(name, '_')[[1]][1]
  plate <- str_split(name, '_')[[1]][4] %>% substr(1,nchar(.)-4)

  stream <- read_lines(file)
  first_data_row <- grep(stream, pattern = "Time\tT")
  
  read_tsv(
      file, col_types=cols(),
      col_names = T,
      skip = (first_data_row - 1),
      locale = locale(encoding = 'windows-1250')) %>% 
    mutate(time=as.numeric(hms(Time))) %>%
    select(-Time) %>% 
    rename(temp=`T° 295,335`) %>% 
    pivot_longer(
      cols = -c(time, temp),
      names_to = 'well',
      values_to = 'fluorescence') %>% 
    mutate(date = factor(date), plate=as.numeric(plate)) %>% 
    select(date, plate, time, temp, well, fluorescence)
}

read_setup_file <- function(file) {
  sample_order <- c('WT','F28V','F28Y','G35D','G35A',
                    'H50N','H50R','F54A','F54W','N156A','N156W',
                    'F159L','F159W','F163L','F163Y')  
  read_tsv(file, col_types=cols()) %>%
    mutate(
      date=factor(date),
      sample=factor(sample, levels=sample_order),
      nucleotide=factor(nucleotide, levels=c('mant-dGTP','mant-dGDP')),
      replicate=factor(replicate),
      GEF_conc=ifelse(replicate==0, 0, 2.5) %>% as.factor()
    )
}

setup_paths <- dir(path = 'data/enzymology', pattern = '20211213_GEF_setup*|20220118_GEF_setup*', full.names = T, recursive = T)
setup_df <- lapply(setup_paths, function (x) read_setup_file(x)) %>% bind_rows

data_paths <- dir(path = 'data/enzymology', pattern = '*20211213_GEF_assay*|20220118_GEF_assay*', full.names = T, recursive = T)
data_df <- lapply(data_paths, function (x) read_assay_file(x)) %>% 
  bind_rows %>% 
  mutate(plate = ifelse(is.na(plate),1,plate))

# read in parameters from data fitting using GEF_assay_shiny_app
df_fits <- 
  read_delim('data/enzymology/GEF_assay_data_fits.txt', delim='\t', col_types = cols()) %>%   
  separate(id, into=c('date','plate','well','sample'), remove=FALSE) %>% 
  mutate(
      date=factor(date),
      plate=as.numeric(plate),
      sample=factor(sample, levels=sample_order),
      nucleotide=factor(nucleotide, levels=c('mant-dGTP','mant-dGDP')),
      replicate=factor(replicate)) %>% 
  group_by(date, plate, sample, nucleotide) %>% 
  mutate(group_id=cur_group_id()) %>% 
  arrange(group_id) %>% 
  ungroup()

df_fits_noenzyme <-
  df_fits %>% 
  filter(replicate==0) %>% 
  select(date, plate, sample, nucleotide, group_id, 'slope_no_enzyme'=slope)
  
df_fits_to_merge <- 
  df_fits %>% 
  left_join(df_fits_noenzyme, by=c('date','plate','sample','nucleotide','group_id')) %>% 
  select(id, intercept, slope, slope_no_enzyme, t1, t2, fittype)

# merge parameters from fits into main dataframe
df <-
  left_join(setup_df, data_df, by=c('date','plate','well')) %>% 
  mutate(id = paste(date, plate, well, sample, sep='_')) %>% 
  left_join(df_fits_to_merge, by='id')

write_csv(df, 'data/GEF_assay_data_and_fits.csv')

