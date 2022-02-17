# SET WORKING DIRECTORY
setwd('/Users/cjmathy/gdrive/gsp1_dms/')

library(tidyverse)
suppressPackageStartupMessages(library(circlize))

# set theme for plotting
# font sizes
SMALL_SIZE = 5
MEDIUM_SIZE = 6
BIG_SIZE = 7

theme_custom <- theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = BIG_SIZE),
    axis.title = element_text(size = BIG_SIZE),
    axis.text = element_text(size = BIG_SIZE),
    axis.ticks = element_line(size = 0.2),
    axis.ticks.length = unit(0.05, 'cm'),
    axis.line = element_line(size = 0.1),
    strip.text.x = element_text(size = SMALL_SIZE)
  )

# set theme for drafting plots
theme_custom_drafting <- theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = BIG_SIZE+5),
    axis.title = element_text(size = BIG_SIZE+5),
    axis.text = element_text(size = BIG_SIZE+5),
    axis.ticks = element_line(size = 0.2),
    axis.ticks.length = unit(0.05, 'cm'),
    axis.line = element_line(size = 0.1),
    strip.text.x = element_text(size = SMALL_SIZE+5)
  )

# GLOBAL VARIABLES
AA_ORDERING <- c('W','Y','F','M','L','I','V','A','C','G',
                'P','T','S','Q','N','E','D','H','R','K','*')
BIN_ORDERING <- c('STOP','WT','toxic','STOP-like',
                  'intermediate','WT-like','beneficial')
BIN_COLORS <- setNames(as.list(
    c('black','forestgreen',RColorBrewer::brewer.pal(7, 'RdBu')[1:5])),
    BIN_ORDERING)
