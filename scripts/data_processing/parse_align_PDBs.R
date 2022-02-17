#!/usr/bin/env RScript

# run from command line using 'Rscript parse_align_PDBs.R' from
# the scripts directory of the project

## Author: Christopher Mathy
## Date: 2020-02-04
## Email: cjmathy@gmail.com
## Email: chris.mathy@ucsf.edu
## Description:
##  This script preprocesses structures of Ran GTPase downloaded
##  using the script 'download_data.sh'. It parses and aligns
##  PDBs using the packages bio3d

# load modules
library(tidyverse)
library(bio3d)

# read structural info file
df <- read_delim('data/pdb_structures_info.txt', delim = '\t', col_types = cols())

### ---------------------------------------------------------------------------
### Clean raw PDBs to have complexes with Ran as chain A and partner as chain B
### Also write out PDBs of monomeric Ran

# list of raw PDBs downloaded from the web (using download_data.sh)
files <- list.files(path = 'data/pdbs_raw', full.names = T)

# set output directories
cmplxdir <- 'data/pdbs_complexes'
randir <- 'data/pdbs_ran'
dir.create(cmplxdir, showWarnings = FALSE)
dir.create(randir, showWarnings = FALSE)

# iterate through each PDB using the mapping function purrr::pwalk
# pmap and pwalk iterate through the rows of a dataframe and perform a function
# pwalk just avoids returning NULL, since we aren't keeping the returned values
# of the anonymous function
df %>%
  dplyr::select(pdb_id, ran_chain, partner_chain, partner_name) %>%
  purrr::pwalk(.f = function (pdb_id, ran_chain, partner_chain, partner_name) {

    # read in file
    file <- grep(pdb_id, files, value=T)
    print(paste0('Processing ', file))
    raw_pdb <- read.pdb(file)

    # split pdb (using bio3d functions)
    ran <- trim.pdb(raw_pdb, chain=ran_chain)
    partner <- trim.pdb(raw_pdb, chain=partner_chain)
    complex <- cat.pdb(ran, partner, rechain=T)  # warnings about chain breaks OK to ignore

    # write ran pdb
    write.pdb(pdb=ran, file = paste0(randir, '/', pdb_id, '.pdb'))

    # write complex pdb
    outfile <- ifelse(!is.na(partner_name),
                      paste0(cmplxdir, '/', pdb_id, '_', partner_name, '.pdb'),
                      paste0(cmplxdir, '/', pdb_id, '.pdb'))
    write.pdb(pdb=complex, file = outfile)
    }
  )

### ---------------------------------------------------------------------------
### Align the structures of the complexes using bio3d and MUSCLE
### Then do the same for just the Ran monomers

structure_align <- function(fdir) {
  files <- list.files(fdir, full.names = T)           # read all the pdb files in one list object
  pdbs <- pdbaln(files, outfile = 'data/Ran_aln.fa')  # multiple sequuence alignment
  core <- core.find(pdbs)                             # find the conserved core residues (which don't move much)
  core.inds <- print(core, vol = 0.5)                 # indices of residues to be aligned
  pdbfit(pdbs, core.inds, outpath=fdir)               # structural superposition and write to folder

  # delete unaligned files and rename aligned files
  unaligned_files <- grep(list.files(path=fdir, full.names=T), pattern='pdb_flsq', inv=T, value=T)
  unlink(unaligned_files)
  aligned_files <- list.files(path=fdir, full.names=T)
  new_filenames <- sapply(aligned_files, gsub, pattern = '.pdb_flsq', replacement = '')
  file.rename(from=aligned_files, to=new_filenames)
}

# align and write both the complexes and the monomeric Ran
structure_align(cmplxdir)
structure_align(randir)

