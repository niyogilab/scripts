#!/usr/bin/env Rscript

# TODO: put antibiotic + concentration in one col from now on,
#       mainly because that's how you apply them to plates.

suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(dplyr))

### helpers for reading/writing files ###

read_text <- function(fileName)
  # from https://stackoverflow.com/questions/9068397
  readChar(fileName, file.info(fileName)$size)

read_plates <- function(filename)
  read.csv(filename, colClasses='character') %>% select(-notes)

read_map <- function(filename, colname) {
  # TODO: check that wells match the generated list
  tmp <- tempfile(fileext='.csv')
  system(paste('tidymap', colname, filename, tmp))
  read.csv(tmp, colClasses='character')
}

write_wells <- function(df, filename)
  write.csv(df, file=filename, quote=FALSE, row.names=FALSE)

### helpers for matching up data ###

plate_dimensions <- function(ptype) {
  # Look up dimensions of known plates
  # TODO: is there any better way?
  dimensions = list(
    '06well' = c(2,3),
    '24well' = c(4,6),
    '48well' = c(6,8),
    '96well' = c(8,12)
  )
  stopifnot(ptype %in% names(dimensions))
  dimensions[[ptype]]
}

list_wells <- function(ptype) {
  # Lists all the wells that should be present in a plate of this type
  dims <- plate_dimensions(ptype)
  rows <- LETTERS[1:dims[1]]
  c(sapply(rows, function(r) paste0(r, 1:dims[2])))
}

matches_keywords <- function(filename, keywords)
  all(sapply(keywords, function(w) grepl(w, filename)))

matching_files <- function(paths, keywords)
  paths[sapply(paths, function(p) matches_keywords(p, keywords))]

### merge algorithm ###

merge_reference <- function(maps, colname, value) {
  # Finds and loads a referenced file by colname + value.
  colname <- gsub("s$", "", colname)
  value   <- gsub("^@", "", value)
  matches <- matching_files(maps, c(colname, value))
  stopifnot(length(matches) == 1)
  read_map(matches, colname) # TODO: error here if wrong number of wells
}

merge_literal <- function(wells, colname, cell) {
  # Contorts a single value into a data frame for easier merging
  df <- data.frame(well=wells)
  df[[colname]] <- cell
  df
}

merge_cell <- function(maps, wells, colname, cell) {
  # Expands one cell of the plate frame into a column indexed by well.
  # References (start with '@') are loaded, and literals recycled to fit.
  if (length(strsplit(cell, '@')[[1]]) > 1)
    merge_reference(maps, colname, cell)
  else
    merge_literal(wells, colname, cell)
}

merge_row <- function(maps, row) {
  # Expands one row of the plate frame into its own frame with rows for wells
  # stopifnot(length(row) == length(colnames))
  each_n <- function(n)
    merge_cell(maps, list_wells(row[['type']]), colnames(row)[n], row[,n])
  seq(row) %>%
    lapply(each_n) %>%
    merge_recurse(by='well')
}

merge_plate <- function(maps, plates)
  # Expands the plate frame into a larger one with a row for every well
  seq(nrow(plates)) %>%
    lapply(function(n) merge_row(maps, plates[n,])) %>%
    bind_rows()

### command line interface ###

main <- function(args) {
  maps   <- args[['map']]
  plates <- read_plates(args[['plate']])
  wells  <- merge_plate(maps, plates)
  if (is.null(args[['out']]))
    out <- stdout()
  else
    out <- args[['out']]
  write_wells(wells, out)
}

# the @var@ gets replaced by Nix at build time
main(docopt(read_text('@usageTxt@')))
