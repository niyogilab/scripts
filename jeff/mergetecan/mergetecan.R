#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(dplyr))

read_text <- function(filename)
  # Read a file as plain text, which is awkward in R.
  # See https://stackoverflow.com/questions/9068397
  readChar(filename, file.info(filename)$size)

read_df <- function(filename)
  # Read a data frame from a CSV file.
  read.csv(filename, colClasses='character') %>% tbl_df()

write_df <- function(df, filename)
  # Write the merged CSV file.
  write.csv(df, file=filename, quote=FALSE, row.names=FALSE)

read_many <- function(filenames)
  # Read a list of csv files and bind them into one big data frame.
  # Note that they should all have the same columns!
  filenames %>%
    lapply(read_df) %>%
    bind_rows()

merge_tecan <- function(meta, tecan)
  # Keep all the tecan rows, but add metadata
  left_join(tecan, meta, by=c('plate', 'well'))

main <- function(args) {
  meta   <- read_df(args[['meta']])
  tecan  <- read_many(args[['tecan']])
  merged <- merge_tecan(meta, tecan)
  if (is.null(args[['--out']]))
    out <- stdout()
  else
    out <- args[['--out']]
  write_df(merged, out)
}

# to debug, comment out the last line and run main interactively:
# main(list(meta='test/meta.csv',
#           tecan=c('test/tecan1.csv', 'test/tecan2.csv')))

# the @var@ gets replaced by Nix at build time
main(docopt(read_text('@usageTxt@')))
