# DATA CLEANING
# Kaitlyn Gaynor
# updated 2019-11-20

library(tidyverse)
library(here)


# Import data -------------------------------------------------------------

# Import genotypes and subset to triplicate records
genotypes <- read_csv(here::here("data", "20191112-triplicates-genotypes.csv")) %>%
  subset(trip == "yes")

# Import qpcr results
qpcr <- read_csv(here::here("data", "20191112-triplicates-qpcr-results.csv"))

# Merge genotypes and qpcr
genotypes_qpcr <- left_join(genotypes, qpcr)


# Add additional columns --------------------------------------------------

# Columns that indicates whether or not the sample worked for 8 or more alleles, and whether or not it worked for all 10

genotypes_qpcr$working <- "NA"
genotypes_qpcr$perfect <- "NA"

for (i in 1:nrow(genotypes_qpcr)) {
  if (genotypes_qpcr$n[[i]] > 7) {
    genotypes_qpcr$working[[i]] <- "yes"
  } else {
    genotypes_qpcr$working[[i]] <- "no"
  }
  if (genotypes_qpcr$n[[i]] == 10) {
    genotypes_qpcr$perfect[[i]] <- "yes"
  } else {
    genotypes_qpcr$perfect[[i]] <- "no"
  }
}


# Export cleaned data -----------------------------------------------------

write_csv(genotypes_qpcr, here::here("data", "cleaned-triplicates-genotypes-qpcr.csv"))
