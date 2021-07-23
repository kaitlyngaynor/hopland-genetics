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



# Export cleaner version for publication ----------------------------------

genotypes_publication <- genotypes_qpcr %>% 
  select(lab_id, cond, stor, n, working, extracted_conc) %>% 
  rename(qpcr_conc = extracted_conc,
         alleles_amplified = n,
         condition = cond,
         storage = stor) %>% 
  # create coarser condition column
  mutate(condition_coarse = fct_collapse(as.factor(condition),
         Slimy = c("1", "1.5"),
         Wet = c("2", "2.5"),
         Shiny = c("3", "3.5"),
         Dull = c("4")))

write_csv(genotypes_publication, here::here("for-dryad", "odocoileus-fecal-genotype-data.csv"))
