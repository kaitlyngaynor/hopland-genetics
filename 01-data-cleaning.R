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
  select(lab_id, cond, stor, n, extracted_conc) %>% 
  rename(alleles_working = n,
         condition = cond,
         storage = stor)

# remove qpcr outliers
# change outliers to NA
genotypes$extracted_conc_nooutlier <- NA
for(i in 1:nrow(genotypes_publication)) {
  if(is.na(genotypes_publication$extracted_conc[i])) {
    genotypes_publication$extracted_conc_nooutlier[i] <-  NA
  }
  else if(genotypes_publication$extracted_conc[i] >= min(boxplot.stats(genotypes_publication$extracted_conc)$out)) {
    genotypes_publication$extracted_conc_nooutlier[i] <- NA
  } else {
    genotypes_publication$extracted_conc_nooutlier[i] <- genotypes_publication$extracted_conc[i]
  }
}

genotypes_publication <- genotypes_publication %>% 
  rename(qpcr_conc = extracted_conc_nooutlier) %>% 
  select(-extracted_conc)

write_csv(genotypes_publication, here::here("for-dryad", "odocoileus-fecal-genotype-data.csv"))
