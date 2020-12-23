library(here)
library(tidyverse)
library(allelematch)

genotypes <- read_csv(here::here("data", "cleaned-triplicates-genotypes-qpcr.csv"))

# determine which genotype has most alleles, use this
genotypes_keep <- genotypes %>% 
  select(lab_id, stor, n) %>% 
  pivot_wider(names_from = stor, values_from = n) %>% 
  gather(stor, n, swab:EtOH) %>% 
  group_by(lab_id) %>% 
  slice(which.max(n)) %>% 
  mutate(keep = "TRUE")

genotypes_fewer <- left_join(genotypes, genotypes_keep) %>% 
  filter(keep == "TRUE") %>% 
  filter(n > 7) # take only the ones with at least 8 working

nrow(genotypes_fewer) # 99 samples with usable (across storage types)

# run allelematch
genotypes_allelematch <- genotypes_fewer %>% 
  select(lab_id, `1_TGLA94`:`2_CELB9`, -`SRY`) %>% # select all disomic loci and lab ID
  amDataset(indexColumn = "lab_id") %>%  # get into alleleMatch format
  amUnique(alleleMismatch = 4) # assign a unique individual ID to each of the samples based on genotype overlap

summary(genotypes_allelematch, csv="data/allelematch_output.csv")

# read back in (can't figure out better way?)
allelematch_output <- read.csv("data/allelematch_output.csv")

# how many unique individuals?
length(unique(allelematch_output$uniqueIndex)) # 68


# determine distribution of individuals
id_count <- allelematch_output %>% 
  count(uniqueIndex) %>% 
  arrange(-n) 

(id_count_count <- id_count %>% 
  rename(samples = n) %>% 
  count(samples))


# only dry ----------------------------------------------------------------

genotypes_fewer_dry <- genotypes %>% 
  filter(stor == "dry") %>% 
  filter(n > 7) # take only the ones with at least 8 working

nrow(genotypes_fewer_dry) # 67 samples with usable genotypes

# run allelematch
genotypes_allelematch_dry <- genotypes_fewer_dry %>% 
  select(lab_id, `1_TGLA94`:`2_CELB9`, -`SRY`) %>% # select all disomic loci and lab ID
  amDataset(indexColumn = "lab_id") %>%  # get into alleleMatch format
  amUnique(alleleMismatch = 4) # assign a unique individual ID to each of the samples based on genotype overlap

summary(genotypes_allelematch_dry, csv="data/allelematch_dry_output.csv")

# read back in (can't figure out better way?)
allelematch_output_dry <- read.csv("data/allelematch_dry_output.csv")

# how many unique individuals?
length(unique(allelematch_output_dry$uniqueIndex)) # 47

# only etoh ----------------------------------------------------------------

genotypes_fewer_etoh <- genotypes %>% 
  filter(stor == "EtOH") %>% 
  filter(n > 7) # take only the ones with at least 8 working

nrow(genotypes_fewer_etoh) # 63 samples with usable genotypes

# run allelematch
genotypes_allelematch_etoh <- genotypes_fewer_etoh %>% 
  select(lab_id, `1_TGLA94`:`2_CELB9`, -`SRY`) %>% # select all disomic loci and lab ID
  amDataset(indexColumn = "lab_id") %>%  # get into alleleMatch format
  amUnique(alleleMismatch = 4) # assign a unique individual ID to each of the samples based on genotype overlap

summary(genotypes_allelematch_etoh, csv="data/allelematch_etoh_output.csv")

# read back in (can't figure out better way?)
allelematch_output_etoh <- read.csv("data/allelematch_etoh_output.csv")

# how many unique individuals?
length(unique(allelematch_output_etoh$uniqueIndex)) # 44

# only swab ----------------------------------------------------------------

genotypes_fewer_swab <- genotypes %>% 
  filter(stor == "swab") %>% 
  filter(n > 7) # take only the ones with at least 8 working

nrow(genotypes_fewer_swab) # 91 samples with usable genotypes

# run allelematch
genotypes_allelematch_swab <- genotypes_fewer_swab %>% 
  select(lab_id, `1_TGLA94`:`2_CELB9`, -`SRY`) %>% # select all disomic loci and lab ID
  amDataset(indexColumn = "lab_id") %>%  # get into alleleMatch format
  amUnique(alleleMismatch = 4) # assign a unique individual ID to each of the samples based on genotype overlap

summary(genotypes_allelematch_swab, csv="data/allelematch_swab_output.csv")

# read back in (can't figure out better way?)
allelematch_output_swab <- read.csv("data/allelematch_swab_output.csv")

# how many unique individuals?
length(unique(allelematch_output_swab$uniqueIndex)) # 63
