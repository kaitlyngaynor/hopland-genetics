### EXPLORE 2017 AND 2018 SWAB SAMPLES FOR MORE INFO


# Import and organize data ------------------------------------------------

# import raw genotype data for 2017 and 2018 samples
genotypes_17 <- read_csv("data/swab_transects/2017_genotypes.csv")
genotypes_18 <- read_csv("data/swab_transects/2018_genotypes.csv")

# each sample was run multiple times, AND there are consensus genotypes in there
# just take max n (alleles working) for a given lab_id (had some issues filtering only consensus, this is easier)
genotypes <- bind_rows(genotypes_17, genotypes_18) %>% 
  group_by(lab_id) %>% 
  slice(which.max(n)) %>% 
  select(lab_id, n) 

# add column for whether it is working
genotypes$working <- "no"
genotypes$working_01 <- 0
for(i in 1:nrow(genotypes)) {
  if(genotypes$n[i] >=8) {
    genotypes$working[i] <- "yes"
    genotypes$working_01[i] <- 1
  }
}

# import metadata
metadata_17 <- read_csv("data/swab_transects/2017_metadata.csv")
metadata_18 <- read_csv("data/swab_transects/2018_metadata.csv")

# just take the pre-fire 2018 samples (kind of a hack here)
metadata_18 <- metadata_18[1:746,]

# combine, just select lab ID and condition
metadata <- bind_rows(metadata_17, metadata_18) %>% 
  select(lab_id, cond_rank)

# join
genotypes <- left_join(genotypes, metadata)


# Fine-scale sample condition ---------------------------------------------

# summarize
(working.cond <- genotypes %>%
    drop_na(cond_rank) %>%
    group_by(as.factor(cond_rank)) %>%
    count(working))
names(working.cond) <- c("cond_rank", "working", "n")

ggplot(working.cond, aes(x = cond_rank, y = n, fill = working)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Percent of samples") +
  xlab("Sample condition") +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())


# Coarser sample condition (USE THIS) ---------------------------------------

# combine genotypes into fewer categories
genotypes$cond_coarse <- fct_collapse(as.factor(genotypes$cond_rank), 
                                      Slimy = c("1", "1.5"),
                                      Wet = c("2", "2.5"),
                                      Shiny = c("3", "3.5"),
                                      Dull = c("4", "4.5", "5"))

# summarize
(working.cond.coarse <- genotypes %>%
  drop_na(cond_coarse) %>%
  group_by(as.factor(cond_coarse)) %>%
  count(working))
names(working.cond.coarse) <- c("cond_coarse", "working", "n")

pdf(here::here("figures", "swab1718-condition-genotype-usable.pdf"), width = 4, height = 5)
ggplot(working.cond.coarse, aes(x = desc(cond_coarse), y = n, fill = working)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Percent of samples working") +
  xlab("Sample condition") +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
dev.off()

# summarize total number in each coarse condition (sample sizes)
working.cond.coarse %>% 
  group_by(cond_coarse) %>% 
  summarise(total = sum(n))


# GLM -------------------------------------------------------------------

# categorical
fit1 <- glm(working_01 ~ cond_coarse, data = genotypes, family = "binomial")
summary(fit1)

# continuous
fit2 <- glm(working_01 ~ cond_rank, data = genotypes, family = "binomial")
summary(fit2)
