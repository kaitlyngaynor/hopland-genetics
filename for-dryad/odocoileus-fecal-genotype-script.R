# Script for analyses in Bach et al. Mammalian Biology

library(lme4)
library(ggplot2)
library(dplyr)
library(rstatix)
library(AICcmodavg)


# Sample summary ---------------------------------------------------------

odo_data <- read.csv("odocoileus-fecal-genotype-data.csv")
odo_data$working <- as.factor(odo_data$working)
head(odo_data)

# subset data frame so there is one row per unique pellet pile (collected in triplicate)
samples <- odo_data %>% 
  select(lab_id, condition, condition_coarse, individual_id) %>% 
  unique

# summarize sample conditions
mean(samples$condition, na.rm = T)
sd(samples$condition, na.rm = T)
count(samples, condition_coarse)

# summarize number that are identifiable
nrow(filter(samples, is.na(individual_id) == FALSE)) # 99

# summarize number of unique deer
length(unique(samples$individual_id)) # 69

# summarize number of samples per deer
samples %>% 
  count(individual_id) %>% 
  arrange(n)



# Genotyping success analysis ---------------------------------------------

# alleles amplified across storage methods
odo_data %>% 
  group_by(storage) %>% 
  summarise(mean_n = mean(alleles_amplified, na.rm = TRUE),
            sd_n = sd(alleles_amplified, na.rm = TRUE))

# usable samples across storage methods
odo_data %>% 
  group_by(storage, working) %>% 
  count %>% 
  pivot_wider(names_from = working, values_from = n) %>% 
  mutate(percent_working = yes / (yes+no))

# GLMMs of predictors of usable genotypes

# no variables—just random effect
fit0 <- glmer(working ~ (1 | lab_id), family = binomial("logit"), data = odo_data)
AICc(fit0) # AICc = 405.8

# effect of condition only
fit1 <- glmer(working ~ condition + (1 | lab_id), family = binomial("logit"), data = odo_data)
AICc(fit1) # AICc = 386.6

# effect of condition only (factor, coarse categories) 
fit1.5 <- glmer(working ~ condition_coarse + (1 | lab_id), family = binomial("logit"), data = odo_data)
AICc(fit1.5) # AICc = 387.1

# effect of storageage only
fit2 <- glmer(working ~ storage + (1 | lab_id), family = binomial("logit"), data = odo_data)
AICc(fit2) # AICc = 379.0

# effect of storageage and condition - BEST!!!
fit3 <- glmer(working ~ storage + condition + (1 | lab_id), family = binomial("logit"), data = odo_data)
AICc(fit3) # AICc = 361.5
summary(fit3)

# effect of storageage and condition (interacting)
fit4 <- glmer(working ~ storage * condition + (1 | lab_id), family = binomial("logit"), data = odo_data)
AICc(fit4) # AICc = 363.1



# qPCR analysis -----------------------------------------------------------

# calculate summary statistics
odo_data %>%
   group_by(storage) %>%
   summarize(mean = mean(qpcr_conc, na.rm = TRUE),
             median = median(qpcr_conc, na.rm = TRUE),
             sd = sd(qpcr_conc, na.rm = TRUE),
             min = min(qpcr_conc, na.rm = TRUE),
             max = max(qpcr_conc, na.rm = TRUE),
             quantile25 = quantile(qpcr_conc, 0.25, na.rm = TRUE),
             quantile75 = quantile(qpcr_conc, 0.75, na.rm = TRUE))

# run repeated-measures ANOVA
odo_data %>% 
  select(lab_id, storage, qpcr_conc) %>% 
  drop_na() %>% 
  mutate(qpcr_conc_log = log(qpcr_conc)) %>% # log-transform data
  anova_test(dv = qpcr_conc_log,
             wid = lab_id,
             between = storage)

# explore relationship between qPCR analysis and genotyping success

# dry
odo_data_dry_work <- subset(odo_data, storage == "dry" & working == "yes")
odo_data_dry_workno <- subset(odo_data, storage == "dry" & working == "no")
t.test(log(odo_data_dry_work$qpcr_conc), log(odo_data_dry_workno$qpcr_conc))
# ethanol
odo_data_EtOH_work <- subset(odo_data, storage == "EtOH" & working == "yes")
odo_data_EtOH_workno <- subset(odo_data, storage == "EtOH" & working == "no")
t.test(log(odo_data_EtOH_work$qpcr_conc), log(odo_data_EtOH_workno$qpcr_conc))
# swab
odo_data_swab_work <- subset(odo_data, storage == "swab" & working == "yes")
odo_data_swab_workno <- subset(odo_data, storage == "swab" & working == "no")
t.test(log(odo_data_swab_work$qpcr_conc), log(odo_data_swab_workno$qpcr_conc))



# Modeling with 2017-2018 swab samples ------------------------------------

odo_data_1718swab <- read.csv("odocoileus-fecal-genotype-data-20172018swab.csv")
odo_data_1718swab$working <- as.factor(odo_data_1718swab$working)

# null - AIC = 2835.7
null_model <- glm(working ~ 1, data = odo_data_1718swab, family = "binomial")
summary(null_model)

# time in storage - AIC = 2830.1
timesat_model <- glm(working ~ days_storage, data = odo_data_1718swab, family = "binomial")
summary(timesat_model)

# condition + time in storage - AIC = 2674.4
timesat_cond_model4 <- glm(working ~ days_storage + condition, data = odo_data_1718swab, family = "binomial")
summary(timesat_cond_model4)

# condition - AIC = 2709.1
cond_model2 <- glm(working ~ condition, data = odo_data_1718swab, family = "binomial")
summary(cond_model2)

# condition * time in storage - AIC = 2656.7 - BEST
timesat_cond_model2 <- glm(working ~ days_storage * condition, data = odo_data_1718swab, family = "binomial")
summary(timesat_cond_model2)
