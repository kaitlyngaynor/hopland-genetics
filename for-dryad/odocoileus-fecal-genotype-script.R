# Script for analyses in Bach et al. Mammalian Biology

# Data import and set-up -------------------------------------------------

library(lme4)
library(ggplot2)
library(dplyr)
library(rstatix)

odo_data <- read.csv("for-dryad/odocoileus-fecal-genotype-data.csv")
head(odo_data)

# create column corresponding to whether or not sample yielded 
# a working genotype (at least 8 alleles amplifying)

odo_data %>% 
  count(lab_id)




# Genotyping success analysis ---------------------------------------------






# qPCR analysis -----------------------------------------------------------

# Figure 4: DNA concentration of samples across storage methods (dry, ethanol, 
# and swab), as determined through qPCR
ggplot(odo_data, aes(x = log(qpcr_conc), col = storage)) +
  geom_density() +
  xlab("DNA concentration (log ng/µL)") +
  ylab("Density") +
  theme_classic() +
  theme(legend.title = element_blank())

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

