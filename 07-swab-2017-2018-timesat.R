### EXPLORE 2017 AND 2018 SWAB SAMPLES FOR MORE INFO

library(tidyverse)

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


# Format 2017 metadata ----------------------------------------------------

metadata_17 <- read.csv("data/swab_transects/2017_metadata.csv") %>%
  select(lab_id, date_collected, cond_rank)

# import 2017 extraction dates and join with collection dates
extraction_17 <- read_csv("data/swab_transects/2017_extraction.csv")
extraction_17$lab_id <- substr(extraction_17$extraction_id, start = 1, stop = 8) # generate new column for lab ID based on extraction ID

metadata_17 <- left_join(metadata_17, extraction_17) %>% 
  drop_na() %>% 
  select(lab_id, date_collected, ext_date, cond_rank) %>% 
  mutate(date_collected = as.Date(date_collected, format = "%m/%d/%y"),
         ext_date = as.Date(ext_date, format = "%m/%d/%y"))


# Format 2018 metadata ----------------------------------------------------

metadata_18 <- read.csv("data/swab_transects/2018_metadata.csv") %>% 
  select(lab_id, date_collected, ext_date, cond_rank) %>% 
  mutate(date_collected = as.Date(date_collected, format = "%m/%d/%y"),
         ext_date = as.Date(ext_date, format = "%Y-%m-%d"))

# just take the pre-fire 2018 samples (kind of a hack here)
metadata_18 <- metadata_18[1:746,]


# Calculate time sat ------------------------------------------------------

# combine, just select lab ID and dates
metadata <- bind_rows(metadata_17, metadata_18) 

# join and drop NA (will get rid of post fire)
genotypes_dates <- left_join(genotypes, metadata) %>% 
  drop_na()

# calculate difference in dates
genotypes_dates$time_sat <- as.numeric(genotypes_dates$ext_date - genotypes_dates$date_collected)

hist(genotypes_dates$time_sat)

ggplot(genotypes_dates, aes(x = working, y = time_sat)) +
  geom_boxplot()

# combine genotypes into fewer categories
genotypes_dates$cond_coarse <- fct_collapse(as.factor(genotypes_dates$cond_rank), 
                                      Slimy = c("1", "1.5"),
                                      Wet = c("2", "2.5"),
                                      Shiny = c("3", "3.5"),
                                      Dull = c("4", "4.5", "5"))


# GLM -------------------------------------------------------------------

# models listed from worst to best

# null - AIC = 2835.7
null_model <- glm(working_01 ~ 1, data = genotypes_dates, family = "binomial")
summary(null_model)

# time sat - AIC = 2830.1
timesat_model <- glm(working_01 ~ time_sat, data = genotypes_dates, family = "binomial")
summary(timesat_model)

# condition (continuous) + time sat - AIC = 2674.4
timesat_cond_model4 <- glm(working_01 ~ time_sat + cond_rank, data = genotypes_dates, family = "binomial")
summary(timesat_cond_model4)

# condition (continuous) - AIC = 2709.1
cond_model2 <- glm(working_01 ~ cond_rank, data = genotypes_dates, family = "binomial")
summary(cond_model2)

# condition (continuous) * time sat - AIC = 2656.7 - BEST
timesat_cond_model2 <- glm(working_01 ~ time_sat * cond_rank, data = genotypes_dates, family = "binomial")
summary(timesat_cond_model2)

# don't use these - with categorical condition (does not align with previous methods)

## condition (categorical) - AIC = 2693.3
#cond_model1 <- glm(working_01 ~ cond_coarse, data = genotypes_dates, family = "binomial")
#summary(cond_model1)
#
## condition (categorical) + time sat - AIC = 2668.8
#timesat_cond_model3 <- glm(working_01 ~ time_sat + cond_coarse, data = genotypes_dates, family = "binomial")
#summary(timesat_cond_model3)
#
## condition (categorical) * time sat - AIC = 2647.4
#timesat_cond_model1 <- glm(working_01 ~ time_sat * cond_coarse, data = genotypes_dates, family = "binomial")
#summary(timesat_cond_model1)


# Figures -----------------------------------------------------------------

pdf(here::here("figures", "swab-storage-time.pdf"), width = 4, height = 3)
ggplot(data=genotypes_dates, mapping=aes(x=time_sat,y=working_01)) + 
  geom_point(size = .5, position = position_jitter(.5,0.05)) +
  stat_smooth(method="lm", method.args=list(family=binomial)) +
  theme_classic() +
  xlab("Time in Storage (Days)") +
  ylab("Probability of Success")
dev.off()

# what is the slope of this line?
summary(timesat_model) # from above
dummy_data <- data.frame(time_sat = seq(0, max(genotypes_dates$time_sat, na.rm=T), by = 0.001)) 
dummy_predictions <- predict(timesat_model, newdata = dummy_data, type = "response")
summary(dummy_predictions)

# change per day
(max(dummy_predictions) - min(dummy_predictions))/max(genotypes_dates$time_sat)

# by condition
ggplot(data=genotypes_dates, mapping=aes(x=time_sat,y=working_01, col = cond_coarse)) + 
  geom_point(size = .5, position = position_jitter(.5,0.05)) +
  stat_smooth(method="lm", method.args=list(family=binomial)) +
  theme_classic() +
  xlab("Time in Storage (Days)") +
  ylab("Probability of Success")
