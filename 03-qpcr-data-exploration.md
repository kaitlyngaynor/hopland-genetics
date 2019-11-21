Set working directory and load data.

    library(here)
    library(tidyverse)

    # import cleaned data file
    genotypes <- read_csv(here::here("data", "cleaned-triplicates-genotypes-qpcr.csv"))

qPCR by storage method
----------------------

    ggplot(genotypes, aes(x = stor, y = extracted_conc, fill = stor)) +
      geom_boxplot() +
      ylab("DNA Concentration") +
      xlab("Storage Method") +
      ggtitle("DNA Concentration Across Storage Methods")

    ## Warning: Removed 53 rows containing non-finite values (stat_boxplot).

![](03-qpcr-data-exploration_files/figure-markdown_strict/concentration%20by%20storage-1.png)

    ggplot(genotypes, aes(x = stor, y = log(extracted_conc), fill = stor)) +
      geom_boxplot() +
      ylab("DNA Concentration (log)") +
      xlab("Storage Method") +
      ggtitle("DNA Concentration Across Storage Methods")

    ## Warning: Removed 53 rows containing non-finite values (stat_boxplot).

![](03-qpcr-data-exploration_files/figure-markdown_strict/concentration%20by%20storage-2.png)

qPCR vs condition, by storage method
------------------------------------

    ggplot(genotypes, aes(x = cond, y = log(extracted_conc), col = stor)) +
      geom_jitter() +
      facet_wrap(~ stor) +
      ylab("DNA Concentration (log)") +
      xlab("Sample Condition") 

    ## Warning: Removed 66 rows containing missing values (geom_point).

![](03-qpcr-data-exploration_files/figure-markdown_strict/concentration%20by%20storage%20*%20condition-1.png)

Calculate some summary stats

    qpcr.summary <- genotypes %>%
      group_by(stor) %>%
      summarize(mean = mean(extracted_conc, na.rm = TRUE),
                median = median(extracted_conc, na.rm = TRUE),
                sd = sd(extracted_conc, na.rm = TRUE),
                min = min(extracted_conc, na.rm = TRUE),
                max = max(extracted_conc, na.rm = TRUE)
                )

    qpcr.summary

    ## # A tibble: 3 x 6
    ##   stor   mean median    sd       min   max
    ##   <chr> <dbl>  <dbl> <dbl>     <dbl> <dbl>
    ## 1 dry   18.7    4.73 42.9  0.0101    300. 
    ## 2 EtOH   7.09   3.39  9.80 0.0000709  54.3
    ## 3 swab  15.4    1.94 66.4  0.0000709 574.

Quick ANOVA to see if there is a difference in DNA concentration across
storage methods (there isn't).

    fit <- aov(extracted_conc ~ stor, data = genotypes)
    summary(fit)

    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## stor          2   6706    3353    1.58  0.208
    ## Residuals   280 594184    2122               
    ## 53 observations deleted due to missingness

Relationship between genotyping and qPCR
----------------------------------------

Here, plotting the number of amplifying genotypes vs. DNA concentration
(qPCR). It doesn't look like there is much of a relationship, even
though we might expect a positive relationship.

    ggplot(genotypes, aes(log(extracted_conc), n, col = stor)) +
      geom_jitter()

    ## Warning: Removed 53 rows containing missing values (geom_point).

![](03-qpcr-data-exploration_files/figure-markdown_strict/genotyping%20vs%20qpcr-1.png)
