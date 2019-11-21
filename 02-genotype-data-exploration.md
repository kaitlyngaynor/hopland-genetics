Set working directory and load data.

    library(here)
    library(tidyverse)
    library(scales)

    # import cleaned data file
    genotypes <- read_csv(here::here("data", "cleaned-triplicates-genotypes-qpcr.csv"))

Storage method comparison
-------------------------

Take a look at the effect of storage method on genotypes.

A histogram showing the number of samples for that amplified at each
number of alleles.

    # wrapped
    ggplot(genotypes, aes(x = n, fill = stor)) +
      geom_histogram() +
      facet_wrap(~ stor) +
      ylab("Number of Samples") +
      xlab("Number of Alleles Amplifying") +
      ggtitle("Breakdown of Markers Working Across Storage Methods")

![](02-genotype-data-exploration_files/figure-markdown_strict/storage%20and%20number%20worked-1.png)

    # side-by-side
    ggplot(genotypes, aes(x = n, fill = stor)) +
      geom_histogram(position = "dodge") +
      ylab("Number of Samples") +
      xlab("Number of Alleles Amplifying") +
      ggtitle("Breakdown of Markers Working Across Storage Methods")

![](02-genotype-data-exploration_files/figure-markdown_strict/storage%20and%20number%20worked-2.png)

Plot whether sample worked against storage method.

    # summarize
    working.storage <- genotypes %>%
      group_by(stor) %>%
      count(working)

    # histogram
    ggplot(working.storage, aes(x = working, y = n, fill = stor)) +
      geom_bar(stat = "identity", position = "dodge") +
      ggtitle("Samples Yielding Usable Genotypes (> 8/10 markers)") +
      ylab("Number of Samples")

![](02-genotype-data-exploration_files/figure-markdown_strict/storage%20method%20working-1.png)

    # another way of viewing it
    ggplot(working.storage, aes(x = stor, y = n, fill = working)) +
      geom_bar(stat = "identity") +
      ggtitle("Samples Yielding Usable Genotypes (> 8/10 markers)") +
      ylab("Number of Samples")

![](02-genotype-data-exploration_files/figure-markdown_strict/storage%20method%20working-2.png)

    # summarize
    perfect.storage <- genotypes %>%
      group_by(stor) %>%
      count(perfect)

    # histogram
    ggplot(perfect.storage, aes(x = perfect, y = n, fill = stor)) +
      geom_bar(stat = "identity", position = "dodge") +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers)") +
      ylab("Number of Samples")

![](02-genotype-data-exploration_files/figure-markdown_strict/storage%20method%20working-3.png)

    # another way of viewing it
    ggplot(perfect.storage, aes(x = stor, y = n, fill = perfect)) +
      geom_bar(stat = "identity") +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers)") +
      ylab("Number of Samples")

![](02-genotype-data-exploration_files/figure-markdown_strict/storage%20method%20working-4.png)

As per Bryan's 9/17/19 suggestion, a look at how individual samples
performed with different storage methods.

First, get data from long format to wide format.

    ## subset only columns consistent across samples (name, condition), then spread
    genotypes.spread <- genotypes %>%
      subset(select=c(lab_id, cond, stor, n)) %>% 
      spread(key = stor, value = n)  
      
    head(genotypes.spread) # confirm it worked, yay!

    ## # A tibble: 6 x 5
    ##   lab_id    cond   dry  EtOH  swab
    ##   <chr>    <dbl> <int> <int> <int>
    ## 1 odo_0124   3.5     1     1    10
    ## 2 odo_0125   3       7     0    10
    ## 3 odo_0126   3      10     3    10
    ## 4 odo_0127   3       0     6    10
    ## 5 odo_0128   3      10    10    10
    ## 6 odo_0129   3       3     0     8

Now plot. I never did figure out how to get a legend for the
geom\_dumbbell plot, so just added manually in Illustrator. This link
seems promising, though I didn't investigate fully:
<https://stackoverflow.com/questions/44653597/adding-a-traditional-legend-to-dumbbell-plot-in-ggaltgeom-dumbbell-in-r>

So, for now, just know that red = dry, and blue = swab.

    # need package ggalt for this plot
    library(ggalt)

    # get samples in order by their dry values
    genotypes.spread$lab_id <- factor(genotypes.spread$lab_id, levels = genotypes.spread$lab_id[order(genotypes.spread$dry)])

    ggplot(genotypes.spread, aes(x = dry, 
                                      xend = swab, 
                                      y = lab_id, 
                                      group = lab_id)) + 
      
      geom_dumbbell(color="gray", # color of bar
                    size=.5,
                    size_x=2,
                    size_xend=2,
                    colour_x = "red",
                    colour_xend="blue",
                    show.legend = T) + # color of end point (swab)
      labs(x=NULL, 
           y=NULL) +
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            plot.background=element_rect(fill="#f7f7f7"),
            panel.background=element_rect(fill="#f7f7f7"),
            panel.grid.minor=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.major.x=element_line(),
            axis.ticks=element_blank(),
            legend.position="top",
            panel.border=element_blank()) +
      xlab("Number of Markers Working") +
      ggtitle("Dry (red) vs Swab (blue)")

![](02-genotype-data-exploration_files/figure-markdown_strict/dumbbell%20plot-1.png)

In this case, green = ethanol.

    # get samples in order by their EtOH values
    genotypes.spread$lab_id <- factor(genotypes.spread$lab_id, levels = genotypes.spread$lab_id[order(genotypes.spread$EtOH)])

    ggplot(genotypes.spread, aes(x = EtOH, 
                                      xend = swab, 
                                      y = lab_id, 
                                      group = lab_id)) + 
      
      geom_dumbbell(color="gray", # color of bar
                    size=.5,
                    size_x=2,
                    size_xend=2,
                    colour_x = "darkgreen",
                    colour_xend="blue",
                    show.legend = T) + # color of end point (swab)
      labs(x=NULL, 
           y=NULL) +
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            plot.background=element_rect(fill="#f7f7f7"),
            panel.background=element_rect(fill="#f7f7f7"),
            panel.grid.minor=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.major.x=element_line(),
            axis.ticks=element_blank(),
            legend.position="top",
            panel.border=element_blank()) +
      xlab("Number of Markers Working") +
      ggtitle("Ethanol (green) vs Swab (blue)")

![](02-genotype-data-exploration_files/figure-markdown_strict/dumbbell%20plot%202-1.png)

Effect of condition on genotypes
--------------------------------

And now a look at the effect of condition.

    ggplot(genotypes, aes(x = cond, y = n, col = stor)) +
      geom_jitter() +
      geom_smooth() +
      facet_wrap(~ stor)

![](02-genotype-data-exploration_files/figure-markdown_strict/condition%20effect-1.png)

    ggplot(genotypes, aes(x = perfect, y = cond, fill = stor)) +
      geom_boxplot() +
    #  geom_jitter() +
      facet_wrap(~ stor) +
      coord_flip()

![](02-genotype-data-exploration_files/figure-markdown_strict/condition%20effect-2.png)

    # summarize
    working.storage.cond <- genotypes %>%
      group_by(stor, as.factor(cond)) %>%
      count(working)
    names(working.storage.cond) <- c("stor", "cond", "working", "n")

    # histogram
    ggplot(working.storage.cond, aes(x = cond, y = n, fill = working)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Usable Genotypes (> 8/10 markers) by Condition & Storage") +
      ylab("Number of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/condition%20effect-3.png)

    # scale to percentage so it's easier to compare across sample types
    ggplot(working.storage.cond, aes(x = cond, y = n, fill = working)) +
      geom_bar(stat = "identity", position = "fill") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Usable Genotypes (> 8/10 markers) by Condition & Storage") +
      ylab("Percent of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/condition%20effect-4.png)

    # summarize
    perfect.storage.cond <- genotypes %>%
      group_by(stor, as.factor(cond)) %>%
      count(perfect)
    names(perfect.storage.cond) <- c("stor", "cond", "perfect", "n")

    # histogram
    ggplot(perfect.storage.cond, aes(x = cond, y = n, fill = perfect)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers) by Condition & Storage") +
      ylab("Number of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/condition%20effect-5.png)

    # scale to percentage so it's easier to compare across sample types
    ggplot(perfect.storage.cond, aes(x = cond, y = n, fill = perfect)) +
      geom_bar(stat = "identity", position = "fill") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers) by Condition & Storage") +
      ylab("Percent of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/condition%20effect-6.png)

### Effect of condition - collapsed categories

Try combining 1 and 1.5, 2 and 2.5, etc to better facilitate comparisons
across conditions (since sample size is low for each condition at the
fine scale). Combining 4 and NA but I'm not sure how legit that is?

    genotypes$cond.coarse <- fct_collapse(as.factor(genotypes$cond), 
                                               Slimy = c("1", "1.5"),
                                               Wet = c("2", "2.5"),
                                               Shiny = c("3", "3.5"),
                                               Dull = c("4"))
    genotypes$cond.coarse <- fct_explicit_na(genotypes$cond.coarse, na_level = "Dull")

    # summarize
    working.storage.cond.coarse <- genotypes %>%
      group_by(stor, as.factor(cond.coarse)) %>%
      count(working)
    names(working.storage.cond.coarse) <- c("stor", "cond.coarse", "working", "n")

    # histogram
    ggplot(working.storage.cond.coarse, aes(x = cond.coarse, y = n, fill = working)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Usable Genotypes (> 8/10 markers) by Condition & Storage") +
      ylab("Number of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/combine%20condition%20scores-1.png)

    # scale to percentage so it's easier to compare across sample types
    ggplot(working.storage.cond.coarse, aes(x = cond.coarse, y = n, fill = working)) +
      geom_bar(stat = "identity", position = "fill") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Usable Genotypes (> 8/10 markers) by Condition & Storage") +
      ylab("Percent of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/combine%20condition%20scores-2.png)

    # summarize
    perfect.storage.cond.coarse <- genotypes %>%
      group_by(stor, cond.coarse) %>%
      count(perfect)
    names(perfect.storage.cond) <- c("stor", "cond.coarse", "perfect", "n")

    # histogram
    ggplot(perfect.storage.cond.coarse, aes(x = cond.coarse, y = n, fill = perfect)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers) by Condition & Storage") +
      ylab("Number of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/combine%20condition%20scores-3.png)

    # scale to percentage so it's easier to compare across sample types
    ggplot(perfect.storage.cond.coarse, aes(x = cond.coarse, y = n, fill = perfect)) +
      geom_bar(stat = "identity", position = "fill") +
      facet_wrap(~ stor) +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers) by Condition & Storage") +
      ylab("Percent of Samples") 

![](02-genotype-data-exploration_files/figure-markdown_strict/combine%20condition%20scores-4.png)

    # same thing, but organized differently
    ggplot(perfect.storage.cond.coarse, aes(x = stor, y = n, fill = perfect)) +
      geom_bar(stat = "identity", position = "fill") +
      facet_wrap(~ cond.coarse) +
      ggtitle("Samples Yielding Complete Genotypes (10/10 markers) by Condition & Storage") +
      ylab("Percent of Samples")

![](02-genotype-data-exploration_files/figure-markdown_strict/combine%20condition%20scores-5.png)
