Set working directory and load data.

    library(lme4)
    library(tidyverse)
    library(here)

    # import cleaned file
    genotypes <- read_csv(here::here("data", "cleaned-triplicates-genotypes-qpcr.csv"))

Generalized Linear Mixed Models
-------------------------------

We have a tricky modeling situation with the dependent variable—it is
effectively a proportion (number of markers working out of 10 markers),
and is bounded between 0 and 10 (or 0 and 1, if calculated as a true
proportion). Furthermore, it is not normally distributed, but instead 0-
and 1-inflated (most of the values are 0 or 1). When Googling, I ran
into these ZOIB models (zero-one inflated beta models) which are a
little confusing to me but seems promising??
<https://journal.r-project.org/archive/2015/RJ-2015-019/RJ-2015-019.pdf>
They entail Bayesian inference and are a little tricky. I didn't go down
this route.

We could instead treat the dependent variable as a binary variable
(worked vs not worked, with the "worked" cut-off being 8 or 10), which
simplifies our analysis greatly, and also may actually be more
meaningful, since we don't actually care about the difference between 4
worked and 5 worked. So this is the route that I took. I'm open to
feedback, and of course open to people playing around with the ZOIB
models on their own (though am not sure I have the capacity to open that
can of worms myself, ha).

I used sample (lab\_id) as a random effect, to control for the fact that
we tested the same samples in triplicate across methods. I tested both
storage method and condition as fixed effects, in addition to running a
null model (with only the random effect).

#### Working samples (&gt;= 8 markers)

    genotypes$working <- as.factor(genotypes$working)

    # no variables—just random effect
    fit0 <- glmer(working ~ (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit0) # AIC = 409.2

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: working ~ (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    405.7    413.4   -200.9    401.7      334 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.4985 -0.5768  0.3642  0.3642  1.0652 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 2.836    1.684   
    ## Number of obs: 336, groups:  lab_id, 112
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.0234     0.2437     4.2 2.67e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # effect of condition only
    fit1 <- glmer(working ~ cond + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit1) # AIC = 389.8

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: working ~ cond + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    386.5    397.9   -190.3    380.5      318 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6504 -0.6246  0.3807  0.3920  1.0449 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 2.351    1.533   
    ## Number of obs: 321, groups:  lab_id, 107
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)   1.6695     0.8413   1.984   0.0472 *
    ## cond         -0.2104     0.2899  -0.726   0.4679  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## cond -0.959

    # effect of condition only (coarse categories) - NEED TO UPDATE THIS
    #fit1.5 <- glmer(working ~ cond.coarse + (1 | lab_id), family = binomial("logit"), data = genotypes)
    #summary(fit1.5) # AIC = 404.2

    # effect of storage only
    fit2 <- glmer(working ~ stor + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit2) # AIC = 382.9

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: working ~ stor + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    378.9    394.2   -185.4    370.9      332 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1231 -0.5982  0.3202  0.4076  1.6718 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 4.708    2.17    
    ## Number of obs: 336, groups:  lab_id, 112
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.7546     0.3539   2.133    0.033 *  
    ## storEtOH     -0.2602     0.3617  -0.719    0.472    
    ## storswab      1.8280     0.4448   4.110 3.96e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) strEOH
    ## storEtOH -0.532       
    ## storswab -0.304  0.377

    # effect of storage and condition
    fit3 <- glmer(working ~ stor + cond + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit3) # AIC = 365.1

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: working ~ stor + cond + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    361.3    380.2   -175.6    351.3      316 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1202 -0.6238  0.3205  0.4271  1.6332 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 3.951    1.988   
    ## Number of obs: 321, groups:  lab_id, 107
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.5348     1.0384   1.478    0.139    
    ## storEtOH     -0.2617     0.3628  -0.721    0.471    
    ## storswab      1.8007     0.4508   3.994 6.49e-05 ***
    ## cond         -0.2523     0.3522  -0.716    0.474    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) strEOH strswb
    ## storEtOH -0.189              
    ## storswab -0.062  0.374       
    ## cond     -0.942  0.007 -0.040

    # effect of storage and condition (interacting)
    fit4 <- glmer(working ~ stor * cond + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit4) # AIC = 366

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: working ~ stor * cond + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    362.8    389.2   -174.4    348.8      314 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9792 -0.5603  0.3177  0.4108  1.5587 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 4.083    2.021   
    ## Number of obs: 321, groups:  lab_id, 107
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)     1.4164     1.3039   1.086   0.2774  
    ## storEtOH       -0.9517     1.4231  -0.669   0.5037  
    ## storswab        3.8846     1.9272   2.016   0.0438 *
    ## cond           -0.2055     0.4530  -0.454   0.6500  
    ## storEtOH:cond   0.2495     0.4981   0.501   0.6165  
    ## storswab:cond  -0.7239     0.6439  -1.124   0.2609  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) strEOH strswb cond   stEOH:
    ## storEtOH    -0.564                            
    ## storswab    -0.385  0.350                     
    ## cond        -0.963  0.544  0.390              
    ## strEtOH:cnd  0.544 -0.967 -0.342 -0.562       
    ## storswb:cnd  0.402 -0.361 -0.971 -0.428  0.376

#### Perfect samples (all 10 markers)

    genotypes$perfect <- as.factor(genotypes$perfect)

    # no variables—just random effect
    fit0 <- glmer(perfect ~ (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit0) # AIC = 428.5

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: perfect ~ (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    425.7    433.3   -210.8    421.7      334 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3820 -0.5197  0.4064  0.4064  1.1531 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 2.925    1.71    
    ## Number of obs: 336, groups:  lab_id, 112
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)   0.5566     0.2171   2.564   0.0104 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # effect of condition only
    fit1 <- glmer(perfect ~ cond + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit1) # AIC = 409.6

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: perfect ~ cond + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    407.1    418.4   -200.5    401.1      318 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5536 -0.8677  0.4188  0.4346  1.1525 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 2.549    1.597   
    ## Number of obs: 321, groups:  lab_id, 107
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)   1.4855     0.8461   1.756   0.0791 .
    ## cond         -0.2955     0.2938  -1.006   0.3145  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## cond -0.964

    # effect of condition only (coarse categories)
    #fit1.5 <- glmer(perfect ~ cond.coarse + (1 | lab_id), family = binomial("logit"), data = genotypes)
    #summary(fit1.5) # AIC = 421.8

    # effect of storage only
    fit2 <- glmer(perfect ~ stor + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit2) # AIC = 404.1

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: perfect ~ stor + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    400.8    416.1   -196.4    392.8      332 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7013 -0.5724  0.1761  0.4594  1.7469 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 4.611    2.147   
    ## Number of obs: 336, groups:  lab_id, 112
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.2221     0.3339   0.665    0.506    
    ## storEtOH     -0.2559     0.3586  -0.714    0.475    
    ## storswab      1.6619     0.4144   4.010 6.07e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) strEOH
    ## storEtOH -0.542       
    ## storswab -0.425  0.403

    # effect of storage and condition
    fit3 <- glmer(perfect ~ stor + cond + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit3) # AIC = 385.2

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: perfect ~ stor + cond + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    382.2    401.1   -186.1    372.2      316 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7480 -0.5729  0.1787  0.4812  1.7901 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 4.116    2.029   
    ## Number of obs: 321, groups:  lab_id, 107
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.3597     1.0344   1.314    0.189    
    ## storEtOH     -0.3249     0.3621  -0.897    0.370    
    ## storswab      1.6567     0.4233   3.913  9.1e-05 ***
    ## cond         -0.3543     0.3538  -1.002    0.317    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) strEOH strswb
    ## storEtOH -0.190              
    ## storswab -0.078  0.390       
    ## cond     -0.947  0.012 -0.055

    # effect of storage and condition (interacting)
    fit4 <- glmer(perfect ~ stor * cond + (1 | lab_id), family = binomial("logit"), data = genotypes)
    summary(fit4) # AIC = 386.4

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: perfect ~ stor * cond + (1 | lab_id)
    ##    Data: genotypes
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    383.8    410.2   -184.9    369.8      314 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6083 -0.5909  0.1873  0.4667  1.8807 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  lab_id (Intercept) 4.198    2.049   
    ## Number of obs: 321, groups:  lab_id, 107
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)     1.0423     1.2945   0.805    0.421  
    ## storEtOH       -0.7429     1.4125  -0.526    0.599  
    ## storswab        3.6698     1.7410   2.108    0.035 *
    ## cond           -0.2374     0.4531  -0.524    0.600  
    ## storEtOH:cond   0.1517     0.4963   0.306    0.760  
    ## storswab:cond  -0.7105     0.5900  -1.204    0.229  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) strEOH strswb cond   stEOH:
    ## storEtOH    -0.552                            
    ## storswab    -0.434  0.386                     
    ## cond        -0.966  0.534  0.427              
    ## strEtOH:cnd  0.532 -0.966 -0.377 -0.551       
    ## storswb:cnd  0.443 -0.392 -0.969 -0.464  0.408

The best models (lowest AIC) for both working genotypes (with at least 8
amplifying markers) and perfect genotypes (all 10 markers working)
include both storage and condition (but NOT the interaction between the
two), as well as sample ID as a random effect. This suggests that both
storage method and sample condition have important (but independent)
effects on whether or not a sample works.

I tested two versions of condition score: one numerical (with 0.5
intervals), and one as a factor (coarser, with four categories), and the
numerical one was better in both models.
