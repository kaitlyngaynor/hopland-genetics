# Hopland genetics manuscript
Code and data for manuscript comparing DNA storage methods.

**01-data-cleaning.R**: Takes genotype results and qPRC results, merges together, restricts to only triplicate samples, and adds columns for whether each sample was "working" (>=8 alleles) or "perfect." Exports file for later analysis.

**02-genotype-data-exploration.Rmd**: Exploratory figures of genotyping results.

**03-qpcr-data-exploration.Rmd**: Exploratory figures of qPCR results, and qPCR * genotyping results.

**04-data-analysis.Rmd**: Analysis of the effects of storage method and sample condition on genotyping and qPCR results.

**05-figures.R**: Make final figures for publication.
