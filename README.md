cdsrplots
================

cdsrplots is an R toolkit for common CDS plotting functions

## Install

``` r
library(devtools)
devtools::install_github("broadinstitute/cdsr_plots")
```

The package can then be loaded by calling

``` r
library(cdsrplots)
```

## Ploting function

make_volcano is carried over from the cdsr package and can be used to plot the output of cdsrmodels::lin_associations() or cdsrmodels::run_lm_stats_limma()

Examples:

(1) 

lin_associations_output <- cdsrmodels::lin_associations(X, Y)

lin_output <- lin_associations_output$res.table %>% tibble::rownames_to_column('Gene')

cdsrplots::make_volcano(lin_output, 'PosteriorMean', 'qvalue', label_var = 'Gene')

(2)

lm_stats_output <- cdsrmodels::run_lm_stats_limma(mat, vec)

cdsrplots::make_volcano(lm_stats_output, 'EffectSize', 'p.value', q_var = 'q.value', label_var = 'Gene')

