
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesnetworks

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/USCbiostats/bayesnetworks.svg?branch=master)](https://travis-ci.org/USCbiostats/bayesnetworks)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/bayesnetworks)](https://cran.r-project.org/package=bayesnetworks)
<!-- badges: end -->

The goal of bayesnetworks is to do MCMC fitting of Bayesian networks,
including source & sink constraints and external network info.

## Installation

<!--
You can install the released version of bayesnetworks from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bayesnetworks")
```
-->

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/bayesnetworks")
```

## Examples

The network is fitted by passing in data and network structure

``` r
library(bayesnetworks)
set.seed(1234)

x <- bn_mcmc(X = network$data, graph = network$dag_info, N = 50000)
```

Some diagnostic plots:

``` r
library(ggplot2)

ggplot(x, aes(iter, globalLL)) +
  geom_line()
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
ggplot(x, aes(iter, FN)) +
  geom_line(aes(color = "FN")) +
  geom_line(aes(y = FP, color = "FP")) +
  labs(y = "value")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
ggplot(x, aes(iter, deletions)) +
  geom_line(aes(color = "deletions")) +
  geom_line(aes(y = additions, color = "additions")) +
  labs(y = "count")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Code of Conduct

Please note that the ‘bayesnetworks’ project is released with a
[Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By
contributing to this project, you agree to abide by its terms.
