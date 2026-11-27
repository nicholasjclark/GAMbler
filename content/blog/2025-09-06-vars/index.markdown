---
title: "State-Space Vector Autoregressions in mvgam"
author: "Nicholas Clark"
subtitle: ""
excerpt: "Vector Autoregressions (VAR models), also known as Multivariate Autoregressions (MAR models), offer a way to model delayed and contemporaneous interactions among sets of multiple time series. These models are widely used in econometrics and psychology, among other fields, where they can be analyzed to ask many interesting questions about potential causality or stability. But software to fit these models to real-world time series, which often present as non-Gaussian counts, proportions or even binary observations with measurement error, is lacking. Here I show how to fit VARs in a State-Space format, and how to interrogate the models to ask meaningful questions about interactions and stability, using the mvgam package in R."
slug: vector-autoregressions
date: 2024-09-06
draft: false
images:
series:
tags:
  - rstats
  - mvgam
  - tutorial
  - forecasting
  - time-series
categories:
  - rstats
  - time-series
  - mvgam
layout: single-sidebar
---


*Updated in November 2026: the code now needs `mvgam` 2.0 or higher. The [package website](https://nicholasjclark.github.io/mvgam/) explains how to install it and what changed.*

[Vector Autoregressions (VAR models)](https://otexts.com/fpp3/VAR.html), also known in the ecological literature as [Multivariate Autoregressions (MAR models)](https://atsa-es.github.io/MARSS/), are a class of Dynamic Linear Models that offer a principled way to model both delayed and contemporaneous interactions among sets of multiple time series. These models are widely used in econometrics and psychology, among other fields, where they can be analyzed to ask many interesting questions about [potential causality](https://www.jstor.org/stable/2951647), [cascading temporal interaction effects](https://www.r-econometrics.com/timeseries/irf/) or [community stability](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/0012-9615%282003%29073%5B0301%3AECSAEI%5D2.0.CO%3B2). A VAR(1) model, where temporal interactions among series can take place at a lag of one time step, is defined as:


`$$y_t \sim{MVNormal}(\alpha + A * y_{t-1}, \Sigma)$$`

Where: 
- `\(y_t\)` is a vector of time series observations at time `\(t\)`
- `\(\alpha\)` is a vector of constants capturing the long-term means of each series in `\(y\)`
- `\(A\)` is a matrix of autoregressive coefficients estimating lagged dependence and cross-dependence between the elements of `\(y\)`
- `\(\Sigma\)` determines the spread (or flexibility) of the process and any contemporaneous correlations among time series errors


Many R packages exist to estimate parameters of VARs, including the [`vars`](https://cran.r-project.org/web/packages/vars/index.html) and [`bsvars`](https://bsvars.github.io/bsvars/) packages. However, as we can see from the model description above, the traditional VAR model makes two big assumptions: first, this model assumes that the time series being analyzed do not show any measurement error. This assumption may hold for certain types of time series that can be nearly perfectly measured, such as sales of items in a store with accurate digital sale tracking. But it is likely to be broken for the majority of series that are measured with imperfect detectors or with time-varying effort. This can be resolved if we instead refer to a [State-Space](https://www.cs.unc.edu/~welch/kalman/media/pdf/Kalman1960.pdf) representation of the model:

<br>

![Illustration of a basic State-Space model, which assumes that  a latent dynamic *process* (X) can evolve independently from the way we take *observations* (Y) of that process](https://raw.githubusercontent.com/nicholasjclark/mvgam/refs/heads/master/vignettes/SS_model.svg)

<br>

A State-Space model is convenient for dealing with many types of real-world time series because they can readily handle observation error, missing values and a broad range of dependency structures. For example, we can modify the VAR(1) described above to a State-Space representation using:

`$$y_t \sim{MVNormal}(\alpha + x_t, \Sigma_{obs})$$`
`$$x_t \sim{MVNormal}(A * x_{t-1}, \Sigma_{process})$$`
Where the observation time series `\(y_t\)` are now considered to be noisy observations of some latent `\(process\)` (labelled as `\(x_t\)`) and it is this `\(process\)` that evolves as a Vector Autoregression. The advantage of this model is that we can separately estimate process and observation errors (modelled using the two `\(\Sigma\)` matrices). This type of formulation is available using the [`MARSS`](https://atsa-es.github.io/MARSS/) package, which has a wonderfully detailed [list of worked examples](https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf) to showcase how flexible and useful these models are. 

But the above approach does not address the second major limitation of traditional VARS, which is that they assume that the observation series can be captured with a Gaussian observation model. This makes it difficult to fit VARs to many real-world time series such as counts of multiple species over time, time series of non-negative rates or proportions of customers that responded to a set of surveys. It turns out that this limitation is a much bigger issue when it comes to fitting VARs to real-world data. At the time of writing I could not find a single R package that makes these models accessible, apart from my [`mvgam`](https://nicholasjclark.github.io/mvgam/) package. Given that the question of how to fit time series models to non-Gaussian observations comes up [again](https://stats.stackexchange.com/questions/161617/software-or-workarounds-for-vector-autoregression-on-count-data) and [again](https://stats.stackexchange.com/questions/47204/how-to-regress-a-time-series-of-proportions) and [again](https://stats.stackexchange.com/questions/244660/modeling-multivariate-time-series-count-data-in-r) (and [again](https://stats.stackexchange.com/questions/70182/time-series-dynamic-poisson-regression)), often with unsatisfying answers, this post will demonstrate how the [`mvgam`](https://nicholasjclark.github.io/mvgam/) package can fill this gap. Herein I will show how the package can be used to fit complex multivariate autoregressive models to non-Gaussian time series and how to interrogate the models to ask insightful questions about temporal interactions, variance decompositions and community stability.

## Environment setup
This tutorial relies on the following packages:



``` r
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(tidyverse)       # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(tidybayes)       # Graceful plotting of Bayesian posterior estimates
library(farver)          # Colour space manipulations
```

I also define a few plotting utilities to streamline my analysis workflow

``` r
theme_set(theme_classic(base_size = 15,
                        base_family = 'serif'))
myhist = function(...){
  geom_histogram(col = 'white',
                 fill = '#B97C7C', ...)
}
hist_theme = function(){
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
}
```

## Load, manipulate and plot the data
We will work with a set of annual American kestrel (*Falco sparverius*) abundance time series, which represent adjusted annual counts of this species taken in three adjacent regions of British Columbia, Canada. These data were collected annually, corrected for changes in observer coverage and detectability, and logged. They can be accessed in the [`MARSS`](https://atsa-es.github.io/MARSS/) package

``` r
load(url('https://github.com/atsa-es/MARSS/raw/master/data/kestrel.rda'))
head(kestrel)
```

```
##      Year British.Columbia Alberta Saskatchewan
## [1,] 1969            0.754   0.460        0.000
## [2,] 1970            0.673   0.899        0.192
## [3,] 1971            0.734   1.133        0.280
## [4,] 1972            0.589   0.528        0.386
## [5,] 1973            1.405   0.789        0.451
## [6,] 1974            0.624   0.528        0.234
```

Arrange the data into a `data.frame` that spreads the time series observations into a 'long' format, [which is required for `mvgam` modelling](https://nicholasjclark.github.io/mvgam/articles/data_in_mvgam.html)

``` r
regions <- c("BC",
             "Alb",
             "Sask")
model_data <- do.call(rbind,
                      lapply(seq_along(regions),
                             function(x){
  data.frame(year = kestrel[, 1],
             # Reverse the logging so that we deal directly 
             # with the detection-adjusted counts
             adj_count = exp(kestrel[, 1 + x]),
             region = regions[x])})) %>%
  # Add series and time indicators for mvgam modelling
  dplyr::mutate(series = as.factor(region),
                time = year)
```

Inspect the time series data structure

``` r
head(model_data)
```

```
##   year adj_count region series time
## 1 1969  2.125485     BC     BC 1969
## 2 1970  1.960109     BC     BC 1970
## 3 1971  2.083398     BC     BC 1971
## 4 1972  1.802185     BC     BC 1972
## 5 1973  4.075527     BC     BC 1973
## 6 1974  1.866379     BC     BC 1974
```

``` r
dplyr::glimpse(model_data)
```

```
## Rows: 120
## Columns: 5
## $ year      <dbl> 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, …
## $ adj_count <dbl> 2.125485, 1.960109, 2.083398, 1.802185, 4.075527, 1.866379, …
## $ region    <chr> "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", …
## $ series    <fct> BC, BC, BC, BC, BC, BC, BC, BC, BC, BC, BC, BC, BC, BC, BC, …
## $ time      <dbl> 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, …
```

``` r
levels(model_data$series)
```

```
## [1] "Alb"  "BC"   "Sask"
```

Plot all three time series together using `mvgam`'s `mvgam_data()`, which also checks the data against the family we plan to use

``` r
mvgam_data(model_data,
           y = 'adj_count',
           family = Gamma(link = 'log'))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" alt="American kestrel adjusted count time series in R" width="60%" style="display: block; margin: auto;" />

Now plot some additional features, including the empirical CDF and estimated autocorrelation functions, for just one series at a time

``` r
mvgam_data(model_data,
           y = 'adj_count',
           family = Gamma(link = 'log'),
           series = 'Alb')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" alt="American kestrel time series features in R" width="60%" style="display: block; margin: auto;" />



``` r
mvgam_data(model_data,
           y = 'adj_count',
           family = Gamma(link = 'log'),
           series = 'BC')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" alt="American kestrel time series features in R" width="60%" style="display: block; margin: auto;" />

These latter plots make it immediately apparent that our outcome variable, i.e. the adjusted counts of kestrels in each region over time, only takes on non-negative real numbers:

``` r
ggplot(model_data,
       aes(x = adj_count)) +
  myhist() +
  labs(x = 'Adjusted count', y = '') +
  hist_theme()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" alt="Empirical distribution of American kestrel adjusted counts" width="60%" style="display: block; margin: auto;" />

## Build a State-Space VAR(1) model
Fitting a multivariate time series model to these adjusted counts is not easy using existing R packages. We'd almost certainly have to log-transform the outcome and hope that we don't lose too much inferential ability in doing so. But, [as I've shown in other blogposts on this site](https://ecogambler.netlify.app/blog/autocorrelated-gams/), there is a better and simpler way to respect the nature of these observations using `mvgam`. Here we will consider a State-Space model that uses a latent Vector Autoregression as a VAR(1) process. This is done by taking advantage of two key features of the workhorse `mvgam()` function: First, we can use the `trend_formula` argument to specify that we want a State-Space model. And second, we can add a `VAR()` term to that formula to specify that we want the latent process to evolve as a Vector Autoregression of order 1. The [`mvgam` cheatsheet](https://nicholasjclark.github.io/mvgam/cheatsheets/mvgam_cheatsheet.pdf) goes into some detail about how each of these arguments work, but your best bet is to read the extensive documentation in `?mvgam::mvgam` to fully understand how these models can be constructed. `mvgam` also takes advantage of [a major development in Bayesian VAR modelling, originally published by Sarah Heaps, that enforces stationarity through a principled prior on the autoregressive coefficients](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2079648). This is extremely helpful for ensuring that we find well-behaving models that don't give nonsensical and unwieldy forecasts where the forecast variance grows without bound. Given that the average counts of these series seems to vary across regions, we will also want to include region-level intercepts in the process model. And finally, given the nature of the response variable (non-negative real values), a suitable observation family would be the Gamma distribution with a log link. First we should inspect default prior distributions for the stochastic parameters involved in our model:

``` r
def_priors <- get_prior(adj_count ~ 0,
                        trend_formula = ~ region + VAR(cor = TRUE),
                        data = model_data,
                        family = Gamma(link = 'log'))
def_priors[, c('prior', 'class', 'coef')]
```

```
##                   prior         class       coef    source
##       gamma(0.01, 0.01)         shape            (unknown)
##      inv_gamma(3, 1.07)   sigma_trend            (unknown)
##  normal(0, sqrt(0.455))     Amu_trend            (unknown)
##  gamma(1.365, 0.071175)  Aomega_trend            (unknown)
##    lkj_corr_cholesky(2) L_Omega_trend            (unknown)
##                  (flat)       b_trend            (unknown)
##                  (flat)       b_trend  regionAlb (unknown)
##                  (flat)       b_trend   regionBC (unknown)
##                  (flat)       b_trend regionSask (unknown)
```

The regional coefficients in the process model (class `b_trend`) have flat priors by default. The process standard deviations (`sigma_trend`, the square roots of the diagonals of the process error matrix `\(\Sigma_{process}\)`) have an inverse Gamma prior. It is generally good practice to update these using domain knowledge if you want a model that is well-behaved and computationally efficient. Given that the adjusted counts are fairly small and that we are using a log link, we don't expect massive effect sizes in the linear predictor so a more suitable prior for the regression coefficients would be Standard Normal (i.e. `\(\beta \sim Normal(0, 1)\)`). We can also change the prior on the process standard deviations to an Exponential prior that places less belief on extremely large values. We will leave the priors for the Gamma shape parameter and for the hyperparameters of the autoregressive partial autocorrelations (`Amu_trend` and `Aomega_trend`, which are described in detail by [Sarah Heaps in her wonderful VAR stationarity paper](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2079648)) as default. So our first model to consider is:


``` r
varmod <- mvgam(
  # Observation formula, empty to only consider the Gamma observation process
  formula = adj_count ~ 0,
  
  # Process model formula with regional intercepts and a VAR(1) process
  # that has a fully parameterized covariance matrix Sigma
  trend_formula = ~ region + VAR(cor = TRUE),
  
  # Modified prior distributions using brms::prior()
  prior = c(prior(std_normal(), class = b_trend),
            prior(exponential(2.5), class = sigma_trend)),
  
  # The time series data in 'long' format
  data = model_data,
  
  # A Gamma observation family with a log link. All three series
  # share the same Gamma shape parameter
  family = Gamma(link = 'log'),
  
  silent = 2,
  seed = 1)
```

### Model diagnostics and inferences
This model only takes a few seconds to fit four parallel Hamiltonian Monte Carlo chains. The sampler does warn that two chains had a low Bayesian fraction of missing information, and the summary shows why: the regional coefficients have low effective sample sizes. With autoregressive coefficients this close to one, the long-run level of each process is hard to separate from the latent states themselves. The estimates are still useful for this illustration, but longer chains would be worthwhile for a formal analysis:

``` r
summary(varmod)
```

```
##  Family: gamma 
##   Links: mu = log; shape = log 
## Formula: adj_count ~ 0 
##    Data: model_data (Number of observations: 120) 
##  Series: 3 
##  Trends: VAR(1); formula: ~region 
##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup draws = 4000
## 
## == Observation Model ==
## Further Distributional Parameters:
##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## shape    35.03      9.60    21.50    58.17 1.01      587     1111
## 
## == Trend Model ==
## Regression Coefficients:
##            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## regionAlb      0.69      0.12     0.41     0.89 1.04       98      144
## regionBC       0.82      0.18     0.43     1.14 1.04       88      166
## regionSask     0.31      0.08     0.12     0.46 1.04       96      221
## 
## Trend Specific Parameters:
##                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sigma_trend[Alb]      0.15      0.05     0.06     0.25 1.00      398      462
## sigma_trend[BC]       0.23      0.05     0.14     0.32 1.00      620     1388
## sigma_trend[Sask]     0.08      0.03     0.04     0.15 1.01      322      819
## Amu_trend[1,1]        0.99      0.37     0.29     1.76 1.01      514     1347
## Amu_trend[2,1]        0.61      0.29     0.08     1.23 1.00      339      694
## Aomega_trend[1,1]    20.17     16.52     1.59    62.49 1.00     1692     1662
## Aomega_trend[2,1]    24.79     17.30     3.63    70.00 1.00     1340     1631
## 
## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
## 
## Next steps:
##   - `pp_check(fit)`: posterior predictive checks
##   - `posterior_transition_matrix(fit)` / `irf(fit)`: lagged effects between series
##   - `residual_cor(fit)`: implied correlations among series
##   - `forecast(fit, newdata = ...)`: out-of-sample forecasts
##   - `lfo_cv(fit)`: leave-future-out model comparison
##   - `conditional_effects(fit)`: covariate effects
## Use `how_to_cite(fit)` for a methods paragraph and `methods_md(fit)` for its equations.
```

We can inspect posterior estimates for some of the key parameters using `bayesplot` functionality. For State-Space models, it is particularly recommended to inspect estimates for process and observation scale parameters (labelled as `sigma` and `shape` in the model's MCMC output):

``` r
mcmc_plot(varmod,
          variable = 'sigma_trend',
          regex = TRUE,
          type = 'trace')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-14-1.png" alt="mvgam posterior estimates using bayesplot" width="60%" style="display: block; margin: auto;" />


``` r
mcmc_plot(varmod,
          variable = 'shape',
          type = 'trace')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-1.png" alt="mvgam posterior estimates using bayesplot" width="60%" style="display: block; margin: auto;" />

Looking at possible coupling between these estimates is also a good idea, which you can do with the `pairs()` function:

``` r
pairs(varmod,
      variable = c('shape', 'sigma_trend'),
      regex = TRUE,
      off_diag_args = list(size = 1,
                           alpha = 0.5))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" alt="mvgam posterior pairs plots" width="60%" style="display: block; margin: auto;" />

Overall these estimates are quite stable and well-behaved. We can also inspect estimates for any linear predictor effects (in both the process and observation models) using `conditional_effects()`, which makes use of the wonderful [`marginaleffects`](https://marginaleffects.com/) package (see [this post on interpreting Generalized Additive Models using `marginaleffects`](https://ecogambler.netlify.app/blog/interpreting-gams/) for a glimpse into the full power of this package). Since we only have the varying intercepts as regressors in our model, we only get one plot from this call:

``` r
conditional_effects(varmod)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" alt="mvgam conditional_effects plots" width="60%" style="display: block; margin: auto;" />

By default `conditional_effects()` will return predictions on the response scale (i.e. considering the full uncertainty in posterior predictions). But we can also plot these on the link or expectation scales if we wish: 

``` r
conditional_effects(varmod,
                    type = 'link')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-18-1.png" alt="mvgam posterior pairs plots" width="60%" style="display: block; margin: auto;" />

Posterior predictive checks are also useful for interrogating a Bayesian model. `mvgam` offers multiple ways to do this. First we can compute unconditional posterior predictive checks, which integrate over the possible (stable) temporal states that *could have* been encountered for each series. In other words, these plots show us the types of predictions that we might make if we had no observations to inform the states of the latent VAR(1) process:

``` r
pp_check(varmod,
         type = "dens_overlay_grouped",
         group = "region",
         ndraws = 50,
         process_error = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-19-1.png" alt="mvgam posterior predictive checks" width="60%" style="display: block; margin: auto;" />

Obviously these prediction densities are more variable than the actual observation densities. But the main shapes seem to have been reasonably well captured. `mvgam` can also create conditional posterior predictive checks, using the actual latent temporal states that were estimated when the model conditioned on the observed data:

``` r
pp_check(varmod,
         type = "dens_overlay_grouped",
         group = "region",
         ndraws = 50,
         incl_autocor = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-20-1.png" alt="mvgam posterior predictive checks" width="60%" style="display: block; margin: auto;" />

These conditional checks also match the data well, which is expected when using a flexible autoregressive process to capture the unobserved dynamics. Next we can compute and plot probabilistic conditional hindcasts for each series:

``` r
hcs <- hindcast(varmod)
```


``` r
plot(hcs,
     series = 1)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-23-1.png" alt="mvgam probabilistic hindcasts" width="60%" style="display: block; margin: auto;" />


``` r
plot(hcs,
     series = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-25-1.png" alt="mvgam probabilistic hindcasts" width="60%" style="display: block; margin: auto;" />


``` r
plot(hcs,
     series = 3)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-27-1.png" alt="mvgam probabilistic hindcasts" width="60%" style="display: block; margin: auto;" />

### Model expansion
Given the subtle variations in the shapes of these observed distributions for each series, there may be support for fitting an equivalent model that allows the Gamma shape parameters to vary among series. Let's fit a second model that allows this:

``` r
varmod2 <- mvgam(
  # Observation formula, with a second formula that gives each
  # series its own Gamma shape parameter
  formula = bf(adj_count ~ 0,
               shape ~ 0 + series),
  
  # Process model formula as before
  trend_formula = ~ region + VAR(cor = TRUE),
  
  # Modified prior distributions, now including a prior for
  # the (log) shape parameters
  prior = c(prior(std_normal(), class = b_trend),
            prior(exponential(2.5), class = sigma_trend),
            prior(normal(2.5, 1), class = b, dpar = shape)),
  
  # The time series data in 'long' format
  data = model_data,
  
  # A Gamma observation family with a log link
  family = Gamma(link = 'log'),
  
  silent = 2,
  seed = 1)
```

This model encounters a few divergent transitions (about 1% of draws) but otherwise samples well. The shape parameters are estimated on the log scale:

``` r
summary(varmod2)
```

```
##  Family: gamma 
##   Links: mu = log; shape = log 
## Formula: adj_count ~ 0 
##          shape ~ 0 + series 
##    Data: model_data (Number of observations: 120) 
##  Series: 3 
##  Trends: VAR(1); formula: ~region 
##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup draws = 4000
## 
## == Observation Model ==
## Regression Coefficients:
##                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## shape_seriesAlb      3.22      0.47     2.51     4.39 1.01      441      686
## shape_seriesBC       3.39      0.55     2.42     4.64 1.01      652      559
## shape_seriesSask     3.58      0.33     3.00     4.30 1.00      779      720
## 
## == Trend Model ==
## Regression Coefficients:
##            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## regionAlb      0.69      0.13     0.40     0.91 1.05       91      273
## regionBC       0.81      0.21     0.37     1.20 1.05       92      196
## regionSask     0.30      0.10     0.07     0.47 1.05       95      217
## 
## Trend Specific Parameters:
##                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sigma_trend[Alb]      0.13      0.05     0.04     0.24 1.04      237      424
## sigma_trend[BC]       0.21      0.05     0.11     0.32 1.01      462      343
## sigma_trend[Sask]     0.08      0.03     0.04     0.16 1.02      306      521
## Amu_trend[1,1]        0.95      0.41     0.18     1.78 1.00      525      795
## Amu_trend[2,1]        0.66      0.31     0.10     1.34 1.02      256      262
## Aomega_trend[1,1]    18.39     16.18     1.15    60.36 1.00      960      645
## Aomega_trend[2,1]    23.92     17.91     2.70    70.60 1.01      690      634
## 
## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
## 
## Next steps:
##   - `pp_check(fit)`: posterior predictive checks
##   - `posterior_transition_matrix(fit)` / `irf(fit)`: lagged effects between series
##   - `residual_cor(fit)`: implied correlations among series
##   - `forecast(fit, newdata = ...)`: out-of-sample forecasts
##   - `lfo_cv(fit)`: leave-future-out model comparison
##   - `conditional_effects(fit)`: covariate effects
## Use `how_to_cite(fit)` for a methods paragraph and `methods_md(fit)` for its equations.
```

But now we have three shape parameters to inspect:

``` r
mcmc_plot(varmod2,
          variable = 'b_shape',
          regex = TRUE,
          type = 'trace')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-30-1.png" alt="mvgam posterior estimates using bayesplot" width="60%" style="display: block; margin: auto;" />

The unconditional posterior predictive checks look similar to those from the first model, with replicates for `'Sask'` matching the sharp peak in its observations a little more closely:

``` r
pp_check(varmod2,
         type = "dens_overlay_grouped",
         group = "region",
         ndraws = 50,
         process_error = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-31-1.png" alt="mvgam posterior predictive checks" width="60%" style="display: block; margin: auto;" />

And so do the probabilistic conditional hindcasts:

``` r
hcs <- hindcast(varmod2)
```


``` r
plot(hcs,
     series = 1)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-34-1.png" alt="mvgam probabilistic hindcasts" width="60%" style="display: block; margin: auto;" />


``` r
plot(hcs,
     series = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-36-1.png" alt="mvgam probabilistic hindcasts" width="60%" style="display: block; margin: auto;" />


``` r
plot(hcs,
     series = 3)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-38-1.png" alt="mvgam probabilistic hindcasts" width="60%" style="display: block; margin: auto;" />

### Model comparison with `loo_compare()`
The plots for the second model look *slightly* better for some of the series, but of course this is a more complex model so we may expect that. Generally we may want to select among these two competing models, and the approximate leave-one-out cross-validation routines provided by the [`loo`](https://mc-stan.org/loo/) package allow us to do this. We simply call `loo_compare()` and feed in the set of models that we'd like to compare:

``` r
loo_compare(varmod, varmod2)
```

```
##    model elpd_diff se_diff p_worse       diag_diff       diag_elpd
##  varmod2       0.0     0.0      NA                 22 k_psis > 0.7
##   varmod      -2.7     1.9    0.93 |elpd_diff| < 4 34 k_psis > 0.7
```

The second model has a higher expected log predictive density than the first (a difference of about 2.7, with a standard error of 1.9), so it appears to be slightly favoured. Notice the warnings about high Pareto `\(k\)` values, though. These are common for State-Space models, because each latent state is informed by very few observations, and they make the leave-one-out approximation less reliable. `lfo_cv()` offers a better comparison by [computing proper scoring rules for out-of-sample probabilistic forecasts](https://nicholasjclark.github.io/mvgam/articles/forecast_evaluation.html). But for simplicity we will stick with the second model for further interrogation.

## Interrogating VAR(1) models in `mvgam`
There are *many* useful questions that can be approached using VAR models, which I eluded to at the beginning of this post. `mvgam` has support for a growing number of investigations that can be performed using built-in functions, and this section will highlight some of those functionalities. 

### Inspecting AR coefficient and covariance estimates
One of the first things we might want to do with a VAR(1) model is to inspect the `\(A\)` matrix, which includes the autoregressive coefficients that estimate lagged dependence and cross-dependence between the elements of the latent process model. `posterior_transition_matrix()` returns posterior draws of this matrix, with rows naming the series that responds and columns naming the series that acts at the previous time step. A little reshaping puts the draws into a format for plotting the full matrix of estimates:

``` r
A_draws <- posterior_transition_matrix(varmod2, summary = FALSE)
A_df <- as.data.frame.table(A_draws, responseName = 'estimate') %>%
  dplyr::rename(draw = Var1, response = Var2, driver = Var3) %>%
  # Keep the central 99% of draws for each coefficient
  dplyr::group_by(response, driver) %>%
  dplyr::filter(estimate > quantile(estimate, 0.005),
                estimate < quantile(estimate, 0.995)) %>%
  dplyr::ungroup()
ggplot(A_df, aes(x = estimate)) +
  myhist(bins = 30) +
  facet_grid(response ~ driver,
             labeller = label_both,
             scales = 'free_x') +
  geom_vline(xintercept = 0,
             col = 'white',
             linewidth = 2) +
  geom_vline(xintercept = 0,
             linewidth = 1) +
  labs(x = 'Autoregressive coefficient', y = '') +
  hist_theme()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-40-1.png" alt="mvgam vector autoregressive coefficient inspection" width="60%" style="display: block; margin: auto;" />

Most of the posterior mass for these coefficients is positive, though the uncertainty is wide. The plot leaves out the most extreme 1% of draws for each coefficient, a few of which reach large values. The clearest lagged interactions run into the process for `'BC'`, from both `'Alb'` (posterior probability of a positive effect 0.89) and `'Sask'` (0.87), and into `'Alb'` from `'Sask'` (0.94). For example, the panel with `response: BC` and `driver: Sask` captures how an *increase* in the process for `'Sask'` at time `\(t\)` is expected to impact the process for `'BC'` at time `\(t+1\)`. This might make biological sense as it could be possible that increasing populations in one region lead to subsequent increases in other regions due to dispersal events.

We might also want to look at evidence for contemporaneous associations among these processes, which are captured in the covariance matrix `\(\Sigma_{process}\)`. The model stores the unique elements of this symmetric matrix, labelled by series:

``` r
mcmc_plot(varmod2,
          variable = c('Sigma_trend[Alb,Alb]',
                       'Sigma_trend[BC,Alb]',
                       'Sigma_trend[Sask,Alb]',
                       'Sigma_trend[BC,BC]',
                       'Sigma_trend[Sask,BC]',
                       'Sigma_trend[Sask,Sask]'),
          type = 'hist') +
  geom_vline(xintercept = 0,
             col = 'white',
             linewidth = 2) +
  geom_vline(xintercept = 0,
             linewidth = 1)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-41-1.png" alt="mvgam vector autoregressive covariance inspection" width="60%" style="display: block; margin: auto;" />

The contemporaneous covariances are small. The posterior probability of a positive covariance is 0.86 between `'BC'` and `'Alb'`, 0.79 between `'Sask'` and `'BC'` and only 0.48 between `'Sask'` and `'Alb'`. Most of the dependence among these processes therefore arises through the lagged interactions in `\(A\)` rather than through shared shocks within a year. The `residual_cor()` function summarises the correlations that the full model implies once both sources are combined.


### Impulse Response Functions
Generalized or Orthogonalized Impulse Response Functions (IRFs) can be computed using the posterior estimates of Vector Autoregressive parameters from `mvgam` models fitted with `VAR()` trend models. This function runs a simulation whereby a positive "shock" is generated for a target process at time `\(t = 0\)`. All else remaining stable, it then monitors how each of the remaining processes in the latent VAR would be expected to respond over the forecast horizon `\(h\)`. The function computes IRFs for all processes in the object and returns them in an array that can be plotted using the `S3` `plot()` function. Here we will use the generalized IRF, which makes no assumptions about the order in which the series appear in the VAR process, and inspect how each process is expected to respond to a sudden, positive *pulse* from the other processes over a horizon of 12 years.

``` r
irfs <- irf(varmod2,
            h = 12,
            orthogonal = FALSE)
```


``` r
plot(irfs,
     series = 1)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-44-1.png" alt="mvgam vector autoregressive impulse response functions" width="60%" style="display: block; margin: auto;" />


``` r
plot(irfs,
     series = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-46-1.png" alt="mvgam vector autoregressive impulse response functions" width="60%" style="display: block; margin: auto;" />


``` r
plot(irfs,
     series = 3)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-48-1.png" alt="mvgam vector autoregressive impulse response functions" width="60%" style="display: block; margin: auto;" />

These plots suggest that, if one region experienced a temporary pulse in the number of American kestrels, the remaining regions would tend to show increases over the following years, and that these effects would fade slowly. But the uncertainty is wide, and the 95% intervals include zero for most cross-region responses. Forty years of annual counts can only tell us so much about these interactions.

### Forecast Error Variance Decompositions
Using the same logic as above, we can inspect forecast error variance decompositions (FEVDs) for each process using the `fevd()` function. This type of analysis asks how orthogonal shocks to all process in the system contribute to the variance of forecast uncertainty for a focal process over increasing horizons. In other words, the proportion of the forecast variance of each latent time series can be attributed to the effects of the other series in the VAR proces. FEVDs are useful because some shocks may not be expected to cause variations in the short-term but may cause longer-term fluctuations

``` r
fevds <- fevd(varmod2, h = 12)
```

The `S3` `plot()` function for these objects returns a `ggplot` object that shows the mean expected contribution to forecast error variance for each process in the VAR model:


``` r
plot(fevds) +
  scale_fill_manual(values = c("#DCBCBC",
                               "#A25050",
                               "#5C0000"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-50-1.png" alt="mvgam vector autoregressive forecast error variance decompositions" width="60%" style="display: block; margin: auto;" />

This plot shows that the variance of forecast uncertainty for each process is initially dominated by contributions from that same process (i.e. self-dependent effects) but that effects from other processes become more important over increasing forecast horizons. Given what we saw from the IRF plots above, these long-term contributions from interactions among the processes makes sense.

### Community stability metrics
`mvgam` can also calculate a range of community stability metrics, which can be used to assess how important inter-series dependencies are to the variability of a multivariate system and to ask how systems are expected to respond to environmental perturbations. The `stability()` function computes a posterior distribution for the long-term stationary forecast distribution of the system, which
has mean `\(\mu_{\infty}\)` and variance `\(\Sigma_{\infty}\)`, to then calculate the following quantities:
- `prop_int`: Proportion of the volume of the stationary forecast distribution that is attributable to lagged interactions (i.e. how important are the autoregressive interaction coefficients in `\(A\)` for explaining the shape of the stationary forecast distribution?): `\(det(A)^2\)`
- `prop_int_adj`: Same as `prop_int` but scaled by the number of series `\(p\)` to facilitate direct comparisons among systems with different numbers of interacting variables: `\(det(A)^{2/p}\)`
- `prop_int_offdiag`: Sensitivity of `prop_int` to inter-series interactions (i.e. how important are the off-diagonals of the autoregressive coefficient matrix `\(A\)` for shaping `prop_int`?), calculated as the relative magnitude of the *off-diagonals* in the partial derivative matrix: `\([2~det(A) (A^{-1})^T]\)`
- `prop_int_diag`: Sensitivity of `prop_int` to intra-series interactions (i.e. how important are the
diagonals of the autoregressive coefficient matrix `\(A\)` for shaping `prop_int`?), calculated as the
relative magnitude of the *diagonals* in the partial derivative matrix: `\([2~det(A) (A^{-1})^T]\)`
- `prop_cov_offdiag`: Sensitivity of `\(\Sigma_{\infty}\)` to inter-series error correlations (i.e. how important are off-diagonal covariances in `\(\Sigma\)` for shaping `\(\Sigma_{\infty}\)`?), calculated as the relative magnitude of the *off-diagonals* in the partial derivative matrix: `\([2~det(\Sigma_{\infty}) (\Sigma_{\infty}^{-1})^T]\)`
- `prop_cov_diag`: Sensitivity of `\(\Sigma_{\infty}\)` to error variances (i.e. how important are diagonal variances in `\(\Sigma\)` for shaping `\(\Sigma_{\infty}\)`?), calculated as the relative magnitude of the *diagonals* in the partial derivative matrix: `\([2~det(\Sigma_{\infty}) (\Sigma_{\infty}^{-1})^T]\)`
- `reactivity`: A measure of the degree to which the system moves away from a stable equilibrium following a perturbation. Values `> 0` suggest the system is reactive, whereby a perturbation of the system in one period can be amplified into the next period. If `\(\sigma_{max}(A)\)` is the largest singular value of `\(A\)`, then reactivity is defined as: `\(log\sigma_{max}(A)\)`
- `mean_return_rate`: Asymptotic (long-term) return rate of the mean of the transition distribution to the stationary mean, calculated using the largest eigenvalue of the matrix `\(A\)`: `\(max(\lambda_{A})\)`. Lower values suggest greater stability of the system
- `var_return_rate`: Asymptotic (long-term) return rate of the variance of the transition distribution to the stationary variance: `\(max(\lambda_{A \otimes{A}})\)`. Again, lower values suggest greater stability of the system

``` r
metrics <- stability(varmod2, summary = FALSE)
```

We can plot the posterior distribution of reactivity, which is [perhaps the most commonly calculated metric of ecological community stability](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/13-0996.1), using `ggplot()`:

``` r
ggplot(metrics, aes(x = reactivity)) +
  myhist() +
  labs(x = 'Reactivity', y = '') +
  geom_vline(xintercept = 0, col = 'white', linewidth = 2.5) +
  geom_vline(xintercept = 0, linewidth = 1.5) +
  hist_theme()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-52-1.png" alt="mvgam vector autoregressive community reactivity estimates" width="60%" style="display: block; margin: auto;" />

This system is considered to be "reactive", whereby this index of the community's short-term response to perturbations is mostly estimated to be `>0`. Again this is not surprising, as both the `irf()` and `fevd()` plots show that the time series of adjusted counts in each region is highly dependent on time-lagged adjusted counts in other regions.

Another useful question to ask from these metrics is what are the relative importances of the diagonals for determining the contributions of the `\(A\)` and `\(\Sigma\)` matrices to the shape of the stationary forecast distribution. In other words, we can decompose the contributions of these two matrices to the system's stochastic equilibrium to understand which components (i.e. the self-dependence vs the inter-dependence effects) seem to have the biggest influences on the variability of the system. What follows is some rather long-winded code that will calculate and plot these contributions, using some helpful example code provided by [Matthew Kay](https://www.mjskay.com/), maintainer of the [`tidybayes`](https://mjskay.github.io/tidybayes/) package

``` r
# Gather posterior estimates of relative contributions for each matrix
# into a tidy format for tidybayes plotting
prop_df <- rbind(metrics %>%
    dplyr::select(prop_cov_diag,
                  prop_cov_offdiag) %>%
    dplyr::rename(diagonals = prop_cov_diag,
                  off_diagonals = prop_cov_offdiag) %>%
    dplyr::mutate(.draw = 1:dplyr::n(),
                  group = 'Process error (Sigma)') %>%
    tidyr::pivot_longer(c(diagonals, off_diagonals),
                        names_to = 'Parameters',
                        values_to = 'proportion'),
  metrics %>%
    dplyr::select(prop_int_diag,
                  prop_int_offdiag) %>%
    dplyr::rename(diagonals = prop_int_diag,
                  off_diagonals = prop_int_offdiag) %>%
    dplyr::mutate(.draw = 1:dplyr::n(),
                  group = 'Autoregressive coefs (A)') %>%
    tidyr::pivot_longer(c(diagonals, off_diagonals),
                        names_to = 'Parameters',
                        values_to = 'proportion'))

# Calculate blurred color maps to graphically represent uncertainty
# in these contributions
p <- seq(0, 1 - .Machine$double.eps, length.out = 200)
df_curves <- prop_df %>%
  dplyr::mutate(cumulative_proportion = cumsum(proportion),
                .by = c(.draw, group)) |>
  dplyr::reframe(
    `P(value)` = (1 - ecdf(cumulative_proportion)(p)) *
      ecdf(cumulative_proportion - proportion)(p),
    proportion = p,
    .by = c(group, Parameters))

# Set the colours and their labels
fill_scale <- scale_fill_manual(values = c("lightgrey", "#A25050"))
fill_scale$train(c("diagonals", "off_diagonals"))
```

Now we can plot these relative contributions in a way that aims to represent their associated posterior uncertainties

``` r
df_curves %>%
  summarise(color = encode_colour(`P(value)` %*%
                            decode_colour(fill_scale$map(Parameters),
                                          to = "oklab"),
                          from = "oklab"),
    .by = c(group,
            proportion)) %>%
  ggplot() +
  geom_blank(aes(fill = Parameters),
             data = data.frame(Parameters = c("diagonals",
                                              "off_diagonals"))) +
  geom_slab(aes(y = proportion,
                x = group,
                fill = I(color),
                group = NA),
    thickness = 1,
    fill_type = "auto",
    side = "both") +
  stat_pointinterval(
    aes(x = group,
        y = proportion,
        group = Parameters),
    data = prop_df |>
      mutate(proportion = cumsum(proportion),
             .by = c(.draw, group)) |>
      filter(Parameters != "off_diagonals")) +
  fill_scale +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(expand = FALSE) +
  labs(x = '',
       y = 'Contribution to stationary distribution')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-54-1.png" alt="mvgam vector autoregressive community reactivity estimates" width="60%" style="display: block; margin: auto;" />

And there we have it. This plot demonstrates that, within the contribution of the `\(A\)` matrix to the stable distribution, the off-diagonals (i.e. the lagged autoregressive interaction effects) play a sizable role. However, within the contribution of the `\(\Sigma\)` matrix to the stable distribution, the diagonals (i.e. the process variances) are much more important than are the off-diagonal covariances. This gives us some insights into the relative importances of lagged vs contemporaneous "interactions" and can hopefully facilitate more in-depth ecological analyses. 

Also, although we haven't covered it in this example, there is of course no reason why you cannot include complex covariate effects in either the process or observation models. Indeed, `mvgam()` can readily handle a wide range of effects including multidimensional penalized smooth functions, Gaussian Process functions, monotonic smooth functions, [distributed lag nonlinear effects](https://ecogambler.netlify.app/blog/distributed-lags-mgcv/), seasonality and [time-varying seasonality](https://ecogambler.netlify.app/blog/time-varying-seasonality/), and hierarchical effects. It is also easy to generate probabilistic forecasts from these models, which you can learn more about by reading the [forecast and forecast evaluation vignette](https://nicholasjclark.github.io/mvgam/articles/forecast_evaluation.html). These options should hopefully make `mvgam` a very attractive resource for tackling real-world multivariate time series when we do not wish to ignore important drivers or use time series transformations that make both inference and prediction difficult. Thanks for your time!

## Further reading
The following papers and resources offer useful material about Vector Autoregressive models and Dynamic GAMs for time series modeling / forecasting
  
Hampton, S. E., Holmes, E. E., Scheef, L. P., Scheuerell, M. D., Katz, S. L., Pendleton, D. E., & Ward, E. J. (2013). [Quantifying effects of abiotic and biotic drivers on community dynamics with multivariate autoregressive (MAR) models](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/13-0996.1). *Ecology*, 94(12), 2663-2669.  
  
Hannaford, N. E., Heaps, S. E., Nye, T. M., Curtis, T. P., Allen, B., Golightly, A., & Wilkinson, D. J. (2023). [A sparse Bayesian hierarchical vector autoregressive model for microbial dynamics in a wastewater treatment plant](https://www.sciencedirect.com/science/article/pii/S0167947322002390). *Computational Statistics & Data Analysis*, 179, 107659.
  
Heaps, S. E. (2023). [Enforcing stationarity through the prior in vector autoregressions](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2079648). *Journal of Computational and Graphical Statistics*, 32(1), 74-83.  
  
Karunarathna, K. A. N. K., Wells, K., & Clark, N. J. (2024). [Modelling nonlinear responses of a desert rodent species to environmental change with hierarchical dynamic generalized additive models](https://www.sciencedirect.com/science/article/pii/S0304380024000371). *Ecological Modelling*, 490, 110648.  
  
Lütkepohl, H (2006). [New introduction to multiple time series analysis](https://link.springer.com/book/10.1007/978-3-540-27752-1). Springer, New York.  
  
Pesaran, P. H., & Yongcheol, S. (1998). [Generalized impulse response analysis in linear multivariate models](https://www.sciencedirect.com/science/article/pii/S0165176597002140) *Economics Letters* 58: 17-29. 
