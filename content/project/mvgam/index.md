---
title: "mvgam"
subtitle: "R 📦 to fit Dynamic Bayesian Generalised Additive Models for time series analysis and forecasting"
excerpt: "R 📦 to fit Dynamic Bayesian Generalised Additive Models for time series analysis and forecasting"
date: 2024-03-12
author: "Nicholas Clark"
draft: false
tags:
  - mvgam
  - package
categories:
  - R
  - package
  - mgcv
  - Stan
  - JAGS
layout: single
links:
- icon: door-open
  icon_pack: fas
  name: website
  url: https://nicholasjclark.github.io/mvgam/
- icon: github
  icon_pack: fab
  name: code
  url: https://github.com/nicholasjclark/mvgam/
---

The goal of mvgam is to use a Bayesian framework to estimate parameters of Dynamic Generalized Additive Models (DGAMs) for time series with dynamic trend components. The package provides an interface to fit Bayesian DGAMs using either [`JAGS`](https://mcmc-jags.sourceforge.io/) or [`Stan`](https://mc-stan.org/) as the backend, but note that users are strongly encouraged to opt for Stan over JAGS. The formula syntax is based on that of the package [`mgcv`](https://cran.r-project.org/web/packages/mgcv/index.html) to provide a familiar GAM modelling interface. There is also built-in support for the increasingly powerful [`marginaleffects` package](https://marginaleffects.com/) to make interpretation easy. The motivation for the package and some of its primary objectives are described in detail by Clark & Wells 2022 (published in [Methods in Ecology and Evolution](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13974)). An introduction to the package and some worked examples are also shown in the below seminar: 

**Ecological Forecasting with Dynamic Generalized Additive Models DGAMs)**

{{% youtube "0zZopLlomsQ" %}}
