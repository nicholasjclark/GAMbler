---
title: "MRFcov"
subtitle: "R 📦 for network analysis using Conditional Random Fields (CRFs)"
excerpt: "R 📦 for network analysis using Conditional Random Fields (CRFs)"
date: 2021-06-21
author: "Nicholas Clark"
draft: false
tags:
  - MRFcov
  - package
categories:
  - R
  - package
  - MRFcov
layout: single
links:
- icon: door-open
  icon_pack: fas
  name: website
  url: https://cran.r-project.org/web/packages/MRFcov/index.html
- icon: github
  icon_pack: fab
  name: code
  url: https://github.com/nicholasjclark/MRFcov
---

This goal of `MRFcov` is to approximate interaction parameters of nodes in undirected Markov Random Fields (MRF) graphical networks. Models can incorporate covariates (a class of models known as Conditional Random Fields; CRFs; following methods developed by [Cheng et al 2014](https://academic.oup.com/biometrics/article/70/4/943/7419933) and [Lindberg 2016](https://www.utupub.fi/handle/10024/124199)), allowing users to estimate how interactions between nodes are predicted to change across covariate gradients. .
  
This package was primarily designed for ecological networks to infer interactions between co-occurring species, which is key to identifying processes governing community assembly. Markov random fields (MRFs) are especially useful for estimating interspecific partial correlations and how these change across environmental gradients. For more details, see the primary publication that describes their implementation ([Clark et al 2016](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.2221)) 
