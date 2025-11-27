---
title: "Time series modeling with Bayesian Dynamic Generalized Additive Models"
excerpt: "In this talk I introduce Bayesian Dynamic Generalized Additive Models (DGAMs) and illustrate their advantages for analyzing and forecasting real-world time series. I discuss mvgam, an open-source R package that can fit DGAMs with nonlinear effects, hierarchical effects and dynamic processes to data from a wide variety of observation distributions. These models are especially useful for analysing multiple series, as they can estimate hierarchical smooth functions while learning complex temporal associations with latent vector autoregressive processes or dimension-reduced dynamic factor processes. Because the package uses Hamiltonian Monte Carlo inference through Stan, it is straightforward to create Stan code and all necessary data structures so that additional stochastic elements can be added to suit the user's bespoke needs. Other key features of {mvgam} are functions to critique models using rolling window forecasts and posterior predictive checks, online data augmentation via a recursive particle filter and graphical tools to visualise probabilistic uncertainties for smooth functions and predictions. I hope show how models that describe real-world complexity, both through nonlinear covariate functions and multi-series dependence, are useful to ask targeted questions about drivers of change."
subtitle: "A talk on timeseries analysis and forecasting with the mvgam R 📦 for the 2023 Australian Statistical Conference"
date: 2023-12-14T14:15:59-06:00
date_end: "2023-12-14T14:45:59-06:00"
show_post_time: false
author: "Nicholas Clark"
location: "Wollongong, Australia"
draft: false
# layout options: single, single-sidebar
layout: single
categories:
- talks
- mvgam
- time-series
links:
- icon: door-open
  icon_pack: fas
  name: website
  url: https://nicholasjclark.github.io/ASC_talk/ASC_talk_slidedeck#1
- icon: github
  icon_pack: fab
  name: code
  url: https://github.com/nicholasjclark/ASC_talk/tree/main?tab=readme-ov-file
---

