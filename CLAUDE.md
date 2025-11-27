# GAMbler Blog - Claude Code Assistant Instructions

## Project Overview
Personal blog focused on statistical modeling, data analysis and visualization, particularly:
- Generalized Additive Models (GAMs) 
- Time series analysis
- Ecological modeling
- Bayesian inference
- R package development (especially mvgam)

## Primary Commands
```bash
# Knit RMarkdown files
Rscript -e "rmarkdown::render('filename.Rmd')"

# Run R scripts
Rscript script.R
```

## R Coding Style
- **Style guide**: Tidyverse style guide
- **Line width**: Maximum 80 characters
- **Function calls**: Multi-argument functions with each argument on a new line
- **Namespace**: Explicit namespace for dplyr functions (e.g., `dplyr::filter()`, NOT `filter()`)
- **Pipes**: Base R pipe `|>` preferred
- **Package loading**: Use `library()` to load packages
- **Variable names**: snake_case
- **Comments**: Liberal use to explain statistical concepts
- **Spacing**: Spaces around operators (`<-`, `+`, `=`)
- **Object naming**: Clear names that reflect content
- **Code progression**: Build from simple to complex

## ggplot2 Styling
```r
# Standard theme
theme_set(
  theme_classic(
    base_size = 15, 
    base_family = 'serif'
  ) +
  theme(
    axis.line.x.bottom = element_line(
      colour = "black",
      size = 1
    ),
    axis.line.y.left = element_line(
      colour = "black",
      size = 1
    )
  )
)

# Color palettes: viridis (options: "D", "plasma", etc.)
```

## Blog Post Structure

### Directory Structure
- **Format**: `YYYY-MM-DD-slug/index.Rmarkdown`
- **Workflow**: index.Rmarkdown → knitted to markdown → blogdown renders site

### YAML Header Structure
```yaml
---
title: "Descriptive title using title case"
author: "Nicholas Clark"
subtitle: ""
excerpt: "Detailed 2-3 sentence description of post content, methods, and key value proposition"
slug: lowercase-hyphenated-slug
date: YYYY-MM-DD
draft: false
images:
series:
tags:
  - rstats
  - [package name: mgcv, mvgam, brms, etc.]
  - tutorial/simulation as appropriate
  - topic tags (forecasting, time-series, ecology, etc.)
categories:
  - rstats
  - [main package or topic]
layout: single-sidebar
---
```

### Content Structure
1. **Headings**: Hierarchical markdown (H1-H3), descriptive and action-oriented
2. **Sections**:
   - Overview/Purpose (introduction)
   - Environment setup (libraries needed)
   - Main content with step-by-step progression
   - Visualization of results
   - Interpretation
   - "Further reading" section (ALWAYS include)
3. **Code blocks**: Well-commented, syntax-highlighted R code
4. **Mathematical notation**: LaTeX for equations

### R Code Chunk Options

# Setup chunk (always included)
```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,   
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 6,
  out.width = "60%",
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)
```

# Common chunk options:
- `cache = TRUE` for computationally intensive chunks
- `include = FALSE` for setup or behind-the-scenes computation
- `eval = FALSE` to show code without running it
- `echo = FALSE` to hide code but show output
- `fig.alt = "Descriptive text"` for accessibility (ALWAYS include)

### Figure Alt Text Pattern
- **Always include** `fig.alt` for accessibility
- **Format**: Describe what the visualization shows + tools/packages used
- **Examples**:
  - `fig.alt = "Visualising distributed lag smooths in mgcv and R"`
  - `fig.alt = "GAM predictions showing nonlinear feature adoption effects by contract tier"`
  - `fig.alt = "Comparison of prediction intervals between Gaussian and Tweedie GAMs for CLV"`

## Writing Style
- **Tone**: Conversational yet technical
- **Perspective**: First-person narrative with occasional humor
- **Structure**: Step-by-step progression from theory to practice
- **Explanations**: Start with theoretical background, implement in code, visualize results, interpret findings
- **Pedagogical approach**: Break down complex concepts, use practical examples
- **Internal links**: Frequently link to other relevant blog posts (e.g., "If you've worked with [GAMs for time series analysis](https://ecogambler.netlify.app/blog/autocorrelated-gams/) like I have...")
- **External links**: Include many hyperlinks to papers, documentation, and resources

## Hyperlink Requirements
- **CRITICAL**: ALL links MUST be verified as correct BEFORE being added to posts
- **Internal links**: Use format `https://ecogambler.netlify.app/blog/[post-slug]/`
- **External links**: Verify URLs work and point to intended resource
- **Academic papers**: Link to DOIs when available
- **Link text**: Use descriptive text that explains what the link points to

## Key Technologies
- **Languages**: R, Stan
- **Packages**: mvgam, mgcv, brms, tidyverse, ggplot2
- **Focus areas**: GAMs, time series, Bayesian modeling, compositional analysis

## Tasks to Assist With
- Generating simulated data for demonstrations
- Writing statistical modeling code
- Creating data visualizations
- Formatting mathematical equations
- Proofreading and editing
- Suggesting relevant literature
- Debugging R code
- Explaining statistical concepts
- Verifying all hyperlinks before inclusion