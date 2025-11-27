# GAMbler Blog Source

This repository contains the source files for the [GAMbler blog](https://ecogambler.netlify.app/), a data science and statistical modeling blog by Nicholas Clark.

## About the Blog

GAMbler is a research-focused blog emphasizing ecological modeling and statistical methods, with particular attention to:

- **Generalized Additive Models (GAMs)** - Comprehensive tutorials and applications
- **Ecological time series analysis** - Methods for analyzing temporal ecological data
- **Bayesian inference** - Applications in ecological modeling
- **R package development** - Particularly the `mvgam` package for multivariate GAMs
- **Statistical visualization** - Creating effective data visualizations in R

## Technical Stack

- **Framework**: Hugo with the Apéro theme
- **Content**: R Markdown posts with reproducible code examples
- **Deployment**: Netlify
- **Languages**: R, Stan

## Author

**Nicholas Clark**  
University of Queensland

## Building Locally

The site uses blogdown for local development. Key files:
- Blog posts: `content/blog/YYYY-MM-DD-slug/index.Rmarkdown`
- Configuration: `config.yaml`
- R build scripts: `R/build.R`, `R/build2.R`

## License

See individual blog posts for specific licensing information.