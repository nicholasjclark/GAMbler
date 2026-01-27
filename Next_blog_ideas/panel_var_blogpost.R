# We'll use WDI data from the World Bank, dplyr and mvgam
library(WDI)
library(tidyverse)
library(mvgam)

# Grab country-level GDP, consumption and investment time series
gdp_cons_inv <- WDI(
  indicator = c(
    "NY.GDP.MKTP.KN",
    "NE.CON.TOTL.KN",
    "NE.GDI.FTOT.KN"
  ),
  start = 1970
)

# Take differenced logs of the target indicators
model_dat <- gdp_cons_inv %>%
  dplyr::mutate(
    dl_gdp = c(NA, diff(log(NY.GDP.MKTP.KN))),
    dl_cons = c(NA, diff(log(NE.CON.TOTL.KN))),
    dl_gfcf = c(NA, diff(log(NE.GDI.FTOT.KN)))
  ) %>%
  dplyr::arrange(country, year) %>%
  tidyr::pivot_longer(
    cols = tidyselect::contains("dl_"),
    names_to = "indicator",
    values_to = "value"
  ) %>%
  dplyr::select(
    country,
    indicator,
    value,
    year
  ) %>%
  dplyr::filter(
    country %in%
      c(
        "United States",
        "United Kingdom",
        "Australia",
        "New Zealand",
        "Chile",
        "Canada",
        "Ireland",
        "South Africa"
      ),
    year > 1970
  ) %>%
  dplyr::mutate(
    time = year,
    country = as.factor(country),
    indicator = factor(
      indicator,
      levels = c(
        "dl_gdp",
        "dl_cons",
        "dl_gfcf"
      )
    )
  )

# Split into training (using same data as Savage, up to 2016)
# and testing (three years ahead, avoiding unpredictable COVID drop)
data_train <- model_dat %>%
  dplyr::filter(year <= 2016)
data_test <- model_dat %>%
  dplyr::filter(
    year > 2016 & year < 2020
  )

# Set up a zero-centred hierarchical panel VAR model
# (a useful extension would be to let the process sigmas vary by type of
# indicator in a hierarchical model, which mvgam currently cannot do)
mod <- mvgam(
  formula = value ~ 0,
  trend_model = VAR(
    gr = country,
    subgr = indicator
  ),
  # Assume observation error is very small
  # and process error is primary source of
  # unmodelled variation
  priors = c(
    prior(
      exponential(100),
      class = sigma_obs
    ),
    prior(
      exponential(5),
      class = sigma
    )
  ),
  family = gaussian(),
  share_obs_params = TRUE,
  data = data_train,
  newdata = data_test,
  control = list(adapt_delta = 0.7)
)

# View the posterior distribution for the correlation pooling parameter
mcmc_plot(
  mod,
  type = "dens",
  variable = "alpha_cor"
) +
  xlim(c(0, 1))

# View some country-level correlation matrices
all_sigmas <- as.matrix(
  mod,
  variable = "Sigma",
  regex = TRUE
)

all_cors <- lapply(
  1:NROW(all_sigmas),
  function(x) {
    cov2cor(
      matrix(
        all_sigmas[x, ],
        nrow = nlevels(mod$obs_data$series),
        ncol = nlevels(mod$obs_data$series)
      )
    )
  }
)

mean_cor <- Reduce("+", all_cors) /
  length(all_cors)

# Australian correlation matrix
round(
  mean_cor[1:3, 1:3],
  3
)

# Chile correlation matrix
round(
  mean_cor[7:9, 7:9],
  3
)

# New Zealand correlation matrix
round(
  mean_cor[13:15, 13:15],
  3
)

# USA correlation matrix
round(
  mean_cor[22:24, 22:24],
  3
)

# View some country-level correlation matrices
all_vars <- as.matrix(
  mod,
  variable = "A",
  regex = TRUE
)
all_vars <- lapply(
  1:NROW(all_vars),
  function(x) {
    matrix(
      all_vars[x, ],
      nrow = nlevels(mod$obs_data$series),
      ncol = nlevels(mod$obs_data$series)
    )
  }
)

mean_var <- Reduce("+", all_vars) /
  length(all_cors)

# Australian VAR matrix
round(
  mean_var[1:3, 1:3],
  3
)

# Chile VAR matrix
round(
  mean_var[7:9, 7:9],
  3
)

# New Zealand VAR matrix
round(
  mean_var[13:15, 13:15],
  3
)

# USA VAR matrix
round(
  mean_var[22:24, 22:24],
  3
)

# FEVDs
fevds <- fevd(mod, h = 6)
fevd_quants <- summary(fevds)
fevd_quants %>%
  dplyr::bind_cols(
    data.frame(
  do.call(
    rbind,
    strsplit(fevd_quants$shock, ' -> ')
  )
) %>%
  dplyr::rename(series = X1,
                target = X2)
) -> fevd_quants

# Function to convert latent process names to indicators
convert_labels = function(input, which_series){
  out <- gsub(
    paste0('Process_',
           which_series[1]),
    'GDP',
    input
  )
  out <- gsub(
    paste0('Process_',
           which_series[2]),
    'CONS',
    out
  )
  out <- gsub(
    paste0('Process_',
           which_series[3]),
    'GFCF',
    out
  )
  out
}

# Function to filter the FEVD summary and plot median contributions to
# forecast variance
plot_fevds = function(which_series = 1:3,
                      title = 'Australia'){
  ggplot2::ggplot(
    fevd_quants %>%
      dplyr::filter(
        series %in% paste0('Process_', which_series)
      ) %>%
      dplyr::mutate(
        target = convert_labels(target, which_series),
        Indicator = convert_labels(series, which_series)
      ) %>%
      dplyr::filter(fevdQ50 != 0),
    ggplot2::aes(
      fill = Indicator, y = fevdQ50,
      x = horizon
    )
  ) +
    ggplot2::geom_bar(
      position = "stack",
      stat = "identity"
    ) +
    ggplot2::facet_wrap(~target) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = 'Forecast horizon',
      y = 'Median contribution to forecast variance',
      title = title
    ) +
    ggplot2::scale_fill_manual(
      values = c("#DCBCBC",
                 "#A25050",
                 "#5C0000")
    )
}

# Australian FEVDs
plot_fevds()

# Chile FEVDs
plot_fevds(
  which_series = 7:9,
  title = 'Chile'
)

# New Zealand FEVDs
plot_fevds(
  which_series = 13:15,
  title = 'New Zealand'
)

# USA FEVDs
plot_fevds(
  which_series = 22:24,
  title = 'USA'
)

# IRFs
irfs <- irf(mod, h = 6)
irf_quants <- summary(irfs)
irf_quants %>%
  dplyr::bind_cols(
    data.frame(
      do.call(
        rbind,
        strsplit(irf_quants$shock, ' -> ')
      )
    ) %>%
      dplyr::rename(series = X1,
                    target = X2)
  ) -> irf_quants

plot_irfs = function(which_series = 1:3,
                     title = 'Australia'){
  ggplot2::ggplot(
    data = irf_quants %>%
      dplyr::filter(
        series %in% paste0('Process_', which_series)
      ) %>%
      dplyr::mutate(
        target = paste0('Response:',
                        convert_labels(target, which_series)),
        Indicator = paste0('Shock:',
                           convert_labels(series, which_series))
      ) %>%
      dplyr::filter(irfQ50 != 0),
    ggplot2::aes(x = horizon, y = irfQ50)
  ) +
    ggplot2::geom_ribbon(
      mapping = ggplot2::aes(
        ymin = irfQ2.5,
        ymax = irfQ97.5
      ),
      fill = "#C79999"
    ) +
    ggplot2::geom_line(
      col = "#8F2727",
      linewidth = 1
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::facet_grid(
      target ~ Indicator,
      scales = 'free_y',
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::labs(
      x = "Horizon",
      y = "Generalized impulse response",
      title = title
    ) +
    ggplot2::theme_bw()
}

# Australian IRFs
plot_irfs()

# Chile IRFs
plot_irfs(
  which_series = 7:9,
  title = 'Chile'
)

# New Zealand IRFs
plot_irfs(
  which_series = 13:15,
  title = 'New Zealand'
)

# USA IRFs
plot_irfs(
  which_series = 22:24,
  title = 'USA'
)

# Look at some hindcasts / forecasts
fcs <- forecast(mod)

# Australian forecasts
patchwork::wrap_plots(
  plot(fcs, series = 1),
  plot(fcs, series = 2),
  plot(fcs, series = 3),
  ncol = 2
)

# Chile forecasts
patchwork::wrap_plots(
  plot(fcs, series = 7),
  plot(fcs, series = 8),
  plot(fcs, series = 9),
  ncol = 2
)

# New Zealand forecasts
patchwork::wrap_plots(
  plot(fcs, series = 13),
  plot(fcs, series = 14),
  plot(fcs, series = 15),
  ncol = 2
)

# USA forecasts
patchwork::wrap_plots(
  plot(fcs, series = 22),
  plot(fcs, series = 23),
  plot(fcs, series = 24),
  ncol = 2
)
