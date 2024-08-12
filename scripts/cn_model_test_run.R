# Load Libraries ---------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(visdat)
library(here)
library(lubridate)
library(readr)
# library(naniar)
library(purrr)

# if (!require(remotes)) {
#   install.packages("remotes")
# }
detach("package:rsofun", unload = TRUE)
remotes::install_github("stineb/rsofun", ref = "cnmodel")
library(rsofun)


# Load Data------------

# FLUXNET
drivers <- readRDS(here("data", "drivers_ch_oe2.rds"))

# Prepare driver object --------------
## N deposition----------


# Reactive N input needs to be in gN per day
# added to forcing time series: specify quantity of N added on which day

# Example N input given a constant rate each day
n_input_test <- function(data) {
  data <- data %>%
    mutate(forcing = purrr::map(
      forcing,
      ~mutate(
        .,
        fharv = 0.0,
        dno3 = 0.002003263,
        dnh4 = 0.002017981))
    )
  return(data)
}

drivers <- n_input_test(drivers)


## Harvesting------------


# The fraction of biomass harvested per day needs to be specified in the forcing time series
# cseed and nseed new seeds added after harvesting
# Example driver update assumes harvesting is 0 and new seeds planted after harvest

fharv_seed <- function(data, use_cseed = 5, cn_seed = 20) {
  use_nseed <- use_cseed / cn_seed

  data <- data %>%
    mutate(
      forcing = purrr::map(
        forcing, ~mutate(
          .,
          fharv = ifelse(month(date) == 7 & mday(date) == 15, 0.0, 0.0),
          cseed = ifelse(month(date) == 2 & mday(date) == 15, use_cseed, 0.0),
          nseed = ifelse(month(date) == 2 & mday(date) == 15, use_nseed, 0.0))
      )
    )

  return(data)
}

drivers <- fharv_seed(drivers)


## Simulation parameters------------

# The spinup of cn_model must be long enough to equilibrate fluxes

# Function to modify specific columns in each dataframe
modify_params <- function(df_list, spinupyears_val, recycle_val) {
  # Map over each dataframe in the list
  df_list <- map(df_list, ~ {
    # Modify specific columns
    mutate(.x,
           spinupyears = spinupyears_val,
           recycle = recycle_val)
  })

  return(df_list)
}

# FLUXNET
drivers$params_siml <- modify_params(drivers$params_siml, 2021, 2)


## Model parameters -------------

pars <- list(
  # P-model
  kphio = 0.07,             # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b = 1.0,
  soilm_thetastar = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim = 30.0,
  kc_jmax = 0.41,

  # Plant
  f_nretain = 0.500000,
  fpc_tree_max = 0.950000,
  growtheff = 0.600000,
  r_root = 2*0.913000,
  r_sapw = 2*0.044000,
  exurate = 0.003000,

  k_decay_leaf = 1.90000,
  k_decay_root = 1.90000,
  k_decay_labl = 1.90000,
  k_decay_sapw = 1.90000,

  r_cton_root = 37.0000,
  r_cton_wood = 100.000,
  r_cton_seed = 15.0000,
  nv_vcmax25 = 0.02 * 13681.77, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.695; previously: 5000.0,
  ncw_min = 0.08 * 1.116222, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.691; previously used: 0.056,
  r_n_cw_v = 0, # assumed that LMA is independent of Vcmax25; previously: 0.1,
  r_ctostructn_leaf = 1.3 * 45.84125, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.699; previously used: 80.0000,
  kbeer = 0.500000,

  # Phenology (should be PFT-specific)
  gddbase = 5.0,
  ramp = 0.0,
  phentype = 2.0,

  # Soil physics (should be derived from params_soil, fsand, fclay, forg, fgravel)
  perc_k1 = 5.0,
  thdiff_wp = 0.2,
  thdiff_whc15 = 0.8,
  thdiff_fc = 0.4,
  forg = 0.01,
  wbwp = 0.029,
  por = 0.421,
  fsand = 0.82,
  fclay = 0.06,
  fsilt = 0.12,

  # Water and energy balance
  kA = 107,
  kalb_sw = 0.17,
  kalb_vis = 0.03,
  kb = 0.20,
  kc = 0.25,
  kCw = 1.05,
  kd = 0.50,
  ke = 0.0167,
  keps = 23.44,
  kWm = 220.0,
  kw = 0.26,
  komega = 283.0,
  maxmeltrate = 3.0,

  # Soil BGC
  klitt_af10 = 1.2,
  klitt_as10 = 0.35,
  klitt_bg10 = 0.35,
  kexu10 = 50.0,
  ksoil_fs10 = 0.021,
  ksoil_sl10 = 7.0e-04,
  ntoc_crit1 = 0.45,
  ntoc_crit2 = 0.76,
  cton_microb = 10.0,
  cton_soil = 9.77,
  fastfrac = 0.985,

  # N uptake
  eff_nup = 0.0001000,
  minimumcostfix = 1.000000,
  fixoptimum = 25.15000,
  a_param_fix = -3.62000,
  b_param_fix = 0.270000,

  # Inorganic N transformations (re-interpreted for simple ntransform model)
  maxnitr =  0.00005,

  # Inorganic N transformations for full ntransform model (not used in simple model)
  non = 0.01,
  n2on = 0.0005,
  kn = 83.0,
  kdoc = 17.0,
  docmax = 1.0,
  dnitr2n2o = 0.01,

  # Additional parameters - previously forgotten
  frac_leaf = 0.5,         # after wood allocation
  frac_wood = 0,           # highest priority in allocation
  frac_avl_labl = 0.1,

  # for development
  tmppar = 9999,

  # simple N uptake module parameters
  nuptake_kc = 600,
  nuptake_kv = 5,
  nuptake_vmax = 0.2
)

# Run the model
# # for these parameters using ch-oe1
#
# # Create output directories
#
# # Define the full path
# full_path <- "/Users/PhillipZywczuk/Documents/2023_2024_postdoc_eth/data/renku/n2o_ssa_lit_review/lit_review_modelling/p_model/data"
#
# # Create the full path and all necessary parent directories
# dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
#
# # Create 'out' directory inside the specified full path
# out_dir <- file.path(full_path, "out")
# dir.create(out_dir, showWarnings = FALSE)
#
# # Create 'vignettes' directory inside the specified full path
# vignettes_dir <- file.path(full_path, "vignettes")
# dir.create(vignettes_dir, showWarnings = FALSE)
#
# # Create 'out' directory inside the 'vignettes' directory
# vignettes_out_dir <- file.path(vignettes_dir, "out")
# dir.create(vignettes_out_dir, showWarnings = FALSE)

# C-only run------------------------
# Define whether to use interactive C-N cycling
drivers$params_siml[[1]]$c_only <- TRUE

# Run the model
output <- runread_cnmodel_f(drivers, par = pars)


## Aggregate outputs -----------
# Extract data and aggregate by mean
df_out <- output$data[[1]] %>%
  as_tibble()

# overall mean
# XXX: may aggregate to growing-season only
df_mean <- df_out |>
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

# mean seasonal cycle
doydf_rsofun <- df_out |>
  as_tibble() |>
  mutate(doy = lubridate::yday(date)) |>
  group_by(doy) |>
  summarise(gpp = mean(gpp), fapar = mean(fapar))

# Benchmark --------------------
## Leaf traits -------------------

# obtain file from https://doi.org/10.5281/zenodo.6831903
df_traits <- read_csv("~/data/leafn_vcmax_ning_dong/data_leafn_vcmax_ning_dong.csv") %>%
  rename(lat = Latitude, lon = longitude)

# Compare to data for non-woody species given that CH-Oe2 is a grassland
df_traits <- df_traits |>
  filter(woody == "non-woody") |>
  rename(lma = LMA, narea = Narea, nmass = Nmass, vcmax = vcmax_obs)

# LMA (rsofun output is in gC m-2-leaf, Ning's data is in "mass" = DM?)
gg_traits_1 <- ggplot() +
  geom_density(aes(lma/2, ..density..), data = df_traits, fill = "#777055ff", color = NA, alpha = 0.5) +
  theme_classic() +
  geom_vline(aes(xintercept = lma), data = df_mean, color = "#29a274ff") +
  labs(x = "LMA") +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(expand = c(0, 0))

gg_traits_2 <- ggplot() +
  geom_density(aes(nmass, ..density..), data = df_traits, fill = "#777055ff", color = NA, alpha = 0.5) +
  theme_classic() +
  geom_vline(aes(xintercept = narea/lma), data = df_mean, color = "#29a274ff") +
  labs(x = expression(italic("N")[mass])) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(expand = c(0, 0))

gg_traits_3 <- ggplot() +
  geom_density(aes(narea, ..density..), data = df_traits, fill = "#777055ff", color = NA, alpha = 0.5) +
  geom_vline(aes(xintercept = narea), data = df_mean, color = "#29a274ff") +
  theme_classic() +
  labs(x = expression(italic("N")[area])) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(expand = c(0, 0))

## Vcmax in micro-mol m-2 s-1
gg_traits_4 <- ggplot() +
  geom_density(aes(vcmax, ..density..), data = df_traits, fill = "#777055ff", color = NA, alpha = 0.5) +
  theme_classic() +
  geom_vline(aes(xintercept =  1e6 * vcmax), data = df_mean, color = "#29a274ff") +
  labs(x = expression(italic("V")[cmax])) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(expand = c(0, 0))


# XXX simulated LMA is high in comparison to the mean across all observations in the Dong et al. dataset.
# May change suitable parameter. Or is it consistent with observations from the site?

## fAPAR mean seasonal cycle -----------
doydf_modis <- drivers |>
  unnest(forcing) |>
  mutate(doy = lubridate::yday(date)) |>
  group_by(doy) |>
  summarise(fapar = mean(fapar, na.rm = TRUE))

gg_fapar_msc <- ggplot() +
  geom_line(data = doydf_modis,  aes(doy, fapar), color = "#777055ff") +
  geom_line(data = doydf_rsofun, aes(doy, fapar), color = "#29a274ff") +
  theme_classic()

# XXX: too high fAPAR (and therefore LAI), adjust parameters

## GPP mean seasonal cycle ---------------
doydf_fluxnet <- drivers |>
  unnest(forcing) |>
  mutate(doy = lubridate::yday(date)) |>
  group_by(doy) |>
  summarise(gpp = mean(gpp))

gg_gpp_msc <- ggplot() +
  geom_line(data = doydf_fluxnet, aes(doy, gpp), color = "#777055ff") +
  geom_line(data = doydf_rsofun, aes(doy, gpp), color = "#29a274ff") +
  theme_classic()

# XXX mean seasonal cycle of observed GPP looks like there is harvesting in early summer

## NEE mean seasonal cycle ---------------

# xxx to be included in rsofun drivers

## Combine plots -----------
panel_traits <- plot_grid(
  gg_traits_1,
  gg_traits_2,
  gg_traits_3,
  gg_traits_4,
  ncol = 1
)

plot_grid(
  panel_traits,
  gg_fapar_msc,
  gg_gpp_msc,
  ncol = 1,
  rel_heights = c(1,0.7,0.7)
)

ggsave(
  here("output/benchmarking_c-only.pdf"),
  width = 5,
  height = 12
  )

# Further workflow ---------------
# - Adjust relevant parameters to achieve good C-only model benchmarking
# - Adjust relevant parameters to reduce N fixation to minimum with un-closed N balance (given that we don't have much N fixation at this site.)
# - Adjust relevant parameters to achieve good C-only model benchmarking with closed N balance


