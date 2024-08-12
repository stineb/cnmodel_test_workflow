library(tidyverse)
library(here)

# FluxDataKit drivers v3.3 Obtain file from Zenodo https://doi.org/10.5281/zenodo.12818273
drivers <- read_rds("~/Downloads/rsofun_driver_data_v3.3.rds")

# subset
drivers_ch_oe2 <- drivers |>
  filter(sitename == "CH-Oe2")

# safe file contained in this repo
write_rds(
  drivers_ch_oe2,
  file = here("data/drivers_ch_oe2.rds")
)


