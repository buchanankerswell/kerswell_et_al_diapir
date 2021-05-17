# Load functions and libraries
source('functions.R')

# Load marker and grid data
paths <- list.files('data', pattern = '_marx.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

# Summarise marker features
purrr::map2(models, paths, ~{
  # Load markers
  load_marx(.y)
  # Compute features
  get(paste0(.x, '.marx')) %>%
  marx_ft()
  # Remove markers and grids
  rm(list = paste0(.x, '.marx'), envir = .GlobalEnv)
  rm(list = paste0(.x, '.grid'), envir = .GlobalEnv)
}) %>%
purrr::set_names(models) -> marx.features

# Save
cat('\nSaving marker features to data/marx_features.RData')
save(marx.features, file = 'data/marx_features.RData')