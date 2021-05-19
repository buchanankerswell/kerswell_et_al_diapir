# Load functions and libraries
source('functions.R')

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data', pattern = '_marx.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Summarise marker features
purrr::map2(models, paths, ~{
  # Load markers
  load_marx(.y)
  # Take the markers dataframe and ...
  get(paste0(.x, '.marx')) %>%
  # Compute features
  marx_ft() -> m
  # Remove markers and grids
  rm(list = paste0(.x, '.marx'), envir = .GlobalEnv)
  rm(list = paste0(.x, '.grid'), envir = .GlobalEnv)
  return(m)
}) %>%
purrr::set_names(models) -> marx.features

# Save
cat('\nSaving marker features to data/marx_features.RData')
save(marx.features, file = 'data/marx_features.RData')

cat('\nDone!')