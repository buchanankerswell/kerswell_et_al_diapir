# Load functions and libraries
source('functions.R')

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data', pattern = '_marx.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Detect number of cores
cores <- parallel::detectCores()
cat('\nParallel computation with', cores, 'cores ...')

# Summarise marker features
fun <- function(model, path) {
  # Load markers
  load_marx(path)
  # Take the markers dataframe and ...
  get(paste0(model, '.marx')) %>%
  # Compute features
  marx_ft(features = c(
    'tsteps',
    'above.fourty.kbar',
    'above.seven.hundred.c',
    'down.dx'
    #'rundown.dx',
    #'rundown.dT'
    )
  ) -> m
  # Remove markers and grids
  rm(list = paste0(model, '.marx'), envir = .GlobalEnv)
  rm(list = paste0(model, '.grid'), envir = .GlobalEnv)
  return(m)
}

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  mc.cores = cores,
  SIMPLIFY = F
) %>%
purrr::set_names(models) -> marx.features

# Print features
print(marx.features)

# Save
cat('\nSaving marker features to data/marx_features.RData')
save(marx.features, file = 'data/marx_features.RData')

cat('\nDone!')