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
  marx.df <- get(paste0(model, '.marx'))
  # Compute features
  marx.df %>%
  marx_ft(features = c(
    #'above.thirty.kbar',
    #'above.four.hundred.c',
    #'runup.dP'
    #'max.T',
    'sum.dP',
    'max.P'
    #'down.dx'
    #'rundown.dx'
    #'rundown.dT'
    )
  ) -> fts.df
  # Classify markers
  marx.class.df <- marx_classify(marx.df, fts.df)
  # Monte Carlo sampling marx_classify
  mc <- monte_carlo(marx.df, fts.df, n = 1000)
  # Save
  return(list(
    'mc' = mc,
    'marx' = marx.class.df$marx,
    'gm' = marx.class.df$mc
  ))
}

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  mc.cores = cores,
  SIMPLIFY = F
) %>%
purrr::set_names(models) -> marx.classified

# Save
cat('\nSaving classifiec markers to data/marx_classified.RData')
save(marx.classified, file = 'data/marx_classified.RData')

cat('\nDone!')