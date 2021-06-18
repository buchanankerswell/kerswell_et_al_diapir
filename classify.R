source('functions.R')

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data', pattern = '_marx.RData', full.names = T)[1:32]
models <- stringr::str_extract(paths, 'cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Detect number of cores
cores <- parallel::detectCores()
cat('\nParallel computation with', cores, 'cores ...')

# Summarise marker features
fun <- function(model, path, n, k) {
  # Load markers
  load_marx(path)
  # Take the markers dataframe and ...
  marx.df <- get(paste0(model, '.marx'))
  # Compute features
  fts.df <- marx.df %>%
  marx_ft(features = c('sum.dP', 'max.P'))
  # Classify markers
  marx.class.df <- marx_classify(marx.df, fts.df, k = k)
  # Monte Carlo sampling marx_classify
  mc <- monte_carlo(marx.df, fts.df, n = n, k = k)
  # Save
  assign(paste0(model, '.marx.classified'), list(
    'mc' = mc,
    'marx' = marx.class.df$marx,
    'gm' = marx.class.df$mc
  ))
  # Save
  cat('\nSaving classified markers to data/', model, '_k', k,  '_marx_classified.RData', sep = '')
  do.call(
    save,
    list(paste0(model, '.marx.classified'),
         file = paste0('data/', model, '_k', k, '_marx_classified.RData')
    )
  )
}

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  k = 6,
  n = 2,
  mc.cores = cores,
  SIMPLIFY = F
) %>% invisible()

cat('\nDone!')
