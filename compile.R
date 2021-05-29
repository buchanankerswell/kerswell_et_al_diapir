
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
  mc <- monte_carlo(marx.df, fts.df, n = 2)
  # Save
  assign(paste0(model, '.marx.classified'), list(
    'mc' = mc,
    'marx' = marx.class.df$marx,
    'gm' = marx.class.df$mc
  ))
  # Save
  cat('\nSaving classified markers to data/', model, '_marx_classified.RData', sep = '')
  do.call(
    save,
    list(paste0(model, '.marx.classified'),
         file = paste0('data/', model, '_marx_classified.RData')
    )
  )
}

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  mc.cores = cores,
  SIMPLIFY = F
) %>% invisible()

cat('\nDone!')