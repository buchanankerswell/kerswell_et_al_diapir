source('functions.R')

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('/Volumes/hd/nmods/kerswell_et_al_marx/data', pattern = '_marx.RData', full.names = T)[1:4]
models <- stringr::str_extract(paths, 'cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Detect number of cores
cores <- parallel::detectCores()
cat('\nParallel computation with', cores, 'cores ...')

# Summarise marker features
fun <- function(model, path, n, thresh, k) {
  # Load markers
  load_marx(path)
  # Crop tsteps after interference timecut
  marx.df <- get(paste0(model, '.marx'))
  tcut <- attr(marx.df, 'tcut')
  # Compute features
  fts.df <- marx.df %>% marx_ft(features = c('sum.dP', 'max.P'))
  # Classify markers
  marx.class.df <- marx_classify(marx.df, fts.df, thresh = thresh, k = k)
  # Monte Carlo sampling marx_classify
  mc <- monte_carlo(marx.df, fts.df, n = n, thresh = thresh, k = k)
  # Save
  assign(paste0(model, '.marx.classified'), list(
    'mc' = mc,
    'marx' = marx.class.df$marx,
    'gm' = marx.class.df$mcl
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
  k = 10,
  n = 100,
  thresh = 1,
  mc.cores = cores,
  SIMPLIFY = F
)

cat('\nDone!')
