source('functions.R')

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('/Volumes/hd/nmods/kerswell_et_al_marx/data', pattern = '_marx.RData', full.names = T)
models <- stringr::str_extract(paths, 'cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Detect number of cores
cores <- parallel::detectCores()
cat('\nParallel computation with', cores, 'cores ...')

# Summarise marker features
fun <- function(model, path, n, p, thresh, k) {
  # Load markers
  load_marx(path)
  # Crop tsteps after interference timecut
  marx <- get(paste0(model, '.marx'))
  tcut <- attr(marx, 'tcut')
  # Compute features
  fts <- marx %>% marx_ft(features = c('sum.dP', 'max.P'))
  # Classify markers
  marx.class <- marx_classify(marx, fts, tcut, thresh, k)
  # Monte Carlo sampling marx_classify
  jk <- jknife(marx, thresh, n, p, k)
  # Save
  assign(paste0(model, '.marx.classified'), list(
    'jk' = jk,
    'marx' = marx.class$marx,
    'mcl' = marx.class$mcl,
    'misclass' = marx.class$misclass,
    'reclass' = marx.class$reclass,
    'lowerB' = marx.class$lowerB
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
  n = 10,
  p = 0.98,
  thresh = 1,
  mc.cores = cores,
  SIMPLIFY = F
)

cat('\nDone!')
