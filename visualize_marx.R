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

fun <- function(model, path) {
  # Load markers
  load_marx(path)
  marx <- get(paste0(model, '.marx'))
  # Marker motion movie
#   marx %>% marx_motion_mov(model)
  # Marker histogram movie
#   marx %>% marx_boxplot_mov(model)
  # Marker PT movie
#   marx %>% marx_PT_mov(model)
  # Marker features plot
  ft <- marx %>% marx_ft()
#   ft %>% marx_features_plot(model)
  # Correlation matrix
  cat('\nPlotting correlation matrix [', model, ']', sep = '')
  cm <- GGally::ggpairs(ft)
  ggsave(
    plot = cm,
    filename = paste0('figs/features/', model, '_correlations.png'),
    device = 'png',
    type = 'cairo',
    width = 48,
    height = 48,
    units = 'in',
    dpi = 72
  )
}

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  mc.cores = cores,
  SIMPLIFY = F
) %>%
invisible()

cat('\nDone!')