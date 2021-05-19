# Load functions and libraries
source('functions.R')

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data', pattern = '_marx.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Summarise marker features
purrr::walk2(models, paths, ~{
  # Load markers
  load_marx(.y)
  marx <- get(paste0(.x, '.marx'))
  # Marker motion movie
#   marx %>% filter(tstep %in% seq(20, 30, 1)) %>% marx_motion_mov(.x)
  # Marker histogram movie
#   marx %>% marx_boxplot_mov(.x)
  # Marker PT movie
#   marx %>% marx_PT_mov(.x)
  # Marker features plot
  ft <- marx %>% marx_ft()
  ft %>% marx_features_plot(.x)
  # Correlation matrix
  cat('\nPlotting correlation matrix [', .x, ']', sep = '')
  cm <- GGally::ggpairs(ft)
  ggsave(
    plot = cm,
    filename = paste0('figs/features/', .x, '_correlations.png'),
    device = 'png',
    type = 'cairo',
    width = 48,
    height = 48,
    units = 'in',
    dpi = 300
  )
})

cat('\nDone!')