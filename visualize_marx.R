# Load functions and libraries
source('functions.R')

# Read Penniston-Dorland et al., 2015 dataset
pd15 <- readr::read_delim('data/PD15.tsv', delim = '\t', col_types = 'cddcccd')
pd15 <- pd15[1:nrow(pd15)-1,]

tibble(T = sort(pd15$temperature)) %>%
mutate(cdf = (row_number()-1)/n()) -> pd15T

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data/k6_classd', pattern = '.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

#cat('\nFound classified markers for:', models, sep = '\n')

paths <- paths[61:64]
models <- models[61:64]

# Detect number of cores
cores <- parallel::detectCores()
cat('\nParallel computation with', cores, 'cores ...')

fun <- function(model, path) {
  # Load markers and grids
  load_marx(paste0('/Volumes/hd/nmods/kerswell_et_al_marx/data/', model, '_marx.RData'))
  load(path)
  # Marker data
  marx <- get(paste0(model, '.marx.classified'))$marx
  # Classification info
  m <- get(paste0(model, '.marx.classified'))$gm
  # Timecut
  tcut <- attr(get(paste0(model, '.marx')), 'tcut')
  # Nodes data
  grid <- get(paste0(model, '.grid'))[[tcut]] %>% mutate(z = z - 18000)
  # Max P summary for CDFs
  get(paste0(model, '.marx.classified'))$mc %>%
  purrr::map_df(~.x$cdfP, .id = 'run') -> maxP
  # Max T summary for CDFs
  get(paste0(model, '.marx.classified'))$mc %>%
  purrr::map_df(~.x$cdfT, .id = 'run') -> maxT
  # Visualize GMM results
  png(
    file = paste0('figs/', model, '_class.png'),
    type = 'cairo',
    res = 300,
    width = 5,
    height = 5,
    units = 'in'
  )
  plot(m, what = 'classification')
  dev.off()
  # Sum of pressure changes recovered
  p1 <- marx %>%
  group_by(id, recovered) %>%
  summarise(sumdP = sum(diff(P)), .groups = 'keep') %>%
  ggplot(aes(x = sumdP/1e4, fill = recovered)) +
  geom_histogram(bins = 100) +
  labs(
    title = paste0('Sum of all dP [', model, ']'),
    x = bquote(sum(dP, '', '')~'[GPa]'),
    y = 'Frequency',
    fill = 'Recovered'
  ) +
  scale_fill_grey(start = 0, end = 0.6) +
  theme_classic()
  ggsave(
    paste0('figs/', model, '_sumdP.png'),
    plot = p1,
    device = 'png',
    type = 'cairo',
    width = 5,
    height = 5,
    dpi = 300
  )
  # Sum of pressure changes classes
  p2 <- marx %>%
  group_by(id, class) %>%
  summarise(sumdP = sum(diff(P)), .groups = 'keep') %>%
  ggplot(aes(x = sumdP/1e4, fill = as.factor(class))) +
  scale_fill_manual(values = wesanderson::wes_palette('IsleofDogs1')) +
  geom_histogram(bins = 100) +
  labs(
    title = paste0('Sum of all dP [', model, ']'),
    x = bquote(sum(dP, '', '')~'[GPa]'),
    y = 'Frequency',
    fill = 'Class'
  ) +
  theme_classic()
  ggsave(
    paste0('figs/', model, '_sumdP_class.png'),
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 5,
    height = 5,
    dpi = 300
  )
  # Max pressure recovered
  p3 <- marx %>%
  group_by(id, recovered) %>%
  summarise(maxP = max(P), .groups = 'keep') %>%
  ggplot(aes(x = maxP/1e4, fill = recovered)) +
  geom_histogram(bins = 100) +
  labs(
    title = paste0('Max P [', model, ']'),
    x = 'Max P [GPa]',
    y = 'Frequency',
    fill = 'Recovered'
  ) +
  scale_fill_grey(start = 0, end = 0.6) +
  theme_classic()
  ggsave(
    paste0('figs/', model, '_maxP.png'),
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 5,
    height = 5,
    dpi = 300
  )
  # Max pressure class
  p4 <- marx %>%
  group_by(id, class) %>%
  summarise(maxP = max(P), .groups = 'keep') %>%
  ggplot(aes(x = maxP/1e4, fill = as.factor(class))) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = wesanderson::wes_palette('IsleofDogs1')) +
  labs(
    title = paste0('Max P [', model, ']'),
    x = 'Max P [GPa]',
    y = 'Frequency',
    fill = 'Class'
  ) +
  theme_classic()
  ggsave(
    paste0('figs/', model, '_maxP_class.png'),
    plot = p4,
    device = 'png',
    type = 'cairo',
    width = 5,
    height = 5,
    dpi = 300
  )
  # Composite 1D classification
  p <- (
    p1 +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()) +
    p3 +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank())
  ) /
  (
    p2 + theme(plot.title = element_blank()) +
    p4 +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank())
  ) +
  plot_layout(guides = 'collect') +
  plot_annotation(title = paste0('Marker classification [', model, ']'), tag_levels = 'a')
  ggsave(
    paste0('figs/', model, '_class_comp.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 5,
    dpi = 300
  )
  # Marker maxP CDF
  p <- maxP %>%
  group_by(run) %>%
  ggplot(aes(x = maxP/1e4, y = cdf)) +
  geom_path(aes(linetype = 'Markers', group = run), show.legend = F) +
  geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
  scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'Markers' = 'solid')) +
  labs(
    title = paste0('Cumulative probability of max P [', model, ']'),
    x = 'Maximum P [GPa]',
    y = 'Probability'
  ) +
  theme_classic()
  ggsave(
    paste0('figs/', model, '_cdfP.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 5,
    height = 5,
    dpi = 300
  )
  # Marker maxT CDF
  p <- maxT %>%
  group_by(run) %>%
  ggplot(aes(x = maxT - 273, y = cdf)) +
  geom_path(aes(linetype = 'Markers', group = run), show.legend = F) +
  geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
  labs(
    title = paste0('Cumulative probability of max T [', model, ']'),
    x = 'Maximum T [C]',
    y = 'Probability'
  ) +
  scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'Markers' = 'solid')) +
  theme_classic()
  ggsave(
    paste0('figs/', model, '_cdfT.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 5,
    height = 5,
    dpi = 300
  )
  # Viscosity
  p <- grid %>%
  draw_grid(
    model = model,
    marx = marx,
    class = 'recovered',
    box = c(up = -18, down = 200, left = 500, right = 1800),
    time = tcut,
    bk.alpha = 0.9,
    mk.alpha = 0.5,
    mk.size = 0.2,
    leg.dir = 'horizontal',
    leg.dir.rec = 'horizontal',
    sub.col = 'deeppink',
    rec.col = 'white',
    leg.pos = 'bottom',
    p.type = 'viscosity',
    v.pal = 'viridis',
    transparent = F)
  ggsave(
    paste0('figs/', model, '_viscosity_profile.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 3,
    dpi = 300
  )
  # Draw cross section at timecut
  p1 <- grid %>%
  draw_grid(
    model = model,
    marx = marx,
    class = 'recovered',
    time = tcut,
    box = c(up = -18, down = 200, left = 500, right = 1800),
    bk.alpha = 0.9,
    mk.size = 0.6,
    sub.col = 'deeppink',
    rec.col = 'white',
    leg.pos = 'bottom',
    leg.dir = 'horizontal',
    p.type = 'stream',
    v.pal = 'viridis',
    transparent = F)
  p2 <- grid %>%
  draw_grid(
    model = model,
    class = 'recovered',
    time = tcut,
    box = c(up = -18, down = 200, left = 500, right = 1800),
    mk.alpha = 1,
    bk.alpha = 0.9,
    leg.pos = 'bottom',
    leg.dir = 'horizontal',
    p.type = 'temperature',
    v.pal = 'magma',
    transparent = F)
  p3 <- grid %>%
  draw_grid(
    model = model,
    class = 'recovered',
    time = tcut,
    box = c(up = -18, down = 200, left = 500, right = 1800),
    bk.alpha = 0.9,
    leg.pos = 'bottom',
    leg.dir = 'horizontal',
    p.type = 'viscosity',
    v.pal = 'viridis',
    transparent = F)
  (p1 + theme(axis.title.x = element_blank())) /
  (p2  + theme(axis.title.x = element_blank())) /
  p3 + plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') -> p
  ggsave(
    paste0('figs/', model, '_comp_profile.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7.5,
    height = 7,
    dpi = 300
  )
  # Marker motion movie
  #marx %>% marx_motion_mov(model, class = T)
  # Marker boxplot movie
  #marx %>% marx_boxplot_mov(model)
  # Marker PT movie
  #marx %>% marx_PT_mov(model, class = T)

  # Clean up environment
  rm(list =
    c(
    paste0(model, '.marx.classified'),
    paste0(model, '.grid'),
    paste0(model, '.marx')
    )
  )
}

#purrr::map2(models, paths, fun)

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  mc.cores = cores,
  SIMPLIFY = F
)

cat('\nDone!')