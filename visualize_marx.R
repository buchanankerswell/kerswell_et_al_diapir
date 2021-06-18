# Load functions and libraries
source('functions.R')

# Read Penniston-Dorland et al., 2015 dataset
pd15 <- readr::read_delim('data/PD15.tsv', delim = '\t', col_types = 'cddcccd')
pd15 <- pd15[1:nrow(pd15)-1,]

tibble(T = sort(pd15$temperature)) %>%
mutate(cdf = (row_number()-1)/n()) -> pd15T

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data/k6', pattern = '.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

#cat('\nFound classified markers for:', models, sep = '\n')

paths <- paths[61:64]
models <- models[61:64]

# path <- paths[5]
# model <- models[5]

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
#   png(
#     file = paste0('figs/', model, '_class.png'),
#     type = 'cairo',
#     res = 300,
#     width = 5,
#     height = 5,
#     units = 'in'
#   )
#   plot(m, what = 'classification')
#   dev.off()
  cat('\nDrawing classification results [', model, ']', sep = '')
  # Sum of pressure changes recovered
  tibble(
    class = seq_len(length(m$parameters$mean[1,])),
    sumdP = m$parameters$mean[1,],
    maxP = m$parameters$mean[2,]
  ) -> cent
  d1 <- marx %>%
  group_by(id, recovered) %>%
  summarise(sumdP = sum(diff(P)), .groups = 'keep')
  d2 <- marx %>%
  group_by(id, class) %>%
  summarise(sumdP = sum(diff(P)), .groups = 'keep')
  d3 <- marx %>%
  group_by(id, recovered) %>%
  summarise(maxP = max(P), .groups = 'keep')
  d4 <- marx %>%
  group_by(id, class) %>%
  summarise(maxP = max(P), .groups = 'keep')
  p1 <-  d1 %>%
  ggplot(aes(x = sumdP/1e4)) +
  geom_histogram(aes(fill = recovered), bins = 100) +
  labs(
    x = bquote('sumdP'~'[GPa]'),
    y = 'Frequency',
    fill = 'Recovered'
  ) +
  scale_fill_grey(start = 0.6, end = 0) +
  theme_classic(base_size = 8)
  # Sum of pressure changes classes
  p2 <- d2 %>%
  ggplot(aes(x = sumdP/1e4)) +
  geom_histogram(aes(fill = as.factor(class)), bins = 100) +
  geom_rug(data = cent, aes(x = sumdP/1e4, color = as.factor(class)), show.legend = F) +
  geom_segment(
    aes(
      linetype = 'centroid threshold',
      x = (median(d2$sumdP) - 3/4*IQR(d2$sumdP))/1e4,
      xend = (median(d2$sumdP) - 3/4*IQR(d2$sumdP))/1e4,
      y = -Inf,
      yend = Inf),
    alpha = 0.5) +
  labs(
    x = bquote('sumdP'~'[GPa]'),
    y = 'Frequency',
    fill = 'Class',
    linetype = NULL
  ) +
  scale_linetype_manual(values = c('reclassify threshold' = 'dotted', 'centroid threshold' = 'solid')) +
  scale_color_manual(values = wesanderson::wes_palette(10, name = 'IsleofDogs1', type = 'continuous')) +
  scale_fill_manual(values = wesanderson::wes_palette(10, name = 'IsleofDogs1', type = 'continuous')) +
  theme_classic(base_size = 8)
  # Max pressure recovered
  p3 <- d3 %>%
  ggplot(aes(x = maxP/1e4)) +
  geom_histogram(aes(fill = recovered), bins = 100) +
  geom_segment(
    aes(
      linetype = 'reclassify threshold',
      x = (median(d3$maxP))/1e4,
      xend = (median(d3$maxP))/1e4,
      y = -Inf,
      yend = Inf),
    alpha = 0.5) +
  labs(
    x = bquote('maxP [GPa]'),
    y = 'Frequency',
    fill = 'Recovered',
    linetype = NULL
  ) +
  scale_linetype_manual(values = c('reclassify threshold' = 'dotted', 'centroid threshold' = 'solid')) +
  scale_fill_grey(start = 0.6, end = 0) +
  theme_classic(base_size = 8)
  # Max pressure class
  p4 <- d4 %>%
  ggplot(aes(x = maxP/1e4)) +
  geom_histogram(aes(fill = as.factor(class)), bins = 100) +
  geom_rug(data = cent, aes(x = maxP/1e4, color = as.factor(class)), show.legend = F) +
  geom_segment(
    aes(
      linetype = 'centroid threshold',
      x = (median(d4$maxP) - 3/4*IQR(d4$maxP))/1e4,
      xend = (median(d4$maxP) - 3/4*IQR(d4$maxP))/1e4,
      y = -Inf,
      yend = Inf),
    alpha = 0.5) +
  scale_linetype_manual(values = c('reclassify threshold' = 'dotted', 'centroid threshold' = 'solid')) +
  scale_color_manual(values = wesanderson::wes_palette(10, name = 'IsleofDogs1', type = 'continuous')) +
  scale_fill_manual(values = wesanderson::wes_palette(10, name = 'IsleofDogs1', type = 'continuous')) +
  labs(
    x = bquote('maxP [GPa]'),
    y = 'Frequency',
    fill = 'Class',
    linetype = NULL
  ) +
  theme_classic(base_size = 8)
  # Composite 1D classification
  p <- (
    p4 +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()) +
    p2 +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank())
  ) /
  (
    p3 +
    p1 +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank())
  ) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a')
  cat('\nSaving classification plot [', model, ']', sep = '')
  ggsave(
    paste0('figs/', model, '_class_comp.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 5,
    dpi = 300
  )
  # Marker maxP CDF
  p1 <- maxP %>%
  group_by(run) %>%
  ggplot(aes(x = maxP/1e4, y = cdf)) +
  geom_path(aes(linetype = 'Markers', group = run), show.legend = F) +
  geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
  scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'Markers' = 'solid')) +
  labs(
    x = 'Maximum P [GPa]',
    y = 'Probability'
  ) +
  theme_classic(base_size = 8)
# Marker maxT CDF
  p2 <- maxT %>%
  group_by(run) %>%
  ggplot(aes(x = maxT - 273, y = cdf)) +
  geom_path(aes(linetype = 'Markers', group = run), show.legend = F) +
  geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
  labs(
    x = 'Maximum T [C]',
    y = 'Probability'
  ) +
  scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'Markers' = 'solid')) +
  theme_classic(base_size = 8)
# Viscosity
  p3 <- grid %>%
  draw_grid(
    model = model,
    marx = marx,
    class = 'recovered',
    box = c(up = -18, down = 200, left = 500, right = 1800),
    time = tcut,
    bk.alpha = 0.7,
    mk.alpha = 1,
    mk.size = 0.25,
    leg.dir = 'horizontal',
    leg.dir.rec = 'horizontal',
    sub.col = 'deeppink',
    rec.col = 'black',
    leg.pos = 'bottom',
    p.type = 'viscosity',
    v.pal = 'viridis',
    base.size = 8,
    transparent = F)
  p <- p3 / (p1 + (p2 +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
  ))) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(heights = c(1, 1.8), guides = 'collect') &
  theme(legend.position = 'bottom')
  ggsave(
    paste0('figs/', model, '_comp.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 5,
    dpi = 300
  )
# Draw cross section at timecut
  p1 <- grid %>%
  draw_grid(
    model = model,
    marx = marx,
    mk.size = 0.25,
    rec.col = 'black',
    sub.col = 'deeppink',
    class = 'recovered',
    time = tcut,
    box = c(up = -18, down = 200, left = 500, right = 1800),
    bk.alpha = 1,
    bk.col = rgb(0, 0, 0, 0.1),
    leg.pos = 'bottom',
    leg.dir = 'horizontal',
    leg.dir.rec = 'horizontal',
    p.type = 'stream',
    iso.alpha = 0.6,
    iso.col = 'black',
    stm.alpha = 0.2,
    v.pal = 'viridis',
    base.size = 8,
    transparent = F)
  p3 <- grid %>%
  draw_grid(
    model = model,
    class = 'recovered',
    time = tcut,
    box = c(up = -18, down = 200, left = 500, right = 1800),
    bk.alpha = 1,
    leg.pos = 'bottom',
    leg.dir = 'horizontal',
    p.type = 'viscosity',
    iso.alpha = 0.6,
    v.pal = 'viridis',
    base.size = 8,
    transparent = F)
  (p1 + theme(axis.title.x = element_blank())) /
  p3 + plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') -> p
  ggsave(
    paste0('figs/', model, '_comp_profile.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 4,
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