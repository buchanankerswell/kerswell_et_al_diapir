# Load functions and libraries
source('functions.R')

# Read Penniston-Dorland et al., 2015 dataset
pd15 <- readr::read_delim('data/PD15.tsv', delim = '\t', col_types = 'cddcccd')
pd15 <- pd15[1:nrow(pd15)-1,]

tibble(T = sort(pd15$temperature)) %>%
mutate(cdf = (row_number()-1)/n()) -> pd15T

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data/k10', pattern = '.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

cat('\nFound classified markers for:', models, sep = '\n')

# Detect number of cores
cores <- parallel::detectCores()
cat('\nParallel computation with', cores, 'cores ...')

fun <- function(model, path) {
  # Load markers and grids
  load_marx(paste0('/Volumes/hd/nmods/kerswell_et_al_marx/data/', model, '_marx.RData'))
  load(path)
  # Classification threshold factor median * thresh
  thresh <- 1
  # Marker data
  marx <- get(paste0(model, '.marx.classified'))$marx
  # Classification info
  mcl <- get(paste0(model, '.marx.classified'))$mcl
  # Timecut
  tcut <- attr(get(paste0(model, '.marx')), 'tcut')
  t2 <- ceiling(tcut/2)
  t3 <- ceiling(tcut/4)
  # Nodes data
  grid <- get(paste0(model, '.grid'))[[tcut]] %>% mutate(z = z - 18000)
  g2 <- get(paste0(model, '.grid'))[[t2]] %>% mutate(z = z - 18000)
  g3 <- get(paste0(model, '.grid'))[[t3]] %>% mutate(z = z - 18000)
  # Max P summary for CDFs
  get(paste0(model, '.marx.classified'))$jk %>%
  purrr::map_df(~.x$cdfP, .id = 'run') -> maxP
  # Max T summary for CDFs
  get(paste0(model, '.marx.classified'))$jk %>%
  purrr::map_df(~.x$cdfT, .id = 'run') -> maxT
  # Cluster centroids
  cent <- marx %>% slice(1) %>% ungroup() %>% count(class) %>%
  mutate(sumdP = mcl$parameters$mean[1,], maxP = mcl$parameters$mean[2,])
  # Summarise maxP and sumdP
  fts.summary <- marx %>%
  group_by(id, recovered, class) %>%
  summarise(sumdP = sum(diff(P)), maxP = max(P), .groups = 'keep')
  # Classification plots
  cat('\nDrawing ... [', model, ']', sep = '')
  p1 <- fts.summary %>%
  ggplot(aes(sumdP/1e4, maxP/1e4)) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = (median(fts.summary$sumdP) - thresh*IQR(fts.summary$sumdP))/1e4,
      ymin = -Inf,
      ymax = Inf),
    color = NA,
    fill = 'grey90',
    alpha = 0.1) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = (median(fts.summary$maxP) - thresh*IQR(fts.summary$maxP))/1e4),
    color = NA,
    fill = 'grey90',
    alpha = 0.1) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = (median(fts.summary$sumdP) - thresh*IQR(fts.summary$sumdP))/1e4,
      ymin = -Inf,
      ymax = (median(fts.summary$maxP) - thresh*IQR(fts.summary$maxP))/1e4),
    color = NA,
    fill = 'grey90',
    alpha = 0.1) +
  geom_point(aes(color = as.factor(class)), size = 0.5, shape = 15, show.legend = F) +
  geom_point(data = cent, aes(fill = as.factor(class), size = n), color = 'white', shape = 22) +
  labs(
    x = bquote('sumdP'~'[GPa]'),
    y = bquote('maxP'~'[GPa]'),
    fill = 'centroid'
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_size_continuous(range = c(2, 7), guide = 'none') +
  scale_color_manual(values = wesanderson::wes_palette(10, name = 'IsleofDogs1', type = 'continuous')) +
  scale_fill_manual(values = wesanderson::wes_palette(10, name = 'IsleofDogs1', type = 'continuous')) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_classic(base_size = 11)
  p2 <- fts.summary %>%
  ggplot(aes(sumdP/1e4, maxP/1e4)) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = (median(fts.summary$sumdP) - thresh*IQR(fts.summary$sumdP))/1e4,
      ymin = -Inf,
      ymax = Inf),
    color = NA,
    fill = 'grey90',
    alpha = 0.1) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = (median(fts.summary$maxP) - thresh*IQR(fts.summary$maxP))/1e4),
    color = NA,
    fill = 'grey90',
    alpha = 0.1) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = (median(fts.summary$sumdP) - thresh*IQR(fts.summary$sumdP))/1e4,
      ymin = -Inf,
      ymax = (median(fts.summary$maxP) - thresh*IQR(fts.summary$maxP))/1e4),
    color = NA,
    fill = 'grey90',
    alpha = 0.1) +
  geom_point(aes(color = recovered), shape = 15, size = 0.5) +
  scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'deeppink')) +
  scale_y_continuous(limits = c(0, NA)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(
    x = 'sumdP [GPa]',
    y = 'Probability'
  ) +
  theme_classic(base_size = 11)
  p <- p1 + (p2 +
    theme(
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
  )) +
  plot_annotation(title = paste0('Classification [', model, ']'), tag_levels = 'a') +
  plot_layout(guides = 'collect')
  cat('\nSaving classification plot [', model, ']', sep = '')
  ggsave(
    paste0('figs/k10/', model, '_class.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 4,
    dpi = 300
  )
  # ECDF plots
  p1 <- maxP %>%
  group_by(run) %>%
  ggplot() +
  geom_ribbon(
    data = maxP[maxP$cdf <= 0.8,],
    aes(x = maxP/1e4, ymin = 0, ymax = cdf, group = run),
    fill = 'grey50',
    alpha = 0.01) +
  geom_path(
    aes(x = maxP/1e4, y = cdf, linetype = 'markers', group = run),
    size = 0.1) +
  geom_ribbon(
    data = pd15[pd15$cumulative <= 0.8,],
    aes(ymin = 0, ymax = cumulative, x = pressure),
    alpha = 0.2) +
  geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
  labs(
    x = 'Maximum P [GPa]',
    y = 'Probability'
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
  theme_classic(base_size = 11)
# Marker maxT CDF
  p2 <- maxT %>%
  group_by(run) %>%
  ggplot() +
  geom_ribbon(
    data = maxT[maxT$cdf <= 0.8,],
    aes(x = maxT - 273, ymin = 0, ymax = cdf, group = run),
    fill = 'grey50',
    alpha = 0.01) +
  geom_path(
    aes(x = maxT - 273, y = cdf, linetype = 'markers', group = run),
    size = 0.1) +
  geom_ribbon(
    data = pd15T[pd15T$cdf <= 0.8,],
    aes(ymin = 0, ymax = cdf, x = T),
    alpha = 0.2) +
  geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
  labs(
    x = 'Maximum T [C]',
    y = 'Probability'
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
  theme_classic(base_size = 11)
  p <- p1 + (p2 +
    theme(
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
  )) +
  plot_annotation(title = paste0('Metamorphic conditions [', model, ']'), tag_levels = 'a') +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')
  cat('\nSaving metamorphic conditions plot [', model, ']', sep = '')
  ggsave(
    paste0('figs/k10/', model, '_meta.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 4,
    dpi = 300
  )

# Profile snapshots
  p1 <- grid %>%
  draw_grid(
    model = model,
    marx = marx,
    class = 'recovered',
    box = c(up = -18, down = 300, left = 500, right = 1800),
    time = tcut,
    bk.alpha = 0.5,
    mk.alpha = 0.3,
    mk.size = 0.25,
    iso.size = 2,
    leg.dir = 'horizontal',
    leg.dir.rec = 'horizontal',
    sub.col = 'deeppink',
    rec.col = 'black',
    leg.pos = 'bottom',
    p.type = 'viscosity',
    v.pal = 'viridis',
    base.size = 11,
    transparent = F)
  p2 <- g2 %>%
  draw_grid(
    model = model,
    marx = marx,
    class = 'recovered',
    box = c(up = -18, down = 300, left = 500, right = 1800),
    time = t2,
    bk.alpha = 0.5,
    mk.alpha = 0.3,
    mk.size = 0.25,
    iso.size = 2,
    leg.dir = 'horizontal',
    leg.dir.rec = 'horizontal',
    sub.col = 'deeppink',
    rec.col = 'black',
    leg.pos = 'bottom',
    p.type = 'viscosity',
    v.pal = 'viridis',
    base.size = 11,
    transparent = F)
  p3 <- g3 %>%
  draw_grid(
    model = model,
    marx = marx,
    class = 'recovered',
    box = c(up = -18, down = 300, left = 500, right = 1800),
    time = t3,
    bk.alpha = 0.5,
    mk.alpha = 0.3,
    mk.size = 0.25,
    iso.size = 2,
    leg.dir = 'horizontal',
    leg.dir.rec = 'horizontal',
    sub.col = 'deeppink',
    rec.col = 'black',
    leg.pos = 'bottom',
    p.type = 'viscosity',
    v.pal = 'viridis',
    base.size = 11,
    transparent = F)
  (p3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())) /
  (p2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())) /
  p1 +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') -> p
  cat('\nSaving snapshot plot [', model, ']', sep = '')
  ggsave(
    paste0('figs/k10/', model, '_snaps.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 7,
    dpi = 300
  )
  # Clean up environment
  rm(
    list = c(
      paste0(model, '.marx.classified'),
      paste0(model, '.grid'),
      paste0(model, '.marx')),
    envir = .GlobalEnv
  )
}

# Parallel computing
parallel::mcmapply(
  fun,
  models,
  paths,
  mc.cores = cores,
  SIMPLIFY = F
)

cat('\nDone!')