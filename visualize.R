# Load functions and libraries
source('functions.R')

# Read Penniston-Dorland et al., 2015 dataset
cat('\nReading Penniston-Dorland 2015 data/')
pd15 <- readr::read_delim('data/PD15.tsv', delim = '\t', col_types = 'cddcccd')
pd15 <- pd15[1:nrow(pd15)-1,]

tibble(T = sort(pd15$temperature)) %>%
mutate(cdf = (row_number()-1)/n()) -> pd15T

# Load marker and grid data
cat('\nReading RData files from data/')
paths <- list.files('data/k6', pattern = '.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

# Load markers and grids
path <- paths[23]
model <- models[23]
load_marx(paste0('data/', model, '_marx.RData'))
load(path)

# Marker data
marx <- get(paste0(model, '.marx.classified'))$marx
# Classification info
m <- get(paste0(model, '.marx.classified'))$gm
# Timecut
tcut <- attr(get(paste0(model, '.marx')), 'tcut')
# Nodes data
grid <- get(paste0(model, '.grid'))[[1]] %>% mutate(z = z - 18000)
# Max P summary for CDFs
get(paste0(model, '.marx.classified'))$mc %>%
purrr::map_df(~.x$cdfP, .id = 'run') -> maxP
# Max T summary for CDFs
get(paste0(model, '.marx.classified'))$mc %>%
purrr::map_df(~.x$cdfT, .id = 'run') -> maxT

# CDFs
p1 <- pd15 %>%
ggplot(aes(x = pressure, y = cumulative)) +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.8,],
  aes(ymin = 0, ymax = cumulative),
  alpha = 0.2) +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.9,],
  aes(ymin = 0, ymax = cumulative),
  alpha = 0.2) +
geom_path() +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
annotate('text', x = pd15$pressure[pd15$cumulative >= 0.8][1], y = 0.1, label = '80%', size = 3, hjust = 1.2) +
annotate('text', x = pd15$pressure[pd15$cumulative >= 0.9][1], y = 0.1, label = '90%', size = 3, hjust = 0.5) +
labs(x = 'Pressure [GPa]', y = 'Probability') +
theme_classic(base_size = 8)
p2 <- pd15T %>%
ggplot(aes(x = T, y = cdf)) +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.8,],
  aes(ymin = 0, ymax = cdf),
  alpha = 0.2) +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.9,],
  aes(ymin = 0, ymax = cdf),
  alpha = 0.2) +
geom_path() +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
annotate('text', x = pd15T$T[pd15T$cdf >= 0.8][1], y = 0.1, label = '80%', size = 3, hjust = 1.2) +
annotate('text', x = pd15T$T[pd15T$cdf >= 0.9][1], y = 0.1, label = '90%', size = 3, hjust = 0.5) +
labs(x = 'Temperature [C]', y = 'Probability') +
theme_classic(base_size = 8)
p <- p1 + p2 +
plot_annotation(tag_levels = 'a')
ggsave(
  paste0('figs/pd15_cdf.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3,
  dpi = 300
)

# Initial conditions
p <- grid %>%
draw_grid(
  model = model,
  marx = marx,
  class = 'type',
  time = 1,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  bk.alpha = 0.3,
  iso.alpha = 0.3,
  iso.col = 'black',
  mk.size = 0.25,
  sub.col = 'deeppink',
  rec.col = 'white',
  leg.pos = 'bottom',
  leg.dir = 'horizontal',
  leg.dir.rec = 'horizontal',
  base.size = 8,
  p.type = 'viscosity',
  v.pal = 'viridis',
  transparent = F) +
theme(panel.background = element_rect(color = 'black'), plot.margin = margin(0, 0.5, 0, 0, 'lines')) +
annotate('rect', xmin = 500, xmax = 1260, ymin = 0, ymax = 11, fill = NA, color = rgb(0,0,0,0.7), size = 0.3) +
annotate('curve', x = 750, xend = 1000, y = 150, yend = 11, curvature = 0.3, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4) +
annotate('text', x = 750, y = 150, label = 'traced markers', hjust = 1, size = 3) +
annotate('text', x = 250, y = 10, label = 'free surface', vjust = 0.2, size = 3) +
annotate('text', x = 1200, y = Inf, label = 'open boundary', vjust = -0.2, size = 3) +
annotate('text', x = -Inf, y = Inf, label = 'free slip', angle = 90, vjust = 1.2, hjust = -0.1, size = 3) +
annotate('text', x = Inf, y = Inf, label = 'free slip', angle = 90, vjust = -0.2, hjust = -0.1, size = 3) +
#
annotate('rect', xmin = 490, xmax = 510, ymin = -4, ymax = 40, fill = NA, color = 'black', size = 0.3) +
annotate('rect', xmin = 1790, xmax = 1810, ymin = -4, ymax = 40, fill = NA, color = 'black', size = 0.3) +
annotate('segment', x = 510, xend = 560, y = 30, yend = 30, arrow = arrow(length = unit(0.05, 'in'), angle = 20), color = 'black', lineend = 'round', linejoin = 'round', size = 1) +
annotate('text', x = 400, y = 150, label = 'convergence region', hjust = 1, size = 3) +
annotate('curve', x = 405, xend = 500, y = 150, yend = 40, curvature = 0.5, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4) +
annotate('text', x = 1800, y = 150, label = 'convergence region', hjust = 1, size = 3) +
annotate('curve', x = 1800, xend = 1800, y = 150, yend = 40, curvature = 0.5, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4) +
annotate('segment', x = 1790, xend = 1750, y = 30, yend = 30, arrow = arrow(length = unit(0.05, 'in'), angle = 20), color = rgb(0, 0, 0, 0.3), lineend = 'round', linejoin = 'round', size = 1) +
annotate('text', x = 1400, y = 100, label = 'trench', hjust = 0, size = 3) +
annotate('curve', x = 1400, xend = 1260, y = 100, yend = 0, curvature = 0.2, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4)
ggsave(
  paste0('figs/', model, '_init.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 1.8,
  dpi = 300
)