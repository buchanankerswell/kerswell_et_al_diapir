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

# Draw cross section at timecut
p1 <- grid %>%
draw_grid(
  model = model,
  marx = marx,
  class = 'type',
  time = 1,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  bk.alpha = 0.4,
  mk.size = 0.6,
  sub.col = 'deeppink',
  rec.col = 'white',
  leg.pos = 'bottom',
  leg.dir = 'horizontal',
  leg.dir.rec = 'horizontal',
  p.type = 'blank',
  v.pal = 'viridis',
  transparent = F) +
theme(panel.background = element_rect(color = 'black')) +
annotate('text', x = 250, y = 10, label = 'free surface', vjust = 0.2, size = 3) +
annotate('text', x = 1200, y = Inf, label = 'open boundary', vjust = -0.2, size = 3) +
annotate('text', x = -Inf, y = Inf, label = 'free slip', angle = 90, vjust = 1.2, hjust = -0.1, size = 3) +
annotate('text', x = Inf, y = Inf, label = 'free slip', angle = 90, vjust = -0.2, hjust = -0.1, size = 3)
p2 <- grid %>%
draw_grid(
  model = model,
  time = 1,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  bk.alpha = 0.9,
  mk.size = 0.6,
  sub.col = 'deeppink',
  rec.col = 'white',
  leg.pos = 'bottom',
  leg.dir = 'horizontal',
  leg.dir.rec = 'horizontal',
  p.type = 'viscosity',
  v.pal = 'viridis',
  transparent = F) +
annotate('rect', xmin = 490, xmax = 510, ymin = -4, ymax = 40, fill = NA, color = 'black', size = 0.3) +
annotate('rect', xmin = 1790, xmax = 1810, ymin = -4, ymax = 40, fill = NA, color = 'black', size = 0.3) +
annotate('segment', x = 510, xend = 560, y = 20, yend = 20, arrow = arrow(length = unit(0.1, 'in'), angle = 20), color = 'black', lineend = 'round', linejoin = 'round', size = 0.3) +
annotate('text', x = 400, y = 150, label = 'convergence region', hjust = 1, size = 3) +
annotate('curve', x = 405, xend = 500, y = 150, yend = 40, curvature = 0.5, arrow = arrow(length = unit(0.1, 'in'), angle = 20), lineend = 'round', size = 0.3) +
annotate('text', x = 1700, y = 150, label = 'weak zone', hjust = 0, size = 3) +
annotate('curve', x = 1695, xend = 1490, y = 150, yend = 75, curvature = -0.2, arrow = arrow(length = unit(0.1, 'in'), angle = 20), lineend = 'round', size = 0.3)
(p1 + theme(axis.title.x = element_blank())) /
p2  + theme(axis.title.x = element_blank()) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
theme(legend.position = 'bottom') -> p
ggsave(
  paste0('figs/', model, '_init.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7.25,
  height = 4,
  dpi = 300
)

marx <- get(paste0(model, '.marx'))
grid <- get(paste0(model, '.grid'))[[tcut+10]] %>% mutate(z = z - 18000)

p1 <- grid %>%
draw_grid(
  model = model,
  marx = marx,
  class = 'type',
  time = tcut+10,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  bk.alpha = 0.4,
  mk.size = 0.6,
  sub.col = 'deeppink',
  rec.col = 'white',
  leg.pos = 'bottom',
  leg.dir = 'horizontal',
  leg.dir.rec = 'horizontal',
  p.type = 'blank',
  v.pal = 'viridis',
  transparent = F) +
theme(panel.background = element_rect(color = 'black'))
p2 <- grid %>%
draw_grid(
  model = model,
  time = tcut+10,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  bk.alpha = 0.9,
  mk.size = 0.6,
  sub.col = 'deeppink',
  rec.col = 'white',
  leg.pos = 'bottom',
  leg.dir = 'horizontal',
  leg.dir.rec = 'horizontal',
  p.type = 'viscosity',
  v.pal = 'viridis',
  transparent = F)
(p1 + theme(axis.title.x = element_blank())) /
p2  + theme(axis.title.x = element_blank()) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
theme(legend.position = 'bottom') -> p
ggsave(
  paste0('figs/', model, '_interference.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7.25,
  height = 4,
  dpi = 300
)

p1 <- pd15 %>%
ggplot(aes(x = pressure, y = cumulative)) +
geom_path() +
labs(
  title = 'Cumulative probability of max P',
  x = 'Maximum P [GPa]',
  y = 'Probability'
) +
theme_classic()
p2 <- pd15T %>%
ggplot(aes(x = T, y = cdf)) +
geom_path() +
labs(
  title = 'Cumulative probability of max T',
  x = 'Maximum T [C]',
  y = 'Probability'
) +
theme_classic()
p1 + p2 +
plot_annotation(tag_levels = 'a', caption = 'from Penniston-Dorland et al. (2015)')
ggsave(
  paste0('figs/pd15_cdfP.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 5,
  height = 5,
  dpi = 300
)