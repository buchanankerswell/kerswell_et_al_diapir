# Load functions and libraries
source('functions.R')

# Read Penniston-Dorland et al., 2015 dataset
cat('\nReading Penniston-Dorland 2015 ...')
pd15 <- readr::read_delim('data/PD15.tsv', delim = '\t', col_types = 'cddcccd')
pd15 <- pd15[1:nrow(pd15)-1,]

tibble(T = sort(pd15$temperature)) %>%
mutate(cdf = (row_number()-1)/n()) -> pd15T

# Load numerical model parameters
load('data/mods.RData')

# Get file paths
cat('\nReading classified marx ...')
paths <- list.files('data/k10', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

# Load classified markers
for (i in paths) load(i)

# Save as list
purrr::map(ls()[grep('classified', ls())], ~get(.x)) %>%
purrr::set_names(models) -> m

# Clean up environment
rm(list = ls()[grep('classified', ls())])

# Time cuts
tcuts <- purrr::map_dbl(m, ~attr(.x$marx, 'tcut'))

# Number of markers by model
purrr::map_df(m, ~{
  .x$marx %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = n())
}, .id = 'model') -> marx.summary

# Combine marx and model data
mods.summary <- mods %>%
select(model, phi, zc, z1100, age, cv) %>%
left_join(marx.summary, by = 'model') %>%
mutate(tcut = tcuts, .before = 'phi')

# Summarise marker stats by model
purrr::map_df(m, ~{
  .x$jk %>%
  purrr::map_df(~purrr::pluck(.x, 'stats'), .id = 'run') %>%
  summarise(
    mean.rec = mean(rec),
    sd.rec = sd(rec),
    mean.sub = mean(sub),
    sd.sub = sd(sub),
    mean.ratio = mean(ratio),
    sd.ratio = sd(ratio)
  )
}, .id = 'model') -> stats.summary

# Combine tables
d <- mods.summary %>% left_join(stats.summary, by = 'model')

# PD15 CDFs
p1 <- pd15 %>%
ggplot(aes(x = pressure, y = cumulative)) +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.8,],
  aes(ymin = 0, ymax = cumulative),
  alpha = 0.2) +
geom_path() +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
annotate('text', x = pd15$pressure[pd15$cumulative >= 0.8][1], y = 0.1, label = '80%', size = 3, hjust = 1.2) +
labs(x = 'Pressure [GPa]', y = 'Probability') +
theme_classic(base_size = 11)
p2 <- pd15T %>%
ggplot(aes(x = T, y = cdf)) +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.8,],
  aes(ymin = 0, ymax = cdf),
  alpha = 0.2) +
geom_path() +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
annotate('text', x = pd15T$T[pd15T$cdf >= 0.8][1], y = 0.1, label = '80%', size = 3, hjust = 1.2) +
labs(x = 'Temperature [C]', y = 'Probability') +
theme_classic(base_size = 11)
p <- p1 +
(p2 +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
)) +
plot_annotation(tag_levels = 'a')
cat('\nSaving PD15 cdf plot')
ggsave(
  paste0('figs/pd15_cdf.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3,
  dpi = 300
)

# Summarise pressure CDF
purrr::map_df(m, ~{
  .x$marx %>%
  filter(recovered == TRUE)
}, .id = 'model') %>%
group_by(id, model) %>%
summarise(maxP = max(P), .groups = 'keep') %>%
left_join(select(mods, model, z1100), by = 'model', .before = 'data') %>%
group_by(z1100) %>%
nest() -> d1

maxP <- d1$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxP) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

# Summarise temperature CDF
purrr::map_df(m, ~{
  .x$marx %>%
  filter(recovered == TRUE)
}, .id = 'model') %>%
group_by(id, model) %>%
summarise(maxT = max(T), .groups = 'keep') %>%
left_join(select(mods, model, z1100), by = 'model', .before = 'data') %>%
group_by(z1100) %>%
nest() -> d2

maxT <- d2$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxT) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

# Biased plots
maxP.bias10 <- d1$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxP) %>%
  slice_head(n = nrow(.x)*0.9) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

maxT.bias10 <- d2$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxT) %>%
  slice_head(n = nrow(.x)*0.9) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

# CDFS summary by lithospheric thickness
p1 <- maxP.bias10 %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.8,],
  aes(ymin = 0, ymax = cumulative, x = pressure),
  alpha = 0.2) +
geom_ribbon(
  data = maxP.bias10[maxP.bias10$cdf <= 0.8,],
  aes(x = maxP/1e4, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxP/1e4, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_x_continuous(breaks = seq(1, 10, 2)) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum P [GPa]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
p2 <- maxT.bias10 %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.8,],
  aes(ymin = 0, ymax = cdf, x = T),
  alpha = 0.2) +
geom_ribbon(
  data = maxT.bias10[maxT.bias10$cdf <= 0.8,],
  aes(x = maxT - 273, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxT - 273, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum T [C]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
pp1 <- p1 + (p2 +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
)) +
plot_annotation(title = 'PT conditions [all models 10% bias]', tag_levels = 'a') +
plot_layout(guides = 'collect') &
theme(legend.position = 'bottom')
cat('\nSaving ECDF bias 10 plot')
ggsave(
  'figs/meta_bias10.png',
  plot = pp1,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 4,
  dpi = 300
)

# Biased plots
maxP.bias20 <- d1$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxP) %>%
  slice_head(n = nrow(.x)*0.7) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

maxT.bias20 <- d2$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxT) %>%
  slice_head(n = nrow(.x)*0.7) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

# CDFS summary by lithospheric thickness
p3 <- maxP.bias20 %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.8,],
  aes(ymin = 0, ymax = cumulative, x = pressure),
  alpha = 0.2) +
geom_ribbon(
  data = maxP.bias20[maxP.bias20$cdf <= 0.8,],
  aes(x = maxP/1e4, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxP/1e4, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_x_continuous(breaks = seq(1, 10, 2)) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum P [GPa]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
p4 <- maxT.bias20 %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.8,],
  aes(ymin = 0, ymax = cdf, x = T),
  alpha = 0.2) +
geom_ribbon(
  data = maxT.bias20[maxT.bias20$cdf <= 0.8,],
  aes(x = maxT - 273, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxT - 273, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum T [C]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
pp2 <- p3 + (p4 +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
)) +
plot_annotation(title = 'PT conditions [all models 20% bias]', tag_levels = 'a') +
plot_layout(guides = 'collect') &
theme(legend.position = 'bottom')
cat('\nSaving ECDF bias 20 plot')
ggsave(
  'figs/meta_bias20.png',
  plot = pp2,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 4,
  dpi = 300
)

# Biased plots
maxP.bias30 <- d1$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxP) %>%
  slice_head(n = nrow(.x)*0.7) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

maxT.bias30 <- d2$data %>%
purrr::map_df(
  ~.x %>%
  arrange(maxT) %>%
  slice_head(n = nrow(.x)*0.7) %>%
  mutate(cdf = (row_number()-1)/n())
) %>%
arrange(model) %>%
left_join(select(mods, model, z1100), by = 'model')

# CDFS summary by lithospheric thickness
p3 <- maxP.bias30 %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.8,],
  aes(ymin = 0, ymax = cumulative, x = pressure),
  alpha = 0.2) +
geom_ribbon(
  data = maxP.bias30[maxP.bias30$cdf <= 0.8,],
  aes(x = maxP/1e4, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxP/1e4, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_x_continuous(breaks = seq(1, 10, 2)) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum P [GPa]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
p4 <- maxT.bias30 %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.8,],
  aes(ymin = 0, ymax = cdf, x = T),
  alpha = 0.2) +
geom_ribbon(
  data = maxT.bias30[maxT.bias30$cdf <= 0.8,],
  aes(x = maxT - 273, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxT - 273, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum T [C]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
pp2 <- p3 + (p4 +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
)) +
plot_annotation(title = 'PT conditions [all models 30% bias]', tag_levels = 'a') +
plot_layout(guides = 'collect') &
theme(legend.position = 'bottom')
cat('\nSaving ECDF bias 30 plot')
ggsave(
  'figs/meta_bias30.png',
  plot = pp2,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 4,
  dpi = 300
)

# CDFS no bias
p1 <- maxP %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15[pd15$cumulative <= 0.8,],
  aes(ymin = 0, ymax = cumulative, x = pressure),
  alpha = 0.2) +
geom_ribbon(
  data = maxP[maxP$cdf <= 0.8,],
  aes(x = maxP/1e4, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxP/1e4, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15, aes(x = pressure, y = cumulative, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_x_continuous(breaks = seq(1, 10, 2)) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum P [GPa]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
p2 <- maxT %>%
group_by(z1100) %>%
ggplot() +
geom_ribbon(
  data = pd15T[pd15T$cdf <= 0.8,],
  aes(ymin = 0, ymax = cdf, x = T),
  alpha = 0.2) +
geom_ribbon(
  data = maxT[maxT$cdf <= 0.8,],
  aes(x = maxT - 273, ymin = 0, ymax = cdf, fill = as.factor(z1100), group = z1100),
  alpha = 0.2) +
geom_line(aes(x = maxT - 273, y = cdf, linetype = 'markers', color = as.factor(z1100), group = z1100), show.legend = F) +
geom_path(data = pd15T, aes(x = T, y = cdf, linetype = 'PD15')) +
guides(fill = guide_legend(override.aes = list(alpha = 1))) +
scale_linetype_manual(name = NULL, values = c('PD15' = 'dotted', 'markers' = 'solid')) +
scale_y_continuous(breaks = seq(0, 1, 0.2)) +
scale_color_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
scale_fill_manual(values = wesanderson::wes_palette(4, name = 'IsleofDogs1')) +
labs(
  x = 'Maximum T [C]',
  y = 'Probability',
  color = bquote(z[1100]~'[km]'),
  fill = bquote(z[1100]~'[km]')
) +
theme_classic(base_size = 11)
p <- p1 + (p2 +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
)) +
plot_annotation(title = 'PT conditions [all models 0% bias]', tag_levels = 'a') +
plot_layout(guides = 'collect') &
theme(legend.position = 'bottom')
cat('\nSaving ECDF plot')
ggsave(
  'figs/meta_all.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 4,
  dpi = 300
)

# Load markers and grids for model cde78
path <- paths[19]
model <- models[19]
load_marx(paste0('/Volumes/hd/nmods/kerswell_et_al_marx/data/', model, '_marx.RData'))
load(path)
# Marker data
marx <- get(paste0(model, '.marx.classified'))$marx
# Classification info
mcl <- get(paste0(model, '.marx.classified'))$mcl
# Timecut
tcut <- attr(get(paste0(model, '.marx')), 'tcut')
# Nodes data
grid <- get(paste0(model, '.grid'))[[1]] %>% mutate(z = z - 18000)
# Max P summary for CDFs
get(paste0(model, '.marx.classified'))$jk %>%
purrr::map_df(~.x$cdfP, .id = 'run') -> maxP
# Max T summary for CDFs
get(paste0(model, '.marx.classified'))$jk %>%
purrr::map_df(~.x$cdfT, .id = 'run') -> maxT
# Initial conditions
p <- grid %>%
draw_grid(
  model = model,
  marx = marx,
  class = 'type',
  time = 1,
  box = c(up = -18, down = 300, left = 0, right = 2000),
  bk.alpha = 0.3,
  iso.alpha = 0.3,
  iso.col = 'black',
  iso.skip = 1,
  mk.size = 0.25,
  sub.col = 'deeppink',
  rec.col = 'white',
  leg.pos = 'bottom',
  leg.dir = 'horizontal',
  leg.dir.rec = 'horizontal',
  base.size = 11,
  p.type = 'viscosity',
  v.pal = 'viridis',
  transparent = F) +
theme(panel.background = element_rect(color = 'black'), plot.margin = margin(0, 0.65, 0, 0, 'lines'), plot.title = element_blank()) +
annotate('rect', xmin = 500, xmax = 1260, ymin = 0, ymax = 11, fill = NA, color = rgb(0,0,0,0.7), size = 0.3) +
annotate('curve', x = 750, xend = 1000, y = 150, yend = 11, curvature = 0.3, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4) +
annotate('text', x = 750, y = 150, label = 'traced markers', hjust = 1, size = 3) +
annotate('text', x = 250, y = 10, label = 'free surface', vjust = 0.2, size = 3) +
annotate('text', x = 1200, y = Inf, label = 'open boundary', vjust = -0.2, size = 3) +
annotate('text', x = -Inf, y = Inf, label = 'free slip', angle = 90, vjust = 1.2, hjust = -0.1, size = 3) +
annotate('text', x = Inf, y = Inf, label = 'free slip', angle = 90, vjust = -0.2, hjust = -0.1, size = 3) +
annotate('rect', xmin = 490, xmax = 510, ymin = -4, ymax = 40, fill = NA, color = 'black', size = 0.3) +
annotate('rect', xmin = 1790, xmax = 1810, ymin = -4, ymax = 40, fill = NA, color = 'black', size = 0.3) +
annotate('segment', x = 510, xend = 560, y = 30, yend = 30, arrow = arrow(length = unit(0.05, 'in'), angle = 20), color = 'black', lineend = 'round', linejoin = 'round', size = 1) +
annotate('text', x = 400, y = 150, label = 'convergence region', hjust = 1, size = 3) +
annotate('curve', x = 405, xend = 500, y = 150, yend = 40, curvature = 0.5, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4) +
annotate('text', x = 1800, y = 150, label = 'convergence region', hjust = 1, size = 3) +
annotate('curve', x = 1800, xend = 1800, y = 150, yend = 40, curvature = 0.5, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4) +
annotate('segment', x = 1790, xend = 1750, y = 30, yend = 30, arrow = arrow(length = unit(0.05, 'in'), angle = 20), color = 'black', lineend = 'round', linejoin = 'round', size = 1) +
annotate('text', x = 1400, y = 100, label = 'trench', hjust = 0, size = 3) +
annotate('curve', x = 1400, xend = 1260, y = 100, yend = 0, curvature = 0.2, arrow = arrow(length = unit(0.05, 'in'), angle = 20), lineend = 'round', size = 0.4)
cat('\nSaving initial conditions plot')
ggsave(
  paste0('figs/', model, '_init.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 2.3,
  dpi = 300
)

# Clean up environment
rm(list = ls()[grep('classified', ls())])

# Diapir examples
path <- paths[c(3, 4, 8, 20, 24)]
model <- models[c(3, 4, 8, 20, 24)]

cat('\nPlotting diapir examples ...')
pp <- purrr::map2(path, model, ~{
  # Load markers and grids
  load_marx(paste0('/Volumes/hd/nmods/kerswell_et_al_marx/data/', .y, '_marx.RData'))
  load(.x)
  # Marker data
  marx <- get(paste0(.y, '.marx.classified'))$marx
  # Classification info
  mcl <- get(paste0(.y, '.marx.classified'))$mcl
  # Timecut
  tcut <- attr(get(paste0(.y, '.marx')), 'tcut')
  # Nodes data
  grid <- get(paste0(.y, '.grid'))[[tcut]] %>% mutate(z = z - 18000)
  # Draw profiles
  p <- grid %>%
  draw_grid(
    model = .y,
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
  # clean up environment
  rm(list = c(paste0(.y, '.marx.classified'), paste0(.y, '.marx'), paste0(.y, '.grid')), envir = .GlobalEnv)
  # return plot
  return(p)
}) %>% purrr::set_names(model)
# Draw diapir plot
(pp[[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) /
(pp[[2]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) /
(pp[[3]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) /
(pp[[4]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) /
pp[[5]] +
plot_annotation(tag_levels = 'a') +
plot_layout(guides = 'collect') &
theme(legend.position = 'bottom') -> p
# Save
cat('\nSaving diapir plot')
ggsave(
  paste0('figs/diapirs.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 11.8,
  dpi = 300
)

# Type I error examples
path <- paths[49]
model <- models[49]

cat('\nPlotting type I error example ...')
# Load markers and grids
load_marx(paste0('/Volumes/hd/nmods/kerswell_et_al_marx/data/', model, '_marx.RData'))
load(path)
# Marker data
marx <- get(paste0(model, '.marx.classified'))$marx
# Classification info
mcl <- get(paste0(model, '.marx.classified'))$mcl
# Timecut
tcut <- attr(get(paste0(model, '.marx')), 'tcut')
# Nodes data
grid <- get(paste0(model, '.grid'))[[tcut]] %>% mutate(z = z - 18000)
# Draw profiles
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
# Max P summary for CDFs
get(paste0(model, '.marx.classified'))$jk %>%
purrr::map_df(~.x$cdfP, .id = 'run') -> maxP
# Max T summary for CDFs
get(paste0(model, '.marx.classified'))$jk %>%
purrr::map_df(~.x$cdfT, .id = 'run') -> maxT
# Marker maxP CDF
p2 <- maxP %>%
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
p3 <- maxT %>%
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
p <- p1 /
(p2 + (p3 + theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
))) +
plot_annotation(tag_levels = 'a') +
plot_layout(guides = 'collect', heights = c(1, 1.5)) &
theme(legend.position = 'bottom')
# Save
cat('\nSaving typeI plot')
ggsave(
  paste0('figs/typeI.png'),
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 6.8,
  dpi = 300
)

# Clean up environment
rm(list = ls()[grep('classified', ls())])

cat('\nDone!')
