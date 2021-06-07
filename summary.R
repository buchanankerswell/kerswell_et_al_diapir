# Load functions and libraries
source('functions.R')

# Load numerical model parameters
load('data/mods.RData')

# Get file paths
paths <- list.files('data/k6_classd', full.names = T)
names <- paths %>% stringr::str_extract('cd.[0-9]+')
marx.files <- list.files('data', pattern = '_marx.RData') %>% stringr::str_extract('cd.[0-9]+')
cat('\nFound models:', names, sep = '\n')

# Load classified markers
for (i in paths[names %in% marx.files]) load(i)

# Save as list
purrr::map(ls()[grep('classified', ls())], ~get(.x)) %>%
purrr::set_names(names[names %in% marx.files]) -> marx.classified

# Number of markers by model
purrr::map_df(marx.classified, ~{
  .x$marx %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = n())
}, .id = 'model') -> marx.summary

mods.summary <- mods %>%
select(model, phi, zc, z1100, age, cv) %>%
left_join(marx.summary, by = 'model')

# Summarise marker stats by model
purrr::map_df(marx.classified, ~{
  .x$mc %>%
  purrr::map_df(~purrr::pluck(.x, 'stats'), .id = 'run') %>%
  summarise(
    mean.rec = mean(recovered),
    sd.rec = sd(recovered),
    med.rec = median(recovered),
    iqr.rec = IQR(recovered),
    mean.sub = mean(subducted),
    sd.sub = sd(subducted),
    med.sub = median(subducted),
    iqr.sub = IQR(subducted),
    mean.ratio = mean(ratio),
    sd.ratio = sd(ratio),
    med.ratio = median(ratio),
    iqr.ratio = IQR(ratio),
    mean.max.P.rec = mean(max.P.rec),
    sd.max.P.rec = sd(max.P.rec),
    med.max.P.rec = median(max.P.rec),
    iqr.max.P.rec = IQR(max.P.rec),
    mean.max.T.rec = mean(max.T.rec),
    sd.max.T.rec = sd(max.T.rec),
    med.max.T.rec = median(max.T.rec),
    iqr.max.T.rec = IQR(max.T.rec)
  )
}, .id = 'model') -> stats.summary

# Combine tables
d <- mods.summary %>% left_join(stats.summary)

p <- d %>%
ungroup() %>%
select(-model, mean.rec, mean.sub, mean.ratio, mean.max.P.rec, mean.max.T.rec) %>%
GGally::ggpairs()
ggsave(plot = p, file = 'figs/corr.png', width = 24, height = 24, dpi = 300, device = 'png', type = 'cairo')
