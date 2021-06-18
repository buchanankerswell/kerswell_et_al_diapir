# Load functions and libraries
source('functions.R')

# Load numerical model parameters
load('data/mods.RData')

# Get file paths
paths <- list.files('data/k6', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')
cat('\nFound models:', models, sep = '\n')

# Load classified markers
for (i in paths) load(i)

# Save as list
purrr::map(ls()[grep('classified', ls())], ~get(.x)) %>%
purrr::set_names(models) -> m

rm(list = ls()[grep('classified', ls())])

# Number of markers by model
purrr::map_df(m, ~{
  .x$marx %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = n())
}, .id = 'model') -> marx.summary

mods.summary <- mods %>%
select(model, phi, zc, z1100, age, cv) %>%
left_join(marx.summary, by = 'model')

# Summarise marker stats by model
purrr::map_df(m, ~{
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

d %>%
ggplot(aes(zc, med.max.P.rec/1e4)) +
geom_point() +
labs(x = 'Coupling Depth [km]', y = 'Max P [GPa]') +
facet_grid(rows = vars(age), cols = vars(cv)) +
theme_grey()

d %>%
ggplot(aes(zc, med.max.P.rec/1e4, color = z1100)) +
geom_point() +
labs(x = 'Coupling Depth [km]', y = 'Max P [GPa]') +
theme_grey()

d %>%
ggplot(aes(z1100, mean.max.P.rec)) +
geom_point()

d %>%
ggplot(aes(zc, mean.ratio)) +
geom_point()

d %>%
ggplot(aes(phi, mean.ratio)) +
geom_point()

d %>%
ggplot(aes(as.factor(age), mean.ratio)) +
geom_boxplot()

d %>%
ggplot(aes(as.factor(cv), mean.ratio)) +
geom_boxplot()

# Read Penniston-Dorland et al., 2015 dataset
pd15 <- readr::read_delim('data/PD15.tsv', delim = '\t', col_types = 'cddcccd')
pd15 <- pd15[1:nrow(pd15)-1,]

# cdfP
purrr::map_df(m, ~{
  .x$marx %>%
  filter(recovered == TRUE) %>%
  summarise(maxP = max(P)) %>%
  arrange(maxP) %>%
  mutate(cdf = (row_number()-1)/n())
}, .id = 'model') -> cdfP

pd <- pd15 %>%
select(pressure, cumulative) %>%
rename(maxP = pressure, cdf = cumulative)

mods %>%
select(zc) %>%
right_join(cdfP, by = 'model') %>%
ggplot(aes(x = maxP/1e4, y = cdf, color = zc, group = model)) +
geom_path(show.legend = F) +
labs(
  x = 'Max P [GPa]',
  y = 'Probability',
  color = bquote(z[1100]),
  title = paste0('Cumulative probability of max P')
) +
theme_classic() +
theme()

# cdfT
purrr::map_df(m, ~{
  .x$marx %>%
  filter(recovered == TRUE) %>%
  summarise(maxT = max(T)) %>%
  arrange(maxT) %>%
  mutate(cdf = (row_number()-1)/n())
}, .id = 'model') -> cdfT

mods %>%
select(z1100) %>%
right_join(cdfT, by = 'model') %>%
ggplot(aes(x = maxT - 273, y = cdf, color = as.factor(z1100), group = model)) +
geom_path() +
labs(
  x = 'Max T [C]',
  y = 'Probability',
  color = bquote(z[1100]),
  title = paste0('Cumulative probability of max T')
) +
scale_color_grey() +
theme_classic()
