# Load functions and libraries
source('functions.R')

# Load classified markers
load('data/mods.RData')
load('data/marx_classified.RData')

# Model summary
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
  .x$mc$stats %>%
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

mods.summary %>% left_join(stats.summary)

purrr::map_df(
  marx.classified,
  ~marx_stats(.x$marx)$cdfP,
  .id = 'model'
) -> maxP

maxP %>%
ggplot(aes(x = maxP, y = cdf, color = model, group = model)) +
geom_path() +
labs(x = 'Maximum Pressure [GPa]', y = 'Probability', color = 'Model') +
theme_classic()
