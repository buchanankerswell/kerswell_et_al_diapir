source('functions.R')

paths <- list.files('data', pattern = '_marx.RData', full.names = T)
models <- paths %>% stringr::str_extract('cd.[0-9]+')

# Summarise marker features
purrr::map2(models, paths, ~{
  # Load markers
  load_marx(.y)
})

df <- cdf78.marx

# Compute features
df %>%
summarise(
  tsteps = n(),
  under.three.kbar = sum(P < 3e3),
  under.ten.kbar = sum(P < 1e4),
  above.thirty.kbar = sum(P > 3e4),
  under.one.hunderd.c = sum(T < 373),
  under.five.hundred.c = sum(T < 773),
  above.seven.hundred.c = sum(T > 973),
  max.P = max(P),
  med.P = median(P),
  mean.P = mean(P),
  iqr.P = IQR(P),
  max.T = max(T),
  med.T = median(T),
  mean.T = mean(T),
  iqr.T = IQR(T),
  up.dP = sum(diff(P) > 0),
  down.dP = sum(diff(P) < 0),
  runup.dP = {rn <- rle(diff(P) > 0); rn$lengths[which(rn$values == TRUE)] %>% max()},
  rundown.dP = {rn <- rle(diff(P) > 0); rn$lengths[which(rn$values == FALSE)] %>% max()},
  sum.dP = sum(diff(P)),
  sumup.dP = sum(diff(P)[which(diff(P) > 0)]),
  sumdown.dP = sum(diff(P)[which(diff(P) < 0)]),
  up.dT = sum(diff(T) > 0),
  down.dT = sum(diff(T) < 0),
  runup.dT = {rn <- rle(diff(T) > 0); rn$lengths[which(rn$values == TRUE)] %>% max()},
  rundown.dT = {rn <- rle(diff(T) > 0); rn$lengths[which(rn$values == FALSE)] %>% max()},
  sum.dT = sum(diff(T)),
  sumup.dT = sum(diff(T)[which(diff(T) > 0)]),
  sumdown.dT = sum(diff(T)[which(diff(T) < 0)]),
  .groups = 'keep') -> d
# Clean up infinites
d[d == -Inf] <- 0
d[d == Inf] <- 0

d %>% GGally::ggpairs() -> p
ggsave(filename = 'figs/features/correlation.png', plot = p, width = 48, height = 48, units = 'in', dpi = 72)