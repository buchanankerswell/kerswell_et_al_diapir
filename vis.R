source('functions.R')

# Load data and glimpse its structure
# load('mods.RData')

print(mods)

# Explore features for individual models
print(mods$vis.oned[[1]])
print(mods$vis.oned[[22]])
# Summarise features by faceting
marx$mark.ft %>%
  f.oned(
    nrow = 4,
    ncol = 4,
    features = 'all',
    xlim = c(-2, 3),
    dlim = 0.3,
    alpha.min = 0.02
  ) %>%
  print()
marx$mark.ft %>%
  f.oned(
    nrow = 4,
    ncol = 4,
    features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
    xlim = c(-2, 3),
    dlim = 0.3,
    alpha.min = 0.02
  ) %>%
  print()
# visualize model parameters and results
# Plot BIC
mods$bic.mod0 %>% p.BIC()
mods$bic.mod1 %>% p.BIC()
mods$bic.mod2 %>% p.BIC()
mods$bic.mod3 %>% p.BIC()
# Plot scatter
# mods$mod0 %>% p.class()
# mods$mod1 %>% p.class()
# 1D facet plots
f.oned(
  mods$mark.ft,
  mods$mod0,
  nrow = 4,
  ncol = 4,
  features = features,
  plot = 'ridge',
  xlim = c(-2, 3),
  dlim = 0.3,
  alpha.min = 0
) %>%
  print()
f.oned(
  mods$mark.ft,
  mods$dr.mod0,
  mod.type = 'dr',
  dr = TRUE,
  nrow = 4,
  ncol = 4,
  features = features,
  plot = 'ridge',
  xlim = NULL,
  dlim = 0.3,
  alpha.min = 0
) %>%
  print()
f.oned(
  mods$mark.ft,
  mods$mod1,
  nrow = 4,
  ncol = 4,
  features = features,
  plot = 'ridge',
  xlim = c(-2, 3),
  dlim = 0.3,
  alpha.min = 0
) %>%
  print()
f.oned(
  mods$mark.ft,
  mods$dr.mod1,
  mod.type = 'dr',
  dr = TRUE,
  nrow = 4,
  ncol = 4,
  features = features,
  plot = 'ridge',
  xlim = NULL,
  dlim = 0.3,
  alpha.min = 0
) %>%
  print()
f.oned(
  mods$mark.ft,
  mods$mod2,
  nrow = 4,
  ncol = 4,
  features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
  plot = 'ridge',
  xlim = c(-2, 3),
  dlim = 0.3,
  alpha.min = 0
) %>%
  print()
f.oned(
  mods$mark.ft,
  mods$dr.mod2,
  mod.type = 'dr',
  dr = TRUE,
  nrow = 4,
  ncol = 4,
  features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
  plot = 'ridge',
  xlim = NULL,
  dlim = 0.3,
  alpha.min = 0
) %>%
  print()
# f.oned(
#   mods$mark.ft,
#   mods$mod3,
#   nrow = 4,
#   ncol = 4,
#   features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
#   plot = 'ridge',
#   xlim = c(-2, 3),
#   dlim = 0.3,
#   alpha.min = 0
# ) %>%
#   print()
# f.oned(
#   mods$mark.ft,
#   mods$dr.mod3,
#   mod.type = 'dr',
#   dr = TRUE,
#   nrow = 4,
#   ncol = 4,
#   features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
#   plot = 'ridge',
#   xlim = NULL,
#   dlim = 0.3,
#   alpha.min = 0
# ) %>%
#   print()
# Close PDF
dev.off()
# visualise fit
f.gif(
  mods$data.tidy,
  mods$mod0,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis.mod0'
)
f.gif(
  mods$data.tidy,
  mods$dr.mod0,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis.dr.mod0'
)
f.gif(
  mods$data.tidy,
  mods$mod1,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis.mod1'
)
f.gif(
  mods$data.tidy,
  mods$dr.mod1,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis.dr.mod1'
)
f.gif(
  mods$data.tidy,
  mods$mod2,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis.mod2'
)
f.gif(
  mods$data.tidy,
  mods$dr.mod2,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis.dr.mod2'
)
# f.gif(
#   mods$data.tidy,
#   mods$mod3,
#   runs = 'cde46',
#   GIF = 'PT',
#   ncol = 4,
#   nrow = 4,
#   grads = TRUE,
#   save = TRUE,
#   fname = 'vis.mod3'
# )
# f.gif(
#   mods$data.tidy,
#   mods$dr.mod3,
#   runs = 'cde46',
#   GIF = 'PT',
#   ncol = 4,
#   nrow = 4,
#   grads = TRUE,
#   save = TRUE,
#   fname = 'vis.dr.mod3'
# )