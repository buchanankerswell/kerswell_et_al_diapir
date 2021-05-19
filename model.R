# Load functions and libraries
source('functions.R')
# Load marker features
load('data/marx_features.RData')
load_marx('data/cdf78_marx.RData')
# Eigenvalue Decomposition Discriminant Analysis (EDDA)
features
ft <- c(
  'tsteps',
  'under.five.hundred.c',
  'above.thirty.kbar',
  'above.seven.hundred.c',
  'under.three.kbar',
  'under.one.hundred.c',
  'under.ten.kbar',
  'max.P',
  'med.P',
  'iqr.T',
  'up.dP',
  'down.dP',
  'sum.dP'
)
# Save IDs
ids <- marx.features[[1]]$id
# Select features
marx.features[[1]] %>%
ungroup() %>%
select(ft) -> X
# Gaussian mixture clustering
mc.bic <- mclustBIC(X, G = seq_len(10))
summary(mc.bic)
plot(mc.bic)
# Clustering
mc <- Mclust(X, verbose = T)
summary(mc)
plot(mc, what = 'classification')
# Dimension reduction
mc.dr <- MclustDR(mc, lambda = 1)
summary(mc.dr)
plot(mc.dr, dimens = 1, what = 'density')
plot(mc.dr, dimens = 2, what = 'density')
plot(mc.dr, dimens = c(1,2), what = 'scatterplot')
# Join classes to markers data
df <- cdf78.marx %>% left_join(tibble(id = ids, class = mc$classification), by = c('id'))
# Animate
marx_PT_mov(df, 'cdf78', class = T)
marx_motion_mov(df, 'cdf78', class = T)
