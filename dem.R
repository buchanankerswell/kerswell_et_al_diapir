source('functions.R')

# Open PDF for plotting
cairo_pdf('vis2.PDF',
          width = 11,
          height = 8.5,
          onefile = TRUE)

# Load data and glimpse its structure
# load('marx.RData')

marx
glimpse(marx)

# vis2ualize data distribution in 1- and 2D
rid <- c(
  'markerID',
  'model'
  # 'TIQR',
  # 'PMedian',
  # 'PIQR',
  # 'dPDownSum',
  # 'dTDownSum',
  # 'dXDownSum',
  # 'dYDownSum',
  # 'deltaTPGrads',
  # 'TMedian',
  # 'maxConsecutiveTDown',
  # 'dXSum',
  # 'TDown',
  # 'TUp',
  # 'PUp',
  # 'maxConsecutivePUp',
  # 'maxConsecutiveTUp',
  # 'dXUpSum',
  # 'minTPGrad',
  # 'maxTPGrad',
  # 'maxConsecutivePDown',
  # 'PDown',
  # 'TMax',
  # 'PMax'
)
features <-
  colnames(marx$mark.ft[[1]])[!(colnames(marx$mark.ft[[1]]) %in% rid)]

# Build models by adding layers
# l.kmean ........... k-means (assume k = # of clusters)
# l.bic ............. Bayesian Information Criterion
# l.gmm ............. Gaussian Mixture Model (assume G = # of components)
# l.gmm.da .......... GMM discriminant analysis (supervis2ed classification)
# l.gmm.dr .......... GMM dimension reduction
# l.den ............. GMM layer suited for density plots
# Add layers of vis2ualizations
# l.oned ............ vis2ualize 1D data
# l.twod ............ vis2ualize 2D data
# l.gif ............. creates gifs of marker position ('xy') and PT
#
# Summarise by faceting plots
# f.summary ........ compiles selected gifs into one

m <- marx %>%
  mutate(
    # 1D vis2 plots
    vis2.oned = mark.ft %>%
      l.oned(
        features = features,
        plot = 'all',
        xlim = c(-2, 3),
        dlim = 0.3,
        alpha.min = 0
      ),
    # BIC
    bic.mod0 = mark.ft %>%
      l.bic(
        features = features,
        init = 'SVD',
        scale = TRUE,
        plot = FALSE
      ),
    # General GMM model
    mod0 = l.gmm(
      mark.ft,
      bic.mod0,
      features = features,
      scale = TRUE
    ),
    # Dimension Reduction
    dr.mod0 = l.gmm.dr(mod0),
    # BIC
    bic.mod1 = mark.ft %>%
      l.bic(
        G = 3,
        init = 'SVD',
        features = features,
        scale = TRUE,
        plot = FALSE
      ),
    # General GMM model
    mod1 = l.gmm(
      mark.ft,
      bic.mod1,
      features = features,
      scale = TRUE
    ),
    # Dimension Reduction
    dr.mod1 = l.gmm.dr(mod1)
  )
# # Explore features for individual models
print(m$vis2.oned[[1]])
# Summarise features by faceting
marx$mark.ft %>%
  f.oned(
    nrow = 4,
    ncol = 4,
    features = features,
    xlim = c(-2, 3),
    dlim = 0.3,
    alpha.min = 0.02
  ) %>%
  print()
# vis2ualize model parameters and results
# Plot BIC
m$bic.mod0 %>% p.BIC()
m$bic.mod1 %>% p.BIC()
# # Plot scatter
# m$mod0 %>% p.class()
# m$mod1 %>% p.class()
# 1D facet plots
f.oned(
  m$mark.ft,
  m$mod0,
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
  m$mark.ft,
  m$dr.mod0,
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
  m$mark.ft,
  m$mod1,
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
  m$mark.ft,
  m$dr.mod1,
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
# Close PDF
dev.off()
# vis2ualise fit
f.summary(
  m$data.tidy,
  m$mod0,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis2.mod0'
) %>%
  print()
f.summary(
  m$data.tidy,
  m$dr.mod0,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis2.dr.mod0'
) %>%
  print()
f.summary(
  m$data.tidy,
  m$mod1,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis2.mod1'
) %>%
  print()
f.summary(
  m$data.tidy,
  m$dr.mod1,
  runs = 'cde46',
  GIF = 'PT',
  ncol = 4,
  nrow = 4,
  grads = TRUE,
  save = TRUE,
  fname = 'vis2.dr.mod1'
) %>%
  print()