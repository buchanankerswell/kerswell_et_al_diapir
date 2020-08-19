source('functions.R')

# Open PDF for plotting
cairo_pdf('vis21.PDF',
          width = 11,
          height = 8.5,
          onefile = TRUE)

# Load data and glimpse its structure
# load('marx.RData')

marx
glimpse(marx)

# Visualize data distribution in 1- and 2D
rid <- c(
  'markerID',
  'model',
  'TIQR',
  'PMedian',
  'PIQR',
  'dPDownSum',
  'dTDownSum',
  'dXDownSum',
  'dYDownSum',
  'deltaTPGrads',
  'TMedian',
  'maxConsecutiveTDown',
  'dXSum',
  'TDown',
  'TUp',
  'PUp',
  'maxConsecutivePUp',
  'maxConsecutiveTUp',
  'dXUpSum',
  'minTPGrad',
  'maxTPGrad',
  'maxConsecutivePDown',
  'PDown',
  'TMax',
  'PMax'
)
features <-
  colnames(marx$mark.ft[[1]])[!(colnames(marx$mark.ft[[1]]) %in% rid)]

# Build models by adding layers
# l.kmean ........... k-means (assume k = # of clusters)
# l.bic ............. Bayesian Information Criterion
# l.gmm ............. Gaussian Mixture Model (assume G = # of components)
# l.gmm.da .......... GMM discriminant analysis (supervised classification)
# l.gmm.dr .......... GMM dimension reduction
# l.den ............. GMM layer suited for density plots
# Add layers of visualizations
# l.oned ............ Visualize 1D data
# l.twod ............ Visualize 2D data
# l.gif ............. creates gifs of marker position ('xy') and PT
#
# Summarise by faceting plots
# f.summary ........ compiles selected gifs into one

m <- marx %>%
  mutate(
    # 1D vis plots
    vis.oned = mark.ft %>%
      l.oned(
        features = features,
        plot = 'all',
        xlim = c(-2, 3),
        dlim = 0.3,
        alpha.min = 0,
        bw = 0.5
      ),
    # BIC
    bic.mod0 = mark.ft %>%
      l.bic(
        features = features,
        scale = TRUE,
        plot = FALSE
      ),
    # General GMM model
    mod0 = l.gmm(
      lst = mark.ft,
      lst.bic = bic.mod0,
      features = features,
      G = NULL,
      scale = TRUE,
      hists = FALSE,
      scatter = FALSE
    ),
    # Dimension Reduction
    dr.mod0 = l.gmm.dr(mod0),
    # BIC
    bic.mod1 = mark.ft %>%
      l.bic(
        G = 3,
        features = features,
        scale = TRUE,
        plot = FALSE
      ),
    # General GMM model
    mod1 = l.gmm(
      lst = mark.ft,
      lst.bic = bic.mod1,
      features = features,
      G = NULL,
      scale = TRUE,
      hists = FALSE,
      scatter = FALSE
    ),
    # Dimension Reduction
    dr.mod1 = l.gmm.dr(mod1)
  )
# # Explore features for individual models
print(m$vis.oned[[1]])
# Summarise features by faceting
marx$mark.ft %>%
  f.oned(
    runs = 'all',
    nrow = 2,
    ncol = 2,
    features = features,
    plot = 'all',
    xlim = c(-2, 3),
    dlim = 0.3,
    alpha.min = 0.05,
    bw = 0.5
  ) %>%
  print()
# Visualize model parameters and results
# Plot BIC
m$bic.mod0 %>% p.BIC()
m$bic.mod1 %>% p.BIC()
# # Plot scatter
# m$mod0 %>% p.class()
# m$mod1 %>% p.class()
# Dimension reduction plots
# m$dr.mod0 %>% p.dr.boundary()
# m$dr.mod1 %>% p.dr.boundary()
# 1D Plots
m$mark.ft %>%
  f.oned(
    mods = m$mod0,
    mod.type = NULL,
    dr = FALSE,
    runs = 'all',
    nrow = 2,
    ncol = 2,
    features = features,
    plot = 'ridge',
    xlim = c(-2, 3),
    dlim = 0.3,
    alpha.min = 0,
    bw = 0.5
  ) %>%
  print()
m$mark.ft %>%
  f.oned(
    mods = m$dr.mod0,
    mod.type = 'dr',
    dr = TRUE,
    runs = 'all',
    nrow = 2,
    ncol = 2,
    features = features,
    plot = 'ridge',
    xlim = NULL,
    dlim = 0.3,
    alpha.min = 0,
    bw = 0.5
  ) %>%
  print()
m$mark.ft %>%
  f.oned(
    mods = m$mod1,
    mod.type = NULL,
    dr = FALSE,
    runs = 'all',
    nrow = 2,
    ncol = 2,
    features = features,
    plot = 'ridge',
    xlim = c(-2, 3),
    dlim = 0.3,
    alpha.min = 0,
    bw = 0.5
  ) %>%
  print()
m$mark.ft %>%
  f.oned(
    mods = m$dr.mod1,
    mod.type = 'dr',
    dr = TRUE,
    runs = 'all',
    nrow = 2,
    ncol = 2,
    features = features,
    plot = 'ridge',
    xlim = NULL,
    dlim = 0.3,
    alpha.min = 0,
    bw = 0.5
  ) %>%
  print()
# Turn off PDF
dev.off()
# # Summarise features by animation
# m$mark.ft %>%
#   a.oned(
#     mods = m$mod1,
#     plot = 'all',
#     features = features,
#     type = 'gif',
#     xlim = c(-2, 3),
#     dlim = 0.3,
#     alpha.min = 0.05,
#     bw = 0.5,
#     save = TRUE,
#     fname = 'oned',
#     fps = 20
#   )
# m$mark.ft %>%
#   a.twod(
#     type = 'gif',
#     features = features,
#     xlim = c(-2, 3),
#     dlim = 0.05,
#     save = TRUE,
#     fname = 'twod',
#     fps = 20
#   )