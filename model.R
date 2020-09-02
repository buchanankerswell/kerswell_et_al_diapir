source('functions.R')

# Load data and glimpse its structure
# load('marx.RData')

print(marx)
glimpse(marx)

# Select features
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

# visualize data distribution in 1- and 2D

# Build models by adding layers
# l.kmean ........... k-means (assume k = # of clusters)
# l.bic ............. Bayesian Information Criterion
# l.gmm ............. Gaussian Mixture Model (assume G = # of components)
# l.gmm.da .......... GMM discriminant analysis (supervised classification)
# l.gmm.dr .......... GMM dimension reduction
# l.den ............. GMM layer suited for density plots
# Add layers of visualizations
# l.oned ............ visualize 1D data
# l.twod ............ visualize 2D data
# l.gif ............. creates gifs of marker position ('xy') and PT
#
# Summarise by faceting plots
# f.gif ........ compiles selected gifs into one

mods <- marx %>%
  mutate(
    # 1D vis plots
    vis.oned = mark.ft %>%
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
        init = NULL,
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
        init = NULL,
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
    dr.mod1 = l.gmm.dr(mod1),
    # BIC
    bic.mod2 = mark.ft %>%
      l.bic(
        G = 3,
        init = NULL,
        features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
        scale = TRUE,
        plot = FALSE
      ),
    # General GMM model
    mod2 = l.gmm(
      mark.ft,
      bic.mod2,
      features = c('tsteps', 'tMax', 'dPSum', 'dPUpSum', 'dTSum', 'dTUpSum'),
      scale = TRUE
    ),
    # Dimension Reduction
    dr.mod2 = l.gmm.dr(mod2),
    # BIC
    bic.mod3 = mark.ft %>%
      l.bic(
        G = 3,
        init = NULL,
        features = c('tsteps', 'tMax'),
        scale = TRUE,
        plot = FALSE
      ),
    # General GMM model
    mod3 = l.gmm(
      mark.ft,
      bic.mod3,
      features = c('tsteps', 'tMax'),
      scale = TRUE
    )
  )

# Clean up environment
rm(rid, marx)
rm(list=lsf.str())
save.image('mods.RData')
