source('functions.R')

# Load data and glimpse its structure
# load('marx.RData')

marx
glimpse(marx)

# Visualize data distribution in 1- and 2D
rid <- c('TIQR',
         +          'PMedian',
         +          'PIQR',
         +          'dPDownSum',
         +          'dTDownSum',
         +          'dXDownSum',
         +          'dYDownSum',
         +          'deltaTPGrads',
         +          'TMedian',
         +          'maxConsecutiveTDown',
         +          'dXSum',
         +          'TDown')
features <- colnames(marx$mark.ft[[1]])[colnames(marx$mark.ft[[1]]) %in% rid]
vis9 <- marx %>%
  mutate(
    # 1D vis plots
    vis.oned = mark.ft %>%
      l.oned(
        features = 'all',
        plot = 'all',
        xlim = c(-3, 3),
        dlim = 0.5
      ),
    # 2D vis plots
    vis.twod = mark.ft %>%
      l.twod(
        features = 'all',
        xlim = c(-3, 3),
        dlim = 0.05
      )
  )

# Open PDF for plotting
cairo_pdf('vis9.PDF',
          width = 11,
          height = 8.5,
          onefile = TRUE)
# # Explore features for individual models
# print(vis9$vis.oned[[1]])
# print(vis9$vis.twod[[1]][[1]])
# Summarise features by faceting
vis9$mark.ft %>%
  f.oned(
    runs = 'all',
    nrow = 4,
    ncol = 4,
    features = 'all',
    plot = 'all',
    xlim = c(-3, 3),
    dlim = 0.5
  ) %>%
  print()

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

# mods1 <- marx %>%
#   mutate(
    # # BIC
    # bic.mod0 = mark.ft %>%
    #   l.bic(
    #     features = c('dTSum', 'dPSum', 'dTUpSum', 'dPUpSum'),
    #     scale = TRUE,
    #     plot = FALSE
    #   ),
    # # General GMM model
    # mod0 = l.gmm(
    #   lst = mark.ft,
    #   lst.bic = bic.mod0,
    #   features = c('dTSum', 'dPSum', 'dTUpSum', 'dPUpSum'),
    #   G = NULL,
    #   scale = TRUE,
    #   hists = FALSE,
    #   scatter = FALSE
    # ),
    # # Dimension Reduction
    # dr.mod0 = l.gmm.dr(mod0),
    # # BIC
    # bic.mod1 = mark.ft %>%
    #   l.bic(
    #     G = 3,
    #     features = c('dTSum', 'dPSum', 'dTUpSum', 'dPUpSum'),
    #     scale = TRUE,
    #     plot = FALSE
    #   ),
    # # General GMM model
    # mod1 = l.gmm(
    #   lst = mark.ft,
    #   lst.bic = bic.mod1,
    #   features = c('dTSum', 'dPSum', 'dTUpSum', 'dPUpSum'),
    #   G = NULL,
    #   scale = TRUE,
    #   hists = FALSE,
    #   scatter = FALSE
    # ),
    # # Dimension Reduction
    # dr.mod1 = l.gmm.dr(mod1)
  # )
# # Visualize model parameters and results
# # Plot BIC
# mods1$bic.mod0 %>% p.BIC()
# mods1$bic.mod1 %>% p.BIC()
# # Plot scatter
# mods1$mod0 %>% p.class()
# mods1$mod1 %>% p.class()
# # Dimension reduction plots
# mods1$dr.mod0 %>% p.dr.boundary()
# mods1$dr.mod1 %>% p.dr.boundary()
# Turn off PDF
dev.off()
# # Summarise features by animation
# mods1$mark.ft %>%
#   a.oned(
#     plot = 'all',
#     features = c('dTSum', 'dPSum', 'dTUpSum', 'dPUpSum'),
#     type = 'gif',
#     xlim = c(-2, 2),
#     dlim = 1,
#     save = TRUE,
#     fname = 'oned',
#     fps = 20
#   )
# mods1$mark.ft %>%
#   a.twod(
#     type = 'gif',
#     features = c('dTSum', 'dPSum', 'dTUpSum', 'dPUpSum'),
#     xlim = c(-2, 2),
#     dlim = 0.05,
#     save = TRUE,
#     fname = 'twod',
#     fps = 20
#   )