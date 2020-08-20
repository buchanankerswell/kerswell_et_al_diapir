if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyverse,
  vroom,
  mclust,
  fs,
  ggforce,
  patchwork,
  gridExtra,
  gganimate,
  ggridges,
  transformr,
  gifski
)
# mclust color palette
c.pal <- RColorBrewer::brewer.pal(n = 12, name = 'Paired')
bic.col <- mclust.options("bicPlotColors")
bic.col[1:14] <- c(c.pal, c.pal[1:2])
mclust.options("bicPlotColors" = bic.col)
mclust.options("classPlotColors" = c.pal)
rm(c.pal, bic.col)
# proportions ----
# Calculate ratios of rocks among classes
proportions <- function(df) {
  if (any(colnames(df) %in% 'cluster')) {
    clusterProps <- df %>%
      group_by(cluster) %>%
      count(markerID) %>%
      summarise(markerCount = n_distinct(markerID), .groups = 'drop') %>%
      mutate(proportion = markerCount / sum(markerCount))
  } else {
    clusterProps <- df %>%
      group_by(class) %>%
      count(markerID) %>%
      summarise(markerCount = n_distinct(markerID), .groups = 'drop') %>%
      mutate(proportion = markerCount / sum(markerCount))
  }
  
  return(clusterProps)
}
# summaryStats ----
# Calculate summary statistics
summaryStats <- function(dflist) {
  # Wrangle data
  subductedMarkers <-
    lapply(dflist, . %>% filter(fate == 'subducted')) %>%
    bind_rows()
  escalatorMarkers <-
    lapply(dflist, . %>% filter(fate == 'escalator')) %>%
    bind_rows()
  elevatorMarkers <-
    lapply(dflist, . %>% filter(fate == 'elevator')) %>%
    bind_rows()
  # n
  nSubducted <- subductedMarkers %>%
    group_by(model) %>%
    summarise(nSubducted = n_distinct(markerID))
  nEscalator <- escalatorMarkers %>%
    group_by(model) %>%
    summarise(nEscalator = n_distinct(markerID))
  nElevator <- elevatorMarkers %>%
    group_by(model) %>%
    summarise(nElevator = n_distinct(markerID))
  n <- list(nSubducted, nEscalator, nElevator) %>%
    reduce(left_join, by = 'model') %>%
    mutate(nTotal = nSubducted + nEscalator + nElevator)
  # Maximum pressures
  maxPEscalator <- escalatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(maxPEscalator = max(Pkbar, na.rm = TRUE))
  maxPElevator <- elevatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(maxPElevator = max(Pkbar, na.rm = TRUE))
  PSummary <-
    list(maxPEscalator,
         maxPElevator) %>%
    reduce(left_join) %>%
    ungroup()
  # Maximum temperatures
  maxTEscalator <- escalatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(maxTEscalator = max(TCel, na.rm = TRUE))
  maxTElevator <- elevatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(maxTElevator = max(TCel, na.rm = TRUE))
  TSummary <-
    list(maxTEscalator,
         maxTElevator) %>%
    reduce(left_join) %>%
    ungroup()
  # Minimum and maximum T/P and T/Depth gradients
  maxTPGradSubducted <- subductedMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    transmute(maxTPGradSubducted = TCel / Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      maxTPGradSubductedMax = max(maxTPGradSubducted),
      maxTPGradSubductedMin = min(maxTPGradSubducted),
      maxTPGradSubductedMean = mean(maxTPGradSubducted),
      maxTPGradSubductedMedian = median(maxTPGradSubducted)
    )
  maxTDepthGradSubducted <- subductedMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    transmute(maxTDepthGradSubducted = TCel / yposkm) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      maxTDepthGradSubductedMax = max(maxTDepthGradSubducted),
      maxTDepthGradSubductedMin = min(maxTDepthGradSubducted),
      maxTDepthGradSubductedMean = mean(maxTDepthGradSubducted),
      maxTDepthGradSubductedMedian = median(maxTDepthGradSubducted)
    )
  minTPGradSubducted <- subductedMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    transmute(minTPGradSubducted = TCel / Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      minTPGradSubductedMax = max(minTPGradSubducted),
      minTPGradSubductedMin = min(minTPGradSubducted),
      minTPGradSubductedMean = mean(minTPGradSubducted),
      minTPGradSubductedMedian = median(minTPGradSubducted)
    )
  minTDepthGradSubducted <- subductedMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    transmute(minTDepthGradSubducted = TCel / yposkm) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      minTDepthGradSubductedMax = max(minTDepthGradSubducted),
      minTDepthGradSubductedMin = min(minTDepthGradSubducted),
      minTDepthGradSubductedMean = mean(minTDepthGradSubducted),
      minTDepthGradSubductedMedian = median(minTDepthGradSubducted)
    )
  maxTPGradEscalator <- escalatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    transmute(maxTPGradEscalator = TCel / Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      maxTPGradEscalatorMax = max(maxTPGradEscalator),
      maxTPGradEscalatorMin = min(maxTPGradEscalator),
      maxTPGradEscalatorMean = mean(maxTPGradEscalator),
      maxTPGradEscalatorMedian = median(maxTPGradEscalator)
    )
  maxTDepthGradEscalator <- escalatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    transmute(maxTDepthGradEscalator = TCel / yposkm) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      maxTDepthGradEscalatorMax = max(maxTDepthGradEscalator),
      maxTDepthGradEscalatorMin = min(maxTDepthGradEscalator),
      maxTDepthGradEscalatorMean = mean(maxTDepthGradEscalator),
      maxTDepthGradEscalatorMedian = median(maxTDepthGradEscalator)
    )
  minTPGradEscalator <- escalatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    transmute(minTPGradEscalator = TCel / Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      minTPGradEscalatorMax = max(minTPGradEscalator),
      minTPGradEscalatorMin = min(minTPGradEscalator),
      minTPGradEscalatorMean = mean(minTPGradEscalator),
      minTPGradEscalatorMedian = median(minTPGradEscalator)
    )
  minTDepthGradEscalator <- escalatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    transmute(minTDepthGradEscalator = TCel / yposkm) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      minTDepthGradEscalatorMax = max(minTDepthGradEscalator),
      minTDepthGradEscalatorMin = min(minTDepthGradEscalator),
      minTDepthGradEscalatorMean = mean(minTDepthGradEscalator),
      minTDepthGradEscalatorMedian = median(minTDepthGradEscalator)
    )
  maxTPGradElevator <- elevatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    transmute(maxTPGradElevator = TCel / Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      maxTPGradElevatorMax = max(maxTPGradElevator),
      maxTPGradElevatorMin = min(maxTPGradElevator),
      maxTPGradElevatorMean = mean(maxTPGradElevator),
      maxTPGradElevatorMedian = median(maxTPGradElevator)
    )
  maxTDepthGradElevator <- elevatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(TCel) %>%
    transmute(maxTDepthGradElevator = TCel / yposkm) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      maxTDepthGradElevatorMax = max(maxTDepthGradElevator),
      maxTDepthGradElevatorMin = min(maxTDepthGradElevator),
      maxTDepthGradElevatorMean = mean(maxTDepthGradElevator),
      maxTDepthGradElevatorMedian = median(maxTDepthGradElevator)
    )
  minTPGradElevator <- elevatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    transmute(minTPGradElevator = TCel / Pkbar) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      minTPGradElevatorMax = max(minTPGradElevator),
      minTPGradElevatorMin = min(minTPGradElevator),
      minTPGradElevatorMean = mean(minTPGradElevator),
      minTPGradElevatorMedian = median(minTPGradElevator)
    )
  minTDepthGradElevator <- elevatorMarkers %>%
    group_by(model, markerID) %>%
    slice_max(Pkbar) %>%
    transmute(minTDepthGradElevator = TCel / yposkm) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(
      minTDepthGradElevatorMax = max(minTDepthGradElevator),
      minTDepthGradElevatorMin = min(minTDepthGradElevator),
      minTDepthGradElevatorMean = mean(minTDepthGradElevator),
      minTDepthGradElevatorMedian = median(minTDepthGradElevator)
    )
  GradSummary <-
    list(
      maxTPGradSubducted,
      maxTDepthGradSubducted,
      minTPGradSubducted,
      minTDepthGradSubducted,
      maxTPGradEscalator,
      maxTDepthGradEscalator,
      minTPGradEscalator,
      minTDepthGradEscalator,
      maxTPGradElevator,
      maxTDepthGradElevator,
      minTPGradElevator,
      minTDepthGradElevator
    ) %>%
    reduce(left_join) %>%
    ungroup()
  # Summary table
  dataSummary <- n %>%
    left_join(PSummary) %>%
    left_join(TSummary) %>%
    left_join(GradSummary)
  return(dataSummary)
}
# data.get ----
# Read and process single data file
data.get <- function(path) {
  data <-
    # Read tab-delimited data file
    vroom(path,
          col_names = FALSE,
          delim = '\t',
          id = "path")
  # Get model names
  fname <- substring(data[1, 1],
                     regexpr("/cd", data[1, 1]) + 1,
                     regexpr("_marks", data[1, 1]) - 1)
  # Number of markers
  nummarkers <- nrow(data) / 6
  # Count timesteps
  ntstep <- ncol(data) - 1
  # Select number of timesteps (columns)
  data <- data[, 2:as.numeric(ntstep + 1)]
  # Add column for marker type, time, position, temp, pressure, and markerID
  data <- data %>%
    mutate(
      marker = rep(c('type', 'time', 'xpos', 'ypos', 'T', 'P'), each = nummarkers),
      markerID = as.factor(rep(c(1:nummarkers), 6)),
      model = as.factor(rep(fname, nrow(data)))
    )
  # Get rid of rows where ALL marker types are NA by removing rows whose sum == 0
  # and filter the rest of the data to select only the markers that are not NAs
  data <- data %>%
    filter((rowSums(.[, seq_along(ntstep)], na.rm = TRUE) > 0) == TRUE)
  assign(fname, data)
}
# data.tidy ----
data.tidy <- function(df) {
  dat <- df %>%
    select(where(is.numeric), marker, markerID, model) %>%
    pivot_longer(-c(markerID, marker, model),
                 names_to = 'timestep',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'marker', values_from = 'value') %>%
    mutate(
      z1100 = as.numeric(substr(df$model[1], 4, nchar(
        as.character(df$model[1])
      ))),
      Pkbar = P / 1000,
      TCel = T - 273.15,
      xposkm = (xpos / 1000) - 20,
      yposkm = (ypos / 1000) - 20
    ) %>%
    drop_na()
  return(dat)
}
# add.ft ----
# Adds features (columns) for ML clustering
add.ft <- function(df) {
  removeCol <- colnames(df)
  removeCol <- removeCol[-c(which(removeCol == 'markerID'),
                            which(removeCol == 'model'))]
  dat <- df %>%
    group_by(markerID) %>%
    mutate(
      PMax = max(Pkbar, na.rm = TRUE),
      PMedian = median(Pkbar, na.rm = TRUE),
      PIQR = IQR(Pkbar, na.rm = TRUE),
      PUp = sum(diff(as.numeric(Pkbar)) > 0, na.rm = TRUE),
      PDown = sum(diff(as.numeric(Pkbar)) < 0, na.rm = TRUE),
      maxConsecutivePUp =
        if (sum(diff(as.numeric(Pkbar)) > 0, na.rm = TRUE) > 0) {
          runs <- rle(diff(as.numeric(Pkbar)) > 0)
          maxRun <-
            max(runs$lengths[runs$values == TRUE], na.rm = TRUE)
        } else{
          maxRun <- 0
        },
      maxConsecutivePDown =
        if (sum(diff(as.numeric(Pkbar)) < 0, na.rm = TRUE) > 0) {
          runs <- rle(diff(as.numeric(Pkbar)) < 0)
          maxRun <-
            max(runs$lengths[runs$values == TRUE], na.rm = TRUE)
        } else{
          maxRun <- 0
        },
      dPSum = sum(diff(as.numeric(Pkbar)), na.rm = TRUE),
      dPUpSum = sum(diff(as.numeric(Pkbar))[which(diff(as.numeric(Pkbar)) > 0)]),
      dPDownSum = sum(diff(as.numeric(Pkbar))[which(diff(as.numeric(Pkbar)) < 0)]),
      TMax = max(TCel, na.rm = TRUE),
      TMedian = median(TCel, na.rm = TRUE),
      TIQR = IQR(TCel, na.rm = TRUE),
      TUp = sum(diff(as.numeric(TCel)) > 0, na.rm = TRUE),
      maxConsecutiveTUp =
        if (sum(diff(as.numeric(TCel)) > 0, na.rm = TRUE) > 0) {
          runs <- rle(diff(as.numeric(TCel)) > 0)
          maxRun <-
            max(runs$lengths[runs$values == TRUE], na.rm = TRUE)
        } else{
          maxRun <- 0
        },
      TDown = sum(diff(as.numeric(TCel)) < 0, na.rm = TRUE),
      maxConsecutiveTDown =
        if (sum(diff(as.numeric(TCel)) < 0, na.rm = TRUE) > 0) {
          runs <- rle(diff(as.numeric(TCel)) < 0)
          maxRun <-
            max(runs$lengths[runs$values == TRUE], na.rm = TRUE)
        } else{
          maxRun <- 0
        },
      dTSum = sum(diff(as.numeric(TCel)), na.rm = TRUE),
      dTUpSum = sum(diff(as.numeric(TCel))[which(diff(as.numeric(TCel)) > 0)]),
      dTDownSum = sum(diff(as.numeric(TCel))[which(diff(as.numeric(TCel)) < 0)]),
      dXSum = sum(diff(as.numeric(xposkm)), na.rm = TRUE),
      dXUpSum = sum(diff(as.numeric(xposkm))[which(diff(as.numeric(xposkm)) > 0)]),
      dXDownSum = sum(diff(as.numeric(xposkm))[which(diff(as.numeric(xposkm)) < 0)]),
      dYSum = sum(diff(as.numeric(yposkm)), na.rm = TRUE),
      dYUpSum = sum(diff(as.numeric(yposkm))[which(diff(as.numeric(yposkm)) > 0)]),
      dYDownSum = sum(diff(as.numeric(yposkm))[which(diff(as.numeric(yposkm)) < 0)])
    )
  # More features
  maxTPGrad <- dat %>%
    slice_max(TCel) %>%
    transmute(maxTPGrad = as.numeric(TCel / Pkbar))
  minTPGrad <- dat %>%
    slice_max(Pkbar) %>%
    transmute(minTPGrad = as.numeric(TCel / Pkbar))
  # Join data
  dat <- list(dat, maxTPGrad, minTPGrad) %>%
    reduce(left_join) %>%
    ungroup() %>%
    mutate(deltaTPGrads = as.numeric(maxTPGrad - minTPGrad)) %>%
    drop_na()
  dat <- dat %>%
    mutate(across(where(is.integer), as.numeric))
  dat <- dat %>%
    mutate(across(where(is.character), as.factor))
  dat <- dat %>%
    mutate(across(contains('timestep'), as.numeric)) %>%
    arrange(markerID, timestep) %>%
    group_by(markerID) %>%
    slice_head(1) %>%
    select(-all_of(removeCol))
  return(dat)
}

# add.couplingft ----
# Add kerswell_et_al_coupling features
add.couplingft <- function(df) {
  dat <- df %>%
    left_join(dataCoupling)
}
# k.mean ----
# k-means clustering
k.mean <- function(df,
                   k,
                   nstart,
                   features,
                   scale = FALSE,
                   scree = FALSE,
                   hists = FALSE) {
  if (features == 'all') {
    features <- df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames()
  } else {
    features <- features
  }
  kmData <- df %>%
    ungroup() %>%
    select(which(colnames(df) %in% features))
  if (scale == FALSE) {
    m <- kmeans(kmData, k, nstart)
  } else {
    m <- kmeans(scale(kmData), k, nstart)
  }
  df$cluster <- as.factor(m$cluster)
  if (scree == TRUE) {
    numCenters <- 1:10
    withinSS <- sapply(numCenters, function(k) {
      m <- kmeans(kmData, k, nstart = nstart)
      return(m$tot.withinss)
    })
    plot(
      numCenters,
      withinSS,
      type = 'b',
      pch = 20,
      ylab = 'Within Sum of Squares',
      xlab = 'Number of Centroids'
    )
  }
  if (hists == TRUE) {
    p <-
      lapply(which(colnames(df) %in% colnames(kmData)), function(i) {
        ggplot(df) +
          geom_histogram(aes_string(colnames(df)[i], fill = 'cluster'), color = 'white') +
          theme_bw(base_size = 14) +
          ggtitle(paste0(df$model[1], ' | ', colnames(df)[i])) +
          scale_fill_brewer(palette = 'Pairs') +
          theme(
            axis.ticks.length = unit(-0.25, "cm"),
            axis.text.x = element_text(margin = unit(c(0.5, 0, 0, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_rect(size = 1, color = 'black'),
            axis.text = element_text(face = 'plain', color = 'black'),
            axis.ticks = element_line(size = 0.5, color = 'black'),
            legend.direction = 'horizontal',
            legend.justification = c(0, 1),
            legend.position = 'none',
            legend.box.margin = margin(c(50, 50, 50, 50)),
            axis.title = element_text(size = 12, face = 'plain'),
            plot.title = element_text(size = 14, face = 'plain')
          )
      })
    grid.arrange(grobs = p[1:length(features)],
                 nrow = ceiling(sqrt(length(features))),
                 ncol = ceiling(sqrt(length(features))))
  }
  return(m)
}
# b.ic ----
# Bayesian Information Criterion
b.ic <- function(df,
                 features = 'all',
                 G = NULL,
                 scale = FALSE,
                 plot = TRUE) {
  if (features == 'all') {
    features <- df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames()
  } else {
    features <- features
  }
  gmmData <- df %>%
    ungroup() %>%
    select(which(colnames(df) %in% features))
  if (scale == FALSE) {
    BIC <- mclustBIC(gmmData, G = G)
  } else {
    BIC <- mclustBIC(scale(gmmData), G = G)
  }
  print(summary(BIC))
  if (plot == TRUE) {
    plot(BIC)
  }
  return(BIC)
}
# gmm ----
# Gaussian mixing model (GMM) classification
gmm <- function(df,
                x = NULL,
                features = 'all',
                G = NULL,
                scale = FALSE,
                hists = FALSE,
                scatter = FALSE) {
  if (features == 'all') {
    features <- df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames()
  } else {
    features <- features
  }
  gmmData <- df %>%
    ungroup() %>%
    select(which(colnames(df) %in% features))
  if (scale == FALSE) {
    m <- Mclust(gmmData, G = G, x = x)
  } else {
    m <- Mclust(scale(gmmData), G = G, x = x)
  }
  print(summary(m))
  df$class <- as.factor(m$classification)
  if (hists == TRUE) {
    p <-
      lapply(which(colnames(df) %in% colnames(gmmData)), function(i) {
        ggplot(df) +
          geom_histogram(aes_string(colnames(df)[i], fill = 'class'), color = 'white') +
          theme_bw(base_size = 14) +
          ggtitle(paste0(df$model[1], ' | ', colnames(df)[i])) +
          scale_fill_brewer(palette = 'Pairs') +
          theme(
            axis.ticks.length = unit(-0.25, "cm"),
            axis.text.x = element_text(margin = unit(c(0.5, 0, 0, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_rect(size = 1, color = 'black'),
            axis.text = element_text(face = 'plain', color = 'black'),
            axis.ticks = element_line(size = 0.5, color = 'black'),
            legend.direction = 'horizontal',
            legend.justification = c(0, 1),
            legend.position = 'none',
            legend.box.margin = margin(c(50, 50, 50, 50)),
            axis.title = element_text(size = 12, face = 'plain'),
            plot.title = element_text(size = 14, face = 'plain')
          )
      })
    grid.arrange(grobs = p[1:length(features)],
                 nrow = ceiling(sqrt(length(features))),
                 ncol = ceiling(sqrt(length(features))))
  }
  if (scatter == TRUE) {
    plot(m, what = "classification")
  }
  return(m)
}
# gmm.da ----
# GMM Discriminant analysis (supervised classification)
gmm.da <- function(df,
                   features = 'all',
                   classCol = 'class',
                   scale = FALSE,
                   scatter = FALSE,
                   hists = TRUE) {
  if (features == 'all') {
    features <- df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames()
  } else {
    features <- features
  }
  cls <- df %>%
    ungroup() %>%
    select(all_of(classCol)) %>%
    mutate(across(where(is.factor), as.character))
  cls <- cls[, 1, drop = TRUE]
  gmmData <- df %>%
    ungroup() %>%
    select(which(colnames(df) %in% features))
  if (scale == FALSE) {
    m <- MclustDA(gmmData, cls, modelType = 'EDDA')
  } else {
    m <- MclustDA(scale(gmmData), cls, modelType = 'EDDA')
  }
  print(summary(m))
  df$class <- m$class
  if (scatter == TRUE) {
    plot(m, what = "scatter")
  }
  if (hists == TRUE) {
    p <-
      lapply(which(colnames(df) %in% colnames(gmmData)), function(i) {
        ggplot(df) +
          geom_histogram(aes_string(colnames(df)[i], fill = 'class'), color = 'white') +
          theme_bw(base_size = 14) +
          ggtitle(paste0(df$model[1], ' | ', colnames(df)[i])) +
          scale_fill_brewer(palette = 'Pairs') +
          theme(
            axis.ticks.length = unit(-0.25, "cm"),
            axis.text.x = element_text(margin = unit(c(0.5, 0, 0, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_rect(size = 1, color = 'black'),
            axis.text = element_text(face = 'plain', color = 'black'),
            axis.ticks = element_line(size = 0.5, color = 'black'),
            legend.direction = 'horizontal',
            legend.justification = c(0, 1),
            legend.position = 'none',
            legend.box.margin = margin(c(50, 50, 50, 50)),
            axis.title = element_text(size = 12, face = 'plain'),
            plot.title = element_text(size = 14, face = 'plain')
          )
      })
    grid.arrange(grobs = p[1:length(features)],
                 nrow = ceiling(sqrt(length(features))),
                 ncol = ceiling(sqrt(length(features))))
  }
  return(m)
}
# gmm.dr ----
# GMM dimension reduction
gmm.dr <- function(obj.mclust,
                   lambda = 1) {
  m <- MclustDR(obj.mclust, lambda)
  return(m)
}
# gmm.den ----
# GMM density model
gmm.den <- function(df,
                    features,
                    G = NULL,
                    plots = FALSE) {
  if (features == 'all') {
    m <-
      df %>% ungroup() %>% select(where(is.numeric)) %>% densityMclust(G = G)
  } else {
    m <-
      df %>% ungroup() %>% select(all_of(features)) %>% densityMclust(G = G)
  }
  print(summary(m))
  if (plots == 'TRUE') {
    p.persp <- plot(m, what = 'density', type = 'persp')
    p.BIC <- plot(m, what = 'BIC')
    return(list(
      model = m,
      p.persp = p.persp,
      p.BIC = p.BIC
    ))
  } else {
    return(list(model = m))
  }
}
# l.kmean ----
# Add k-means layer
l.kmean <- function(lst, ...) {
  purrr::map(lst, try(k.mean)
             , ...)
}

# l.bic ----
# Add BIC layer
l.bic <- function(lst, ...) {
  purrr::map(lst, try(b.ic)
             , ...)
}
# l.gmm ----
# Add GMM layer
l.gmm <- function(lst, lst.bic = NULL, ...) {
  if (!is.null(lst.bic)) {
    purrr::map2(lst, lst.bic, gmm, ...)
  } else {
    purrr::map(lst, gmm, ...)
  }
}
# l.gmm.da ----
# Add GMMDA layer
l.gmm.da <- function(lst, ...) {
  l <- purrr::map(lst, pluck('data'))
  purrr::map(l, gmm.da, ...)
}
# l.gmm.dr ----
# add GMM DR layer
l.gmm.dr <- function(lst, ...) {
  purrr::map(lst, gmm.dr, ...)
}
# l.density ----
# Add GMM density layer
l.density <- function(lst, ...) {
  purrr::map(lst, try(gmm.density)
             , ...)
}
# l.oned ----
# add list column of 1D density plots
l.oned <- function(lst, ...) {
  purrr::map(lst, p.oned, ...)
}
# l.twod ----
# add list column of 2D density plots
l.twod <- function(lst, ...) {
  purrr::map(lst, p.twod, ...)
}
# l.gif ----
# Adds gif layer to dataframe `data`
l.gif <- function(lst, mods, GIF = c('xy', 'PT')) {
  d <- map2(lst, purrr::map(mods, pluck('data')), left_join)
  if (GIF == 'xy') {
    anim <- d %>% purrr::map(try(gif.xy))
  } else {
    anim <- d %>% purrr::map(try(gif.pt))
  }
}
# p.BIC ----
# BIC Plot
p.BIC <- function(lst, ...) {
  purrr::map(lst, plot)
}
# p.class ----
# Classification scatter plot
p.class <- function(lst, ...) {
  purrr::map(lst, plot, what = 'classification')
}
# p.scatter ----
# Scatter plot
p.scatter <- function(lst, ...) {
  purrr::map(lst, plot, what = 'scatter')
}
# p.persp ----
# Perpective density plot
p.persp <- function(lst, ...) {
  purrr::map(lst, plot, what = 'density', type = 'persp')
}
# p.dr.contour ----
# dimension reduction contour plot
p.dr.contour <- function(lst, ...) {
  purrr::map(lst, plot, what = 'contour')
}
# p.dr.boundary ----
# Dimension reduction boundary plot
p.dr.boundary <- function(lst, ...) {
  purrr::map(lst, plot, what = 'boundaries', ngrid = 200)
}
# p.oned ----
# Visualize features in one-dimension
p.oned <- function(df,
                   mod = NULL,
                   mod.type = NULL,
                   dr = FALSE,
                   sigma = 2,
                   features,
                   plot,
                   xlim = c(-2.5, 2.5),
                   dlim = 0.5,
                   alpha.min = 0.02,
                   bw = 0.5) {
  if (sigma == 2) {
    width <- 'twosigma'
  } else if (sigma == 3) {
    width <- 'threesigma'
  } else {
    width <- 'twosigma'
  }
  fname <- df$model[1]
  if (features == 'all') {
    df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames() ->
      features
  } else {
    features <- features
  }
  df %>%
    ungroup() %>%
    select(all_of(features)) %>%
    pivot_longer(cols = everything(),
                 names_to = 'var',
                 values_to = 'val') %>%
    group_by(var) %>%
    filter(!is.na(val)) %>%
    nest() %>%
    mutate(d.scaled = purrr::map(data, ~ as_tibble(scale(.x)))) %>%
    select(var, d.scaled) %>%
    unnest(cols = d.scaled) %>%
    add_column(model = fname) ->
    df.scaled
  if (!is.null(mod)) {
    if (!is.null(mod.type)) {
      if (mod.type == 'dr' & dr == TRUE) {
        ddr <- mod$dir %>%
          as_tibble() %>%
          pivot_longer(cols = everything(),
                       names_to = 'var',
                       values_to = 'val')
        cntr <- mod$dir %>%
          as_tibble() %>%
          add_column(cls = mod$classification) %>%
          pivot_longer(cols = -cls,
                       names_to = 'var',
                       values_to = 'val') %>%
          mutate(across(where(is.character), factor)) %>%
          group_by(cls, var) %>%
          summarise(
            cntr = mean(val),
            sigma = sd(val),
            twosigma = sd(val) * 2,
            threesigma = sd(val) * 3
          )
        c <- ddr %>% left_join(cntr)
      } else {
        cntr <-
          as_tibble(mod$parameters$mean,
                    rownames = 'var',
                    .name_repair = 'unique') %>%
          rename_with( ~ gsub('...', '', .x), .cols = where(is.numeric)) %>%
          pivot_longer(
            cols = where(is.numeric),
            names_to = 'cls',
            values_to = 'cntr'
          ) %>%
          mutate(across(where(is.character), factor))
        sig <- apply(mod$parameters$variance$sigma, 3, diag) %>%
          as_tibble(rownames = 'var') %>%
          rename_with( ~ gsub('V', '', .x), .cols = where(is.numeric)) %>%
          pivot_longer(
            cols = where(is.numeric),
            names_to = 'cls',
            values_to = 'variance'
          ) %>%
          mutate(
            across(where(is.character), factor),
            twosigma = sqrt(variance) * 2,
            threesigma = sqrt(variance) * 3
          )
        c <- cntr %>% left_join(sig)
      }
    }
  }
  if (plot == 'strip') {
    p.strip <- ggplot(df.scaled) +
      stat_density(
        aes(
          x = val,
          fct_reorder(var, val, median),
          fill = ..density..
        ),
        geom = 'raster',
        position = 'identity',
        adjust = bw
      ) +
      scale_fill_viridis_c(
        option = 'magma',
        limits = c(0, dlim),
        begin = 0,
        end = 1,
        na.value = 'lemonchiffon'
      ) +
      coord_cartesian(xlim = xlim) +
      xlab('Scaled Value') +
      ylab('') +
      ggtitle(fname) +
      theme(
        panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.major = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7)),
        panel.grid.minor = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7))
      )
  } else if (plot == 'ridge') {
    if (!is.null(mod)) {
      if (!is.null(mod.type)) {
        if (mod.type == 'dr' & dr == TRUE) {
          cntr <- c %>% group_by(model, cls) %>% select(-val) %>% distinct()
          ddr <-
            c %>% group_by(var) %>% select(model, var, val) %>% distinct()
          p <- ggplot() +
            geom_tile(
              data = cntr,
              aes_string(
                x = 'cntr',
                y = 'var',
                color = 'cls',
                fill = 'cls',
                width = width,
                group = 'cls'
              ),
              alpha = 0.4,
              height = 0.2
            ) +
            geom_point(
              data = cntr,
              aes(
                x = cntr,
                y = var,
                fill = cls,
                color = cls,
                group = cls
              ),
              shape = 15
            ) +
            geom_density_ridges(
              data = ddr,
              aes(x = val, y = fct_reorder(var, val, median)),
              rel_min_height = alpha.min,
              fill = 'grey50',
              alpha = 0.5
            ) +
            coord_cartesian(xlim = xlim) +
            scale_color_brewer(name = 'Group', palette = 'Paired') +
            scale_fill_brewer(name = 'Group', palette = 'Paired') +
            xlab('Scaled Value') +
            ylab('') +
            ggtitle(fname)
        } else {
          p <- ggplot(df.scaled) +
            geom_tile(
              data = c,
              aes_string(
                x = 'cntr',
                y = 'var',
                color = 'cls',
                fill = 'cls',
                width = width,
                group = 'cls'
              ),
              alpha = 0.4,
              height = 0.2
            ) +
            geom_point(
              data = c,
              aes(
                x = cntr,
                y = var,
                fill = cls,
                color = cls,
                group = cls
              ),
              shape = 15
            ) +
            geom_density_ridges(
              aes(x = val, y = fct_reorder(var, val, median)),
              rel_min_height = alpha.min,
              fill = 'grey50',
              alpha = 0.5
            ) +
            coord_cartesian(xlim = xlim) +
            scale_color_brewer(name = 'Group', palette = 'Paired') +
            scale_fill_brewer(name = 'Group', palette = 'Paired') +
            xlab('Scaled Value') +
            ylab('') +
            ggtitle(fname)
        }
      } else {
        p <- ggplot(df.scaled) +
          geom_tile(
            data = c,
            aes_string(
              x = 'cntr',
              y = 'var',
              color = 'cls',
              fill = 'cls',
              width = width,
              group = 'cls'
            ),
            alpha = 0.4,
            height = 0.2
          ) +
          geom_point(
            data = c,
            aes(
              x = cntr,
              y = var,
              fill = cls,
              color = cls,
              group = cls
            ),
            shape = 15
          ) +
          geom_density_ridges(
            aes(x = val, y = fct_reorder(var, val, median)),
            rel_min_height = alpha.min,
            fill = 'grey50',
            alpha = 0.5
          ) +
          coord_cartesian(xlim = xlim) +
          scale_color_brewer(name = 'Group', palette = 'Paired') +
          scale_fill_brewer(name = 'Group', palette = 'Paired') +
          xlab('Scaled Value') +
          ylab('') +
          ggtitle(fname)
      }
    } else {
      p <- ggplot(df.scaled) +
        geom_density_ridges(
          aes(x = val, fct_reorder(var, val, median)),
          rel_min_height = alpha.min,
          fill = 'deeppink4',
          alpha = 0.5
        ) +
        coord_cartesian(xlim = xlim) +
        xlab('Scaled Value') +
        ylab('') +
        ggtitle(fname)
    }
  } else if (plot == 'all') {
    if (!is.null(mod)) {
      if (!is.null(mod)) {
        if (mod.type == 'dr' & dr == TRUE) {
          cntr <- c %>% group_by(model, cls) %>% select(-val) %>% distinct()
          ddr <-
            c %>% group_by(var) %>% select(model, var, val) %>% distinct()
          p.ridge <- ggplot() +
            geom_tile(
              data = cntr,
              aes_string(
                x = 'cntr',
                y = 'var',
                color = 'cls',
                fill = 'cls',
                width = width,
                group = 'cls'
              ),
              alpha = 0.4,
              height = 0.2
            ) +
            geom_point(
              data = cntr,
              aes(
                x = cntr,
                y = var,
                fill = cls,
                color = cls,
                group = cls
              ),
              shape = 15
            ) +
            geom_density_ridges(
              data = ddr,
              aes(x = val, y = fct_reorder(var, val, median)),
              rel_min_height = alpha.min,
              fill = 'grey50',
              alpha = 0.5
            ) +
            coord_cartesian(xlim = xlim) +
            scale_color_brewer(name = 'Group', palette = 'Paired') +
            scale_fill_brewer(name = 'Group', palette = 'Paired') +
            xlab('Scaled Value') +
            ylab('') +
            ggtitle(fname)
        } else {
          p.ridge <- ggplot(df.scaled) +
            geom_tile(
              data = c,
              aes_string(
                x = 'cntr',
                y = 'var',
                color = 'cls',
                fill = 'cls',
                width = width,
                group = 'cls'
              ),
              alpha = 0.4,
              height = 0.2
            ) +
            geom_point(
              data = c,
              aes(
                x = cntr,
                y = var,
                fill = cls,
                color = cls,
                group = cls
              ),
              shape = 15
            ) +
            geom_density_ridges(
              aes(x = val, y = fct_reorder(var, val, median)),
              rel_min_height = alpha.min,
              fill = 'grey50',
              alpha = 0.5
            ) +
            coord_cartesian(xlim = xlim) +
            scale_color_brewer(name = 'Group', palette = 'Paired') +
            scale_fill_brewer(name = 'Group', palette = 'Paired') +
            xlab('Scaled Value') +
            ylab('') +
            ggtitle(fname)
        }
      } else {
        p.ridge <- ggplot(df.scaled) +
          geom_tile(
            data = c,
            aes_string(
              x = 'cntr',
              y = 'var',
              color = 'cls',
              fill = 'cls',
              width = width,
              group = 'cls'
            ),
            alpha = 0.4,
            height = 0.2
          ) +
          geom_point(
            data = c,
            aes(
              x = cntr,
              y = var,
              fill = cls,
              color = cls,
              group = cls
            ),
            shape = 15
          ) +
          geom_density_ridges(
            aes(x = val, y = fct_reorder(var, val, median)),
            rel_min_height = alpha.min,
            fill = 'grey50',
            alpha = 0.5
          ) +
          coord_cartesian(xlim = xlim) +
          scale_color_brewer(name = 'Group', palette = 'Paired') +
          scale_fill_brewer(name = 'Group', palette = 'Paired') +
          xlab('Scaled Value') +
          ylab('') +
          ggtitle(fname)
      }
    } else {
      p.ridge <- ggplot(df.scaled) +
        geom_density_ridges(
          aes(x = val, fct_reorder(var, val, median)),
          rel_min_height = alpha.min,
          fill = 'deeppink4',
          alpha = 0.5
        ) +
        coord_cartesian(xlim = xlim) +
        xlab('Scaled Value') +
        ylab('') +
        ggtitle(fname)
    }
    p.strip <- ggplot(df.scaled) +
      stat_density(
        aes(
          x = val,
          y = fct_reorder(var, val, median),
          fill = ..density..
        ),
        geom = 'raster',
        position = 'identity',
        adjust = bw
      ) +
      scale_fill_viridis_c(
        option = 'magma',
        limits = c(0, dlim),
        begin = 0,
        end = 1,
        na.value = 'lemonchiffon'
      ) +
      coord_cartesian(xlim = xlim) +
      xlab('Scaled Value') +
      ylab('') +
      ggtitle(fname) +
      theme(
        panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.major = element_line(color = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.minor = element_line(color = rgb(0.00146, 0.000466, 0.0139))
      )
    p <- list(p.ridge = p.ridge, p.strip = p.strip)
  }
  return(p)
}
# p.twod ----
# Visualize gaussian density of features in 2D
p.twod <- function(df,
                   features = 'all',
                   xlim = c(-2, 2),
                   dlim = 0.05) {
  fname <- df$model[1]
  if (features == 'all') {
    df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames() ->
      features
  } else {
    features <- features
  }
  f <- vector(mode = 'list', length = length(features))
  for (i in 1:length(features)) {
    v <- features[i]
    df.scaled <- df %>%
      ungroup() %>%
      select(all_of(features)) %>%
      pivot_longer(-all_of(v),
                   names_to = 'var',
                   values_to = 'val') %>%
      group_by(var) %>%
      filter(!is.na(val)) %>%
      nest() %>%
      mutate(d.scaled = purrr::map(data, ~ as_tibble(scale(.x)))) %>%
      select(var, d.scaled) %>%
      unnest(cols = d.scaled) %>%
      add_column(model = fname)
    p <- df.scaled %>%
      ggplot() +
      geom_bin2d(
        aes_string(x = 'val', y = v, fill = '..density..'),
        bins = 30,
        na.rm = TRUE
      ) +
      scale_fill_viridis_c(
        option = 'magma',
        limits = c(0, dlim),
        begin = 0,
        end = 1,
        na.value = 'lemonchiffon'
      ) +
      coord_cartesian(xlim = xlim,
                      ylim = xlim) +
      xlab('Scaled Value') +
      theme(
        panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      facet_wrap(~ var)
    f[[i]] <- p
  }
  names(f) <- features
  w <-
    wrap_plots(f, guides = 'collect') + plot_annotation(title = fname)
  return(w)
}
# f.oned ----
# visualize features in 1D for many models
f.oned <-
  function(lst,
           mods = NULL,
           mod.type = NULL,
           dr = FALSE,
           sigma = 2,
           runs = 'all',
           ncol = 4,
           nrow = 3,
           features = 'all',
           plot = 'all',
           xlim = c(-2, 2),
           dlim = 0.5,
           alpha.min = 0.02,
           bw = 0.5) {
    if (sigma == 2) {
      width <- 'twosigma'
    } else if (sigma == 3) {
      width <- 'threesigma'
    } else {
      width <- 'twosigma'
    }
    if (runs == 'all') {
      lst <- lst
      mods <- mods
    } else {
      lst <- lst %>% keep(!(names(lst) %in% runs))
      if (!is.null(mods)) {
        mods <- mods %>% keep(!(names(mods) %in% runs))
      }
    }
    if (features == 'all') {
      df %>%
        ungroup() %>%
        select(where(is.numeric)) %>%
        colnames() ->
        features
    } else {
      features <- features
    }
    df %>%
      ungroup() %>%
      group_by(model) %>%
      select(all_of(features)) %>%
      pivot_longer(cols = -model,
                   names_to = 'var',
                   values_to = 'val') %>%
      group_by(model, var) %>%
      filter(!is.na(val)) %>%
      nest() %>%
      mutate(d.scaled = purrr::map(data, ~ as_tibble(scale(.x)))) %>%
      select(var, d.scaled) %>%
      unnest(cols = d.scaled) ->
      df.scaled
    if (!is.null(mods)) {
      c <- purrr::map_dfr(
        mods,
        .id = 'model',
        .f = function(mod) {
          if (!is.null(mod.type)) {
            if (mod.type == 'dr' & dr == TRUE) {
              ddr <- mod$dir %>%
                as_tibble() %>%
                pivot_longer(cols = everything(),
                             names_to = 'var',
                             values_to = 'val')
              cntr <- mod$dir %>%
                as_tibble() %>%
                add_column(cls = mod$classification) %>%
                pivot_longer(
                  cols = -cls,
                  names_to = 'var',
                  values_to = 'val'
                ) %>%
                mutate(across(where(is.character), factor)) %>%
                group_by(cls, var) %>%
                summarise(
                  cntr = mean(val),
                  sigma = sd(val),
                  twosigma = sd(val) * 2,
                  threesigma = sd(val) * 3
                )
              c <- ddr %>% left_join(cntr) %>% distinct()
            } else {
              cntr <- as_tibble(mod$mu, rownames = 'var') %>%
                rename_with(~ gsub('V', '', .x), .cols = where(is.numeric)) %>%
                mutate(var = factor(features)) %>%
                pivot_longer(
                  cols = where(is.numeric),
                  names_to = 'cls',
                  values_to = 'cntr'
                ) %>%
                mutate(across(where(is.character), factor))
              sig <- apply(mod$sigma, 3, diag) %>%
                as_tibble(rownames = 'var') %>%
                rename_with(~ gsub('V', '', .x), .cols = where(is.numeric)) %>%
                pivot_longer(
                  cols = where(is.numeric),
                  names_to = 'cls',
                  values_to = 'variance'
                ) %>%
                mutate(
                  across(where(is.character), factor),
                  twosigma = sqrt(variance) * 2,
                  threesigma = sqrt(variance) * 3
                )
              c <- cntr %>% left_join(sig)
            }
          } else {
            cntr <-
              as_tibble(mod$parameters$mean,
                        rownames = 'var',
                        .name_repair = 'unique') %>%
              rename_with(~ gsub('...', '', .x), .cols = where(is.numeric)) %>%
              pivot_longer(
                cols = where(is.numeric),
                names_to = 'cls',
                values_to = 'cntr'
              ) %>%
              mutate(across(where(is.character), factor))
            sig <- apply(mod$parameters$variance$sigma, 3, diag) %>%
              as_tibble(rownames = 'var') %>%
              rename_with(~ gsub('V', '', .x), .cols = where(is.numeric)) %>%
              pivot_longer(
                cols = where(is.numeric),
                names_to = 'cls',
                values_to = 'variance'
              ) %>%
              mutate(
                across(where(is.character), factor),
                twosigma = sqrt(variance) * 2,
                threesigma = sqrt(variance) * 3
              )
            c <- cntr %>% left_join(sig)
          }
        }
      ) %>%
        mutate(across(where(is.character), factor))
    }
    n.pg <- ceiling(length(unique(df.scaled$model)) / (ncol * nrow))
    f <- vector(mode = 'list', length = length(n.pg))
    f.strip <- vector(mode = 'list', length = length(n.pg))
    f.ridge <- vector(mode = 'list', length = length(n.pg))
    if (plot == 'strip') {
      p <- ggplot(df.scaled) +
        stat_density(
          aes(
            x = val,
            fct_reorder(var, val, median),
            fill = ..density..
          ),
          geom = 'raster',
          position = 'identity',
          adjust = bw
        ) +
        scale_fill_viridis_c(
          option = 'magma',
          limits = c(0, dlim),
          begin = 0,
          end = 1,
          na.value = 'lemonchiffon'
        ) +
        coord_cartesian(xlim = xlim) +
        xlab('Scaled Value') +
        ylab('') +
        theme(
          panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
          panel.grid.major = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7)),
          panel.grid.minor = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7))
        )
      for (i in 1:n.pg) {
        f[[i]] <- p + facet_wrap_paginate(~ model,
                                          ncol = ncol,
                                          nrow = nrow,
                                          page = i)
      }
    } else if (plot == 'ridge') {
      if (!is.null(mods)) {
        if (!is.null(mod.type)) {
          if (mod.type == 'dr' & dr == TRUE) {
            cntr <- c %>% group_by(model, cls) %>% select(-val) %>% distinct()
            ddr <-
              c %>% group_by(var) %>% select(model, var, val) %>% distinct()
            p <- ggplot() +
              geom_tile(
                data = cntr,
                aes_string(
                  x = 'cntr',
                  y = 'var',
                  color = 'cls',
                  fill = 'cls',
                  width = width,
                  group = 'cls'
                ),
                alpha = 0.4,
                height = 0.2
              ) +
              geom_point(
                data = cntr,
                aes(
                  x = cntr,
                  y = var,
                  fill = cls,
                  color = cls,
                  group = cls
                ),
                shape = 15
              ) +
              geom_density_ridges(
                data = ddr,
                aes(x = val, y = fct_reorder(var, val, median)),
                rel_min_height = alpha.min,
                fill = 'grey50',
                alpha = 0.5
              ) +
              coord_cartesian(xlim = xlim) +
              scale_color_brewer(name = 'Group', palette = 'Paired') +
              scale_fill_brewer(name = 'Group', palette = 'Paired') +
              xlab('Scaled Value') +
              ylab('')
          } else {
            p <- ggplot(df.scaled) +
              geom_tile(
                data = c,
                aes_string(
                  x = 'cntr',
                  y = 'var',
                  color = 'cls',
                  fill = 'cls',
                  width = width,
                  group = 'cls'
                ),
                alpha = 0.4,
                height = 0.2
              ) +
              geom_point(
                data = c,
                aes(
                  x = cntr,
                  y = var,
                  fill = cls,
                  color = cls,
                  group = cls
                ),
                shape = 15
              ) +
              geom_density_ridges(
                aes(x = val, y = fct_reorder(var, val, median)),
                rel_min_height = alpha.min,
                fill = 'grey50',
                alpha = 0.5
              ) +
              coord_cartesian(xlim = xlim) +
              scale_color_brewer(name = 'Group', palette = 'Paired') +
              scale_fill_brewer(name = 'Group', palette = 'Paired') +
              xlab('Scaled Value') +
              ylab('')
          }
        } else {
          p <- ggplot(df.scaled) +
            geom_tile(
              data = c,
              aes_string(
                x = 'cntr',
                y = 'var',
                color = 'cls',
                fill = 'cls',
                width = width,
                group = 'cls'
              ),
              alpha = 0.4,
              height = 0.2
            ) +
            geom_point(
              data = c,
              aes(
                x = cntr,
                y = var,
                fill = cls,
                color = cls,
                group = cls
              ),
              shape = 15
            ) +
            geom_density_ridges(
              aes(x = val, y = fct_reorder(var, val, median)),
              rel_min_height = alpha.min,
              fill = 'grey50',
              alpha = 0.5
            ) +
            coord_cartesian(xlim = xlim) +
            scale_color_brewer(name = 'Group', palette = 'Paired') +
            scale_fill_brewer(name = 'Group', palette = 'Paired') +
            xlab('Scaled Value') +
            ylab('')
        }
      } else {
        p <- ggplot(df.scaled) +
          geom_density_ridges(
            aes(x = val, fct_reorder(var, val, median)),
            rel_min_height = alpha.min,
            fill = 'deeppink4',
            alpha = 0.5
          ) +
          coord_cartesian(xlim = xlim) +
          xlab('Scaled Value') +
          ylab('')
      }
      for (i in 1:n.pg) {
        f[[i]] <- p + facet_wrap_paginate(~ model,
                                          ncol = ncol,
                                          nrow = nrow,
                                          page = i)
      }
    } else if (plot == 'all') {
      if (!is.null(mods)) {
        if (!is.null(mod.type)) {
          if (mod.type == 'dr' & dr == TRUE) {
            cntr <- c %>% group_by(model, cls) %>% select(-val) %>% distinct()
            ddr <-
              c %>% group_by(var) %>% select(model, var, val) %>% distinct()
            p.ridge <- ggplot() +
              geom_tile(
                data = cntr,
                aes_string(
                  x = 'cntr',
                  y = 'var',
                  color = 'cls',
                  fill = 'cls',
                  width = width,
                  group = 'cls'
                ),
                alpha = 0.4,
                height = 0.2
              ) +
              geom_point(
                data = cntr,
                aes(
                  x = cntr,
                  y = var,
                  fill = cls,
                  color = cls,
                  group = cls
                ),
                shape = 15
              ) +
              geom_density_ridges(
                data = ddr,
                aes(x = val, y = fct_reorder(var, val, median)),
                rel_min_height = alpha.min,
                fill = 'grey50',
                alpha = 0.5
              ) +
              coord_cartesian(xlim = xlim) +
              scale_color_brewer(name = 'Group', palette = 'Paired') +
              scale_fill_brewer(name = 'Group', palette = 'Paired') +
              xlab('Scaled Value') +
              ylab('')
          } else {
            p.ridge <- ggplot(df.scaled) +
              geom_tile(
                data = c,
                aes_string(
                  x = 'cntr',
                  y = 'var',
                  color = 'cls',
                  fill = 'cls',
                  width = width,
                  group = 'cls'
                ),
                alpha = 0.4,
                height = 0.2
              ) +
              geom_point(
                data = c,
                aes(
                  x = cntr,
                  y = var,
                  fill = cls,
                  color = cls,
                  group = cls
                ),
                shape = 15
              ) +
              geom_density_ridges(
                aes(x = val, y = fct_reorder(var, val, median)),
                rel_min_height = alpha.min,
                fill = 'grey50',
                alpha = 0.5
              ) +
              coord_cartesian(xlim = xlim) +
              scale_color_brewer(name = 'Group', palette = 'Paired') +
              scale_fill_brewer(name = 'Group', palette = 'Paired') +
              xlab('Scaled Value') +
              ylab('')
          }
        } else {
          p.ridge <- ggplot(df.scaled) +
            geom_tile(
              data = c,
              aes_string(
                x = 'cntr',
                y = 'var',
                color = 'cls',
                fill = 'cls',
                width = width,
                group = 'cls'
              ),
              alpha = 0.4,
              height = 0.2
            ) +
            geom_point(
              data = c,
              aes(
                x = cntr,
                y = var,
                fill = cls,
                color = cls,
                group = cls
              ),
              shape = 15
            ) +
            geom_density_ridges(
              aes(x = val, y = fct_reorder(var, val, median)),
              rel_min_height = alpha.min,
              fill = 'grey50',
              alpha = 0.5
            ) +
            coord_cartesian(xlim = xlim) +
            scale_color_brewer(name = 'Group', palette = 'Paired') +
            scale_fill_brewer(name = 'Group', palette = 'Paired') +
            xlab('Scaled Value') +
            ylab('') +
            ggtitle(fname)
        }
      } else {
        p.ridge <- ggplot(df.scaled) +
          geom_density_ridges(
            aes(x = val, fct_reorder(var, val, median)),
            rel_min_height = alpha.min,
            fill = 'deeppink4',
            alpha = 0.5
          ) +
          coord_cartesian(xlim = xlim) +
          xlab('Scaled Value') +
          ylab('')
      }
      p.strip <- ggplot(df.scaled) +
        stat_density(
          aes(
            x = val,
            fct_reorder(var, val, median),
            fill = ..density..
          ),
          geom = 'raster',
          position = 'identity',
          adjust = bw
        ) +
        scale_fill_viridis_c(
          option = 'magma',
          limits = c(0, dlim),
          begin = 0,
          end = 1,
          na.value = 'lemonchiffon'
        ) +
        coord_cartesian(xlim = xlim) +
        xlab('Scaled Value') +
        ylab('') +
        theme(
          panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
          panel.grid.major = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7)),
          panel.grid.minor = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7))
        )
      p <- list(p.ridge = p.ridge, p.strip = p.strip)
      for (i in 1:n.pg) {
        f.strip[[i]] <- p$p.strip + facet_wrap_paginate( ~ model,
                                                         ncol = ncol,
                                                         nrow = nrow,
                                                         page = i)
        f.ridge[[i]] <- p$p.ridge + facet_wrap_paginate( ~ model,
                                                         ncol = ncol,
                                                         nrow = nrow,
                                                         page = i)
      }
      f <- list(f.strip = f.strip, f.ridge = f.ridge)
    }
    return(f)
  }
# f.summary ----
# Combine gifs into one faceted summary gif
f.summary <- function(lst,
                      mods,
                      runs = 'all',
                      GIF = c('xy', 'PT'),
                      ncol = 4,
                      nrow = 3,
                      grads = FALSE,
                      save = TRUE,
                      fname) {
  if (runs == 'all') {
    lst <- lst
    mods <- mods
  } else {
    lst <- lst %>% keep(!(names(lst) %in% runs))
    mods <- mods %>% keep(!(names(mods) %in% runs))
  }
  d <- purrr::map2(
    lst,
    mods,
    ~ .y$classification %>%
      as_tibble() %>%
      rename(cls = value) %>%
      add_column(markerID = factor(unique(.x$markerID))) %>%
      left_join(.x) %>%
      mutate(across(contains('cls'), factor)) %>%
      bind_rows
  )
  n.pg <- ceiling(length(unique(df$model)) / (ncol * nrow))
  if (GIF == 'xy') {
    anim <- df %>% gif.xy()
  } else {
    anim <- df %>% gif.pt()
  }
  if (grads == TRUE) {
    for (i in seq_along(n.pg)) {
      p <- anim +
        facet_wrap_paginate(~ model,
                            nrow = 3,
                            ncol = 4,
                            page = i) +
        geom_abline(size = 0.1,
                    intercept = 0,
                    slope = 1 / 18.5) +
        geom_abline(size = 0.1,
                    intercept = 0,
                    slope = 1 / 37.0) +
        geom_abline(size = 0.1,
                    intercept = 0,
                    slope = 1 / 74.0)
    }
  } else {
    for (i in seq_along(n.pg)) {
      p <- anim +
        facet_wrap_paginate(~ model,
                            nrow = 3,
                            ncol = 4,
                            page = i)
    }
  }
  if (save == TRUE) {
    anim_save(
      paste0(fname, '.', i, '.gif'),
      animate(
        p,
        height = 8.5,
        width = 11,
        units = 'in',
        res = 300
      )
    )
  }
  return(p)
}
# gif.xy ----
# Marker Position Plot
gif.xy <- function(df) {
  p.pos <- ggplot(df) +
    geom_point(aes(
      xposkm,
      yposkm,
      alpha = time,
      size = Pkbar,
      group = markerID,
      color = cls
    ),
    na.rm = TRUE) +
    ylim(130, 0) +
    xlab('Distance (km)') +
    ylab('Depth (km)') +
    scale_alpha_continuous(range = c(0.01, 1),
                           guide = 'none') +
    scale_size_continuous(
      range = c(0.01, 0.25),
      limits = c(0, 40),
      guide = 'none'
    ) +
    scale_color_brewer(name = '', palette = 'Paired') +
    labs(title = 'Marker Positions | Time: {round(frame_time/10^6, 2)}Ma') +
    transition_time(time) +
    ease_aes('linear') +
    shadow_wake(wake_length = 0.1, alpha = FALSE) +
    theme_bw(base_size = 14) +
    theme(
      axis.ticks.length = unit(-0.25, "cm"),
      axis.text.x = element_text(margin = unit(c(0.5, 0, 0, 0), "cm")),
      axis.text.y = element_text(margin = unit(c(0, 0.5, 0, 0), "cm")),
      panel.border = element_rect(size = 1, color = 'black'),
      axis.text = element_text(face = 'plain', color = 'black'),
      axis.ticks = element_line(size = 0.5, color = 'black'),
      legend.direction = 'horizontal',
      legend.justification = c(0, 1),
      legend.position = c(0, 0),
      legend.box.margin = margin(c(50, 50, 50, 50)),
      axis.title = element_text(size = 12, face = 'plain'),
      plot.title = element_text(size = 14, face = 'plain')
    )
  return(p.pos)
}

# gif.pt ----
# Marker PT Plot
gif.pt <- function(df) {
  p.pt <- ggplot(df) +
    geom_point(aes(
      TCel,
      Pkbar,
      alpha = time,
      size = Pkbar,
      group = markerID,
      color = cls
    ),
    na.rm = TRUE) +
    xlim(0, 1000) +
    ylim(0, 40) +
    xlab(expression('Temperature' ~ (degree * C))) +
    ylab('Pressure (kbar)') +
    scale_alpha_continuous(range = c(0.01, 1),
                           guide = 'none') +
    scale_size_continuous(
      range = c(0.01, 0.25),
      limits = c(0, 40),
      guide = 'none'
    ) +
    scale_color_brewer(name = '', palette = 'Paired') +
    labs(title = 'Marker PT Paths | Time: {round(frame_time/10^6, 2)}Ma') +
    transition_time(time) +
    ease_aes('linear') +
    shadow_wake(wake_length = 0.1, alpha = FALSE) +
    theme_bw(base_size = 14) +
    theme(
      axis.ticks.length = unit(-0.25, "cm"),
      axis.text.x = element_text(margin = unit(c(0.5, 0, 0, 0), "cm")),
      axis.text.y = element_text(margin = unit(c(0, 0.5, 0, 0), "cm")),
      panel.border = element_rect(size = 1, color = 'black'),
      axis.text = element_text(face = 'plain', color = 'black'),
      axis.ticks = element_line(size = 0.5, color = 'black'),
      legend.direction = 'horizontal',
      legend.justification = c(-0.2, -0.5),
      legend.position = c(0, 0),
      axis.title = element_text(size = 12, face = 'plain'),
      plot.title = element_text(size = 14, face = 'plain')
    )
  return(p.pt)
}
# a.oned ----
# Summarises one-dimension data by flicking through models
a.oned <- function(lst,
                   mods,
                   mod.type = NULL,
                   dr = FALSE,
                   sigma = 2,
                   features = 'all',
                   plot = 'all',
                   xlim = c(-2, 2),
                   dlim = 1,
                   alpha.min = 0.02,
                   bw = 0.5,
                   save = TRUE,
                   fps = 20,
                   type = 'gif',
                   fname = 'oned') {
  df <- lst %>% bind_rows()
  if (features == 'all') {
    df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames() ->
      features
  } else {
    features <- features
  }
  df %>%
    ungroup() %>%
    select(all_of(features)) %>%
    pivot_longer(cols = everything(),
                 names_to = 'var',
                 values_to = 'val') %>%
    group_by(var) %>%
    filter(!is.na(val)) %>%
    nest() %>%
    mutate(d.scaled = purrr::map(data, ~ as_tibble(scale(.x)))) %>%
    select(var, d.scaled) %>%
    unnest(cols = d.scaled) %>%
    add_column(model = fname) ->
    df.scaled
  if (!is.null(mods)) {
    c <- purrr::map_dfr(mods, function(mod) {
      cntr <-
        as_tibble(mod$parameters$mean,
                  rownames = 'var',
                  .name_repair = 'unique') %>%
        rename_with( ~ gsub('...', '', .x), .cols = where(is.numeric)) %>%
        pivot_longer(
          cols = where(is.numeric),
          names_to = 'cls',
          values_to = 'cntr'
        ) %>%
        mutate(across(where(is.character), factor))
      sig <- apply(mod$parameters$variance$sigma, 3, diag) %>%
        as_tibble(rownames = 'var') %>%
        rename_with( ~ gsub('V', '', .x), .cols = where(is.numeric)) %>%
        pivot_longer(
          cols = where(is.numeric),
          names_to = 'cls',
          values_to = 'variance'
        ) %>%
        mutate(
          across(where(is.character), factor),
          twosigma = sqrt(variance) * 2,
          threesigma = sqrt(variance) * 3
        )
      c <- cntr %>% left_join(sig)
    })
  }
  if (plot == 'strip') {
    a <- ggplot(df.scaled) +
      stat_density(
        aes(
          x = val,
          fct_reorder(var, val, median),
          fill = ..density..
        ),
        color = 'lemonchiffon',
        geom = 'raster',
        position = 'identity',
        adjust = bw
      ) +
      scale_fill_viridis_c(
        option = 'magma',
        limits = c(0, dlim),
        begin = 0,
        end = 1,
        na.value = 'lemonchiffon'
      ) +
      coord_cartesian(xlim = xlim) +
      xlab('Scaled Value') +
      ylab('') +
      labs(title = '{closest_state}') +
      theme(
        panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.major = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7)),
        panel.grid.minor = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7))
      ) +
      transition_states(df.scaled$model) +
      ease_aes('linear') +
      enter_fade() +
      exit_fade()
  } else if (plot == 'ridge') {
    if (!is.null(mods)) {
      if (sigma == 2) {
        width <- 'twosigma'
      } else if (sigma == 3) {
        width <- 'threesigma'
      } else {
        width <- 'twosigma'
      }
      a <- ggplot(df.scaled) +
        geom_tile(
          data = c,
          aes_string(
            x = 'cntr',
            y = 'var',
            color = 'cls',
            fill = 'cls',
            width = width,
            group = 'cls'
          ),
          alpha = 1,
          height = 0.2
        ) +
        geom_density_ridges(
          aes(x = val, y = fct_reorder(var, val, median)),
          rel_min_height = alpha.min,
          fill = 'deeppink4',
          alpha = 0.5
        ) +
        coord_cartesian(xlim = xlim) +
        scale_color_brewer(name = 'Group', palette = 'Paired') +
        scale_fill_brewer(name = 'Group', palette = 'Paired') +
        xlab('Scaled Value') +
        ylab('') +
        transition_states(df.scaled$model) +
        labs(title = '{closest_state}') +
        ease_aes('linear') +
        enter_fade() +
        exit_fade()
    } else {
      a <- ggplot(df.scaled) +
        geom_density_ridges(
          aes(x = val, fct_reorder(var, val, median)),
          rel_min_height = alpha.min,
          fill = 'deeppink4',
          alpha = 0.5
        ) +
        coord_cartesian(xlim = xlim) +
        xlab('Scaled Value') +
        ylab('') +
        transition_states(df.scaled$model) +
        labs(title = '{closest_state}') +
        ease_aes('linear') +
        enter_fade() +
        exit_fade()
    }
  } else if (plot == 'all') {
    if (!is.null(mods)) {
      if (sigma == 2) {
        width <- 'twosigma'
      } else if (sigma == 3) {
        width <- 'threesigma'
      } else {
        width <- 'twosigma'
      }
      a.ridge <- ggplot(df.scaled) +
        geom_tile(
          data = c,
          aes_string(
            x = 'cntr',
            y = 'var',
            color = 'cls',
            fill = 'cls',
            width = width,
            group = 'cls'
          ),
          alpha = 1,
          height = 0.2
        ) +
        geom_density_ridges(
          aes(x = val, y = fct_reorder(var, val, median)),
          rel_min_height = alpha.min,
          fill = 'deeppink4',
          alpha = 0.5
        ) +
        coord_cartesian(xlim = xlim) +
        scale_color_brewer(name = 'Group', palette = 'Paired') +
        scale_fill_brewer(name = 'Group', palette = 'Paired') +
        xlab('Scaled Value') +
        ylab('') +
        transition_states(df.scaled$model) +
        labs(title = '{closest_state}') +
        ease_aes('linear') +
        enter_fade() +
        exit_fade()
    } else {
      a.ridge <- ggplot(df.scaled) +
        geom_density_ridges(
          aes(x = val, fct_reorder(var, val, median)),
          rel_min_height = alpha.min,
          fill = 'deeppink4',
          alpha = 0.5
        ) +
        coord_cartesian(xlim = xlim) +
        xlab('Scaled Value') +
        ylab('') +
        transition_states(df.scaled$model) +
        labs(title = '{closest_state}') +
        ease_aes('linear') +
        enter_fade() +
        exit_fade()
    }
    a.strip <- ggplot(df.scaled) +
      stat_density(
        aes(
          x = val,
          fct_reorder(var, val, median),
          fill = ..density..
        ),
        geom = 'raster',
        position = 'identity',
        adjust = bw
      ) +
      scale_fill_viridis_c(
        option = 'magma',
        limits = c(0, dlim),
        begin = 0,
        end = 1,
        na.value = 'lemonchiffon'
      ) +
      coord_cartesian(xlim = xlim) +
      xlab('Scaled Value') +
      ylab('') +
      labs(title = '{closest_state}') +
      theme(
        panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.major = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7)),
        panel.grid.minor = element_line(color = rgb(0.00146, 0.000466, 0.0139, 0.7))
      ) +
      transition_states(df.scaled$model) +
      ease_aes('linear') +
      enter_fade() +
      exit_fade()
    a <- list(a.ridge = a.ridge, a.strip = a.strip)
  }
  if (save == TRUE) {
    if (type == 'gif') {
      if (length(a) > 1) {
        anim_save(
          paste0(fname, '.strip.gif'),
          animate(
            a$a.strip,
            height = 4.25,
            width = 5.5,
            units = 'in',
            res = 300,
            fps = fps
          )
        )
        anim_save(
          paste0(fname, '.ridge.gif'),
          animate(
            a$a.ridge,
            height = 4.25,
            width = 5.5,
            units = 'in',
            res = 300,
            fps = fps
          )
        )
      } else {
        anim_save(
          paste0(fname, '.gif'),
          animate(
            a$a.ridge,
            height = 4.25,
            width = 5.5,
            units = 'in',
            res = 300,
            fps = fps
          )
        )
      }
    } else if (type == 'avi') {
      if (length(a) > 1) {
        anim_save(
          paste0(fname, '.strip.avi'),
          animate(
            a$a.strip,
            height = 4.25,
            width = 5.5,
            units = 'in',
            res = 300,
            fps = fps,
            renderer = av_renderer()
          )
        )
        anim_save(
          paste0(fname, '.ridge.avi'),
          animate(
            a$a.ridge,
            height = 4.25,
            width = 5.5,
            units = 'in',
            res = 300,
            fps = fps,
            renderer = av_renderer()
          )
        )
      } else {
        anim_save(
          paste0(fname, '.avi'),
          animate(
            a$a.ridge,
            height = 4.25,
            width = 5.5,
            units = 'in',
            res = 300,
            fps = fps,
            renderer = av_renderer()
          )
        )
      }
    }
  }
  return(a)
}
# a.twod ----
# Summarises two-dimension data by flicking through models
a.twod <- function(lst,
                   features = 'all',
                   xlim = c(-2, 2),
                   dlim = 0.05,
                   save = TRUE,
                   fps = 20,
                   type = 'gif',
                   fname = 'twod') {
  df <- lst %>% bind_rows()
  if (features == 'all') {
    df %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      colnames() ->
      features
  } else {
    features <- features
  }
  anims <- vector(mode = 'list', length = length(features))
  for (i in 1:length(features)) {
    v <- features[i]
    df.scaled <- df %>%
      group_by(model) %>%
      select(all_of(features)) %>%
      pivot_longer(-c(all_of(v), model),
                   names_to = 'var',
                   values_to = 'val') %>%
      group_by(model, var) %>%
      filter(!is.na(val)) %>%
      nest() %>%
      mutate(d.scaled = purrr::map(data, ~ as_tibble(scale(.x)))) %>%
      select(var, d.scaled) %>%
      unnest(cols = d.scaled) %>%
      ungroup()
    a <- df.scaled %>%
      ggplot() +
      geom_bin2d(
        aes_string(x = 'val', y = v, fill = '..density..'),
        bins = 30,
        na.rm = TRUE
      ) +
      scale_fill_viridis_c(
        option = 'magma',
        limits = c(0, dlim),
        begin = 0,
        end = 1,
        na.value = 'lemonchiffon'
      ) +
      coord_cartesian(xlim = xlim,
                      ylim = xlim) +
      xlab('Scaled Value') +
      labs(title = '{closest_state}') +
      transition_states(df.scaled$model) +
      ease_aes('linear') +
      enter_fade() +
      exit_fade() +
      theme(
        panel.background = element_rect(fill = rgb(0.00146, 0.000466, 0.0139)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      facet_wrap(~ var)
    anims[[i]] <- a
  }
  names(anims) <- features
  w <-
    wrap_plots(a, guides = 'collect')
  if (save == TRUE) {
    if (type == 'gif') {
      anim_save(
        filename = paste0('facet.', .x$labels$y, '.gif'),
        animation = animate(
          w,
          height = 8.5,
          width = 11,
          units = 'in',
          res = 300,
          fps = fps
        )
      )
    } else if (type == 'avi') {
      anim_save(
        filename = paste0('facet.', .x$labels$y, '.avi'),
        animation = animate(
          w,
          height = 8.5,
          width = 11,
          units = 'in',
          res = 300,
          fps = fps,
          renderer = av_renderer()
        )
      )
    }
  }
  return(w)
}
