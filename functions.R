# Load packages
# Quiet loading
sshhh <- function(p){
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, character.only=TRUE)))}

# Package list
c('magrittr', 'tidyr', 'readr', 'purrr',
  'mclust', 'ggforce', 'dplyr',
  'patchwork', 'gridExtra', 'gganimate',
  'ggridges', 'progress',
  'parallel') -> p.list

cat('Loading libraries:', p.list, sep = '\n')

# auto-load quietly
sapply(p.list, sshhh)

cat('Loading functions\n')

# Read binary (.prn) files and trace markers
trace_marx <- function(
  prn.paths,
  marx.est = 20000,
  area = c(500000, 1260000, 17500, 28500),
  markers = TRUE,
  grid = TRUE){
  if(grid == FALSE & markers == FALSE) {
    stop('grid and markers cannot both be FALSE')
  }
  # List prn files
  fpaths <- prn.paths
  forder <- order(
    fpaths %>%
    purrr::map_int(~.x %>% stringr::str_extract('[0-9]+.prn+') %>%
    stringr::str_extract('[0-9]+') %>%
    as.integer()))
  fpaths.ordered <- fpaths[forder]
  fnames <- fpaths.ordered %>% stringr::str_extract('cd[a-z][0-9]+_[0-9]+')

  # Create marker arrays (n.markers x n.tsteps)
  mmm <- matrix(NA, marx.est) # marker global index
  mty <- matrix(NA, marx.est, length(fpaths)) # type
  mti <- matrix(NA, marx.est, length(fpaths)) # time, [yr]
  mxx <- matrix(NA, marx.est, length(fpaths)) # x, [m]
  mzz <- matrix(NA, marx.est, length(fpaths)) # y, [m]
  mtk <- matrix(NA, marx.est, length(fpaths)) # T, [K]
  mpb <- matrix(NA, marx.est, length(fpaths)) # P, [bar]

  # Create grid array
  grids <- vector('list', length(fpaths))

  # Coordinates of the sampling area, [m]
  xmin <- area[1]
  xmax <- area[2]
  zmin <- area[3]
  zmax <- area[4]

  # Counters
  mfind <- 0 # Marker counter
  g.count <- 1 # Grid counter
  tstep <- 1 # .prn (tstep) counter

  for(path in fpaths.ordered){
    # Filename
    f <- path
    # Model name
    f.name <- stringr::str_extract(f, 'cd.[0-9]+.[0-9]+')
    # Open connection
    f.prn <- file(f, 'rb')
    # Read sizes of variables
    readBin(f.prn, 'integer', 4, 1, signed = F)
    # Read model parameters
    # Grid resolution
    xnumx <- readBin(f.prn, 'integer', 1, 8)
    znumz <- readBin(f.prn, 'integer', 1, 8)
    # Markers per cell
    mnumx <- readBin(f.prn, 'integer', 1, 8)
    mnumz <- readBin(f.prn, 'integer', 1, 8)
    # Number of markers
    marknum <- readBin(f.prn, 'integer', 1, 8)
    # Model sizes
    xsize <- readBin(f.prn, 'numeric', 1, 8)
    zsize <- readBin(f.prn, 'numeric', 1, 8)
    # Pressure value
    pinit <- readBin(f.prn, 'numeric', 5, 8)
    # Gravity
    gx <- readBin(f.prn, 'numeric', 1, 8)
    gz <- readBin(f.prn, 'numeric', 1, 8)
    # Number of rocks
    rocknum <- readBin(f.prn, 'integer', 1, 4)
    # Number of Boundary conditions
    boundnum <- readBin(f.prn, 'integer', 1, 8)
    # Stage, time
    stg <- readBin(f.prn, 'integer', 1, 4)
    timesum <- readBin(f.prn, 'numeric', 1, 8)
    # Skip rock properties
    curpos <- 4+2*4+16*8+rocknum*(8*24+4)
    seek(f.prn, curpos, 'start')

    if(grid == FALSE & markers == TRUE) {
      # Initialize matrices
      pr <- matrix(NA, znumz, xnumx)
      # Progress bar
      pb.nodes <- progress_bar$new(
        format = paste0(
          'Reading Pressure Nodes [',
          f.name,
          '] [:bar] :percent in: :elapsed'),
        total = xnumx*znumz,
        clear = FALSE,
        width = 100)
      # Read nodes information
      for(i in seq_len(xnumx)){
        for(j in seq_len(znumz)){
          vbuf <- readBin(f.prn, 'numeric', 3, 4)
          pr[j,i] <- vbuf[1]
          pb.nodes$tick()
        }
      }
    } else if(grid == TRUE){
      # Initialize Matrices
      pr <- matrix(NA, znumz, xnumx)
      vx <- matrix(NA, znumz, xnumx)
      vz <- matrix(NA, znumz, xnumx)
      exx <- matrix(NA, znumz, xnumx)
      ezz <- matrix(NA, znumz, xnumx)
      exz <- matrix(NA, znumz, xnumx)
      sxx <- matrix(NA, znumz, xnumx)
      szz <- matrix(NA, znumz, xnumx)
      sxz <- matrix(NA, znumz, xnumx)
      ro <- matrix(NA, znumz, xnumx)
      nu <- matrix(NA, znumz, xnumx)
      nd <- matrix(NA, znumz, xnumx)
      mu <- matrix(NA, znumz, xnumx)
      ep <- matrix(NA, znumz, xnumx)
      et <- matrix(NA, znumz, xnumx)
      pr0 <- matrix(NA, znumz, xnumx)
      prb <- matrix(NA, znumz, xnumx)
      dv <- matrix(NA, znumz, xnumx)
      tk <- matrix(NA, znumz, xnumx)
      cp <- matrix(NA, znumz, xnumx)
      kt <- matrix(NA, znumz, xnumx)
      ht <- matrix(NA, znumz, xnumx)
      eii <- matrix(1, znumz, xnumx)*1e-16
      sii <- matrix(1, znumz, xnumx)*1e+4
      # Progress bar
      pb.nodes <- progress_bar$new(
        format = paste0(
          'Reading Nodes [',
          f.name,
          '] [:bar] :percent in: :elapsed'),
        total = xnumx*znumz,
        clear = FALSE,
        width = 100)
      # Read nodes information
      for(i in seq_len(xnumx)){
        for(j in seq_len(znumz)){
          vbuf <- readBin(f.prn, 'numeric', 3, 4)
          pr[j,i] <- vbuf[1]
          vx[j,i] <- vbuf[2]
          vz[j,i] <- vbuf[3]
          vbuf1 <- readBin(f.prn, 'integer', 3, 8)
          vbuf2 <- readBin(f.prn, 'numeric', 16, 4)
          exx[j,i] <- vbuf2[1]
          ezz[j,i] <- vbuf2[2]
          exz[j,i] <- vbuf2[3]
          sxx[j,i] <- vbuf2[4]
          szz[j,i] <- vbuf2[5]
          sxz[j,i] <- vbuf2[6]
          ro[j,i] <- vbuf2[7]
          nu[j,i] <- vbuf2[8]
          nd[j,i] <- vbuf2[9]
          mu[j,i] <- vbuf2[10]
          ep[j,i] <- vbuf2[11]
          et[j,i] <- vbuf2[12]
          pr0[j,i] <- vbuf2[13]
          prb[j,i] <- vbuf2[14]
          dv[j,i] <- vbuf2[15]
          tk[j,i] <- vbuf2[16]
          vbuf3 <- readBin(f.prn, 'integer', 1, 8)
          vbuf4 <- readBin(f.prn, 'numeric', 3, 4)
          cp[j,i] <- vbuf4[1]
          kt[j,i] <- vbuf4[2]
          ht[j,i] <- vbuf4[3]
          pb.nodes$tick()
        }
      }
    }

    # Skip all nodes
    curpos2 <- curpos+(4*22+8*4)*xnumx*znumz
    seek(f.prn, curpos2, 'start')
    # Read gridline positions
    gx <- readBin(f.prn, 'numeric', xnumx, 4)
    gz <- readBin(f.prn, 'numeric', znumz, 4)

    if(grid == TRUE) {
      # Progress bar
      pb.eii <- progress_bar$new(
        format =
          paste0('Calc. Stress & Strain [',
                 f.name,
                 '] [:bar] :percent in: :elapsed'),
        total = (xnumx*znumz)-2440,
        clear = FALSE,
        width = 100)
      for (i in seq_len(xnumx-2)){
        for (j in seq_len(znumz-2)){
          eii[j+1,i+1] <-
            (exz[j+1,i+1]^2+((exx[j+1,i+1]+exx[j+2,i+1]+exx[j+1,i+2]+exx[j+2,i+2])/4)^2)^0.5;
          sii[j+1,i+1] <-
            (sxz[j+1,i+1]^2+((sxx[j+1,i+1]+sxx[j+2,i+1]+sxx[j+1,i+2]+sxx[j+2,i+2])/4)^2)^0.5;
          pb.eii$tick()
        }
      }
    }

    if(markers == TRUE){
      # Skip all boundary conditions
      curpos3 <- curpos2+(xnumx+znumz)*4+(4*4+8*3)*(boundnum-1)
      # Pressure points coordinates
      px <- gx
      px[2:xnumx] <- (gx[1:xnumx-1]+gx[2:xnumx])/2
      pz <- gz
      pz[2:znumz] <- (gz[1:znumz-1]+gz[2:znumz])/2
      # Progress bar
      pb.marx <- progress_bar$new(
        format =
          paste0('Reading ',
                 marknum,
                 ' Markers [',
                 f.name,
                 '] [:bar] :percent in: :elapsed'),
        total = marknum,
        clear = FALSE,
        width = 100)
      # Find Markers
      if(mfind == 0){
        seek(f.prn, curpos3, 'start')
        for(m in seq_len(marknum)){
          mbuf <- readBin(f.prn, 'numeric', 9, 4) %>% replace(is.na(.), Inf)
          mt <- readBin(f.prn, 'integer', 1, 1, signed = F)
          mx <- mbuf[1]
          mz <- mbuf[2]
          mk <- mbuf[3]
          # Save markers from the interest area
          if(mx>=xmin & mx<=xmax & mz>=zmin & mz<=zmax & mt>1){
            mfind <- mfind+1
            mmm[mfind] <- m # global index
            mty[mfind,tstep] <- mt # type
            mti[mfind,tstep] <- timesum # time [yr]
            mxx[mfind,tstep] <- mx # x [m]
            mzz[mfind,tstep] <- mz # z [m]
            mtk[mfind,tstep] <- mk # T [K]
            # Interpolate pressure
            # Find indexes for the upper-left pressure node by bisection
            jmin <- 2
            jmax <- znumz
            while((jmax-jmin)>1){
              j <- trunc((jmax+jmin)/2)
              if(pz[j]>mz){
                jmax <- j
              } else {
                jmin <- j
              }
            }
            j <- jmin
            if(j>znumz-1){
              j <- znumz-1
            }
            imin <- 2
            imax <- xnumx
            while((imax-imin)>1){
              i <- trunc((imax+imin)/2)
              if(px[i]>mx){
                imax <- i
              } else {
                imin <- i
              }
            }
            i <- imin
            if(i>xnumx-1){
              i <- xnumx-1
            }
            # Bilinear interpolation
            dxm <- (mx-px[i])/(px[i+1]-px[i])
            dzm <- (mz-pz[j])/(pz[j+1]-pz[j])
            # P [bar]
            mpb[mfind,tstep] <-
              1e-5*((1-dxm)*(1-dzm)*pr[j,i]+
                    (dxm)*(1-dzm)*pr[j,i+1]+
                    (1-dxm)*(dzm)*pr[j+1,i]+
                    (dxm)*(dzm)*pr[j+1,i+1])
          }
          pb.marx$tick()
        }
        if(mfind == 0){
          break
        }
        # Update PT-paths of selected markers
      } else {
        for(mi in seq_len(mfind)){
          m <- mmm[mi]
          curpos4 <- curpos3+(m-1)*(9*4+1)
          seek(f.prn, curpos4, 'start')
          mbuf <- readBin(f.prn, 'numeric', 9, 4) %>% replace(is.na(.), Inf)
          mt <- readBin(f.prn, 'integer', 1, 1, signed = F)
          mx <- mbuf[1]
          mz <- mbuf[2]
          mk <- mbuf[3]
          # Save marker data
          mty[mi,tstep] <- mt # type
          mti[mi,tstep] <- timesum # time [yr]
          mxx[mi,tstep] <- mx # x [m]
          mzz[mi,tstep] <- mz # z [m]
          mtk[mi,tstep] <- mk # T [K]
          # Interpolate pressure
          # Find indexes for the upper-left pressure node by bisection
          jmin <- 2
          jmax <- znumz
          while((jmax-jmin)>1){
            j <- trunc((jmax+jmin)/2)
            if(pz[j]>mz){
              jmax <- j
            } else {
              jmin <- j
            }
          }
          j <-jmin
          if(j>znumz-1){
            j <- znumz-1
          }
          imin <- 2
          imax <- xnumx
          while((imax-imin)>1){
            i <- trunc((imax+imin)/2)
            if(px[i]>mz){
              imax <- i
            } else {
              imin <- i
            }
          }
          i <- imin
          if(i>xnumx-1){
            i <- xnumx-1
          }
          # Bilinear interpolation
          dxm <- (mx-px[i])/(px[i+1]-px[i])
          dzm <- (mz-pz[j])/(pz[j+1]-pz[j])
          # P [bar]
          mpb[mi,tstep] <-
            1e-5*((1-dxm)*(1-dzm)*pr[j,i]+
                  (dxm)*(1-dzm)*pr[j,i+1]+
                  (1-dxm)*(dzm)*pr[j+1,i]+
                  (dxm)*(dzm)*pr[j+1,i+1])
        }
      }
      # Close connection
      close(f.prn)
      if(mfind == 0){
        break
      }
      # .prn (tstep) counter
      tstep <- tstep+1
    }

    # Save grid
    assign(paste0('grid.', f.name),
      # Combine into one tibble
      purrr::map2(
        list(pr, vx, vz, exx, ezz, exz, sxx,
             szz, sxz, ro, nu, nd, mu, ep,
             et, pr0, prb, dv, tk, cp, kt, ht),
        c('pr', 'vx', 'vz', 'exx', 'ezz', 'exz', 'sxx',
          'szz', 'sxz', 'ro', 'nu', 'nd', 'mu', 'ep',
          'et', 'pr0', 'prb', 'dv', 'tk', 'cp', 'kt', 'ht'),
        ~{rownames(.x) <- gz
          colnames(.x) <- gx
          .x %>%
          as_tibble(rownames = 'z') %>%
          tidyr::pivot_longer(-z, names_to = 'x', values_to = .y) %>%
          mutate('x' = as.integer(x), 'z' = as.integer(z))}) %>%
      purrr::reduce(dplyr::left_join, by = c('z', 'x')))

    # Save
    grids[[g.count]] <- get(paste0('grid.', f.name))
    g.count <- g.count + 1
    cat('\nSaving', paste0('grid.', f.name))
  }

  # Save markers
  assign(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')),
    # Combine into one tibble
    purrr::map2(
      list(mti, mxx, mzz, mtk, mpb, mty),
      c('time', 'x', 'z', 'T', 'P', 'type'),
      ~{.x %>%
        as_tibble(rownames = 'id') %>%
        tidyr::pivot_longer(-id, names_to = 'tstep', values_to = .y) %>%
        mutate('tstep' = tstep %>% stringr::str_extract('[0-9]+') %>% as.integer(),
               'id' = as.integer(id))}) %>%
    purrr::reduce(dplyr::left_join, by = c('id', 'tstep')) %>%
    group_by(id) %>%
    tidyr::drop_na())
  cat('\nSaving markers', paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')), '\n')

  if(markers == TRUE & grid == TRUE) {
    # Print grids
    cat('\nSaved ', length(grids), ' grids:')
    cat('\n', fnames, sep = '\n')
    # Print markers
    print(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+'))))
    cat('\nMarker types')
    print(table(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')))$type))

    return(list(grid = grids %>% purrr::set_names(fnames),
                marx = get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')))))

  } else if(markers == TRUE & grid == FALSE) {
    # Print markers
    print(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+'))))
    cat('\nMarker types')
    print(table(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')))$type))

    return(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+'))))
  } else if(markers == FALSE & grid == TRUE) {
    # Print grids
    cat('\nSaved ', length(grids), ' grids:')
    cat('\n', fnames, sep = '\n')
    return(grids)
  }

}

# Save rock type colormap
tibble(type = seq_len(40), r = c(0.792, 0.50588, 1, 0.68235, 1, 0.75294,
    0.50196, 0, 0, 0, 0.3, 0.14118, 0, 0.9, 0.4, 0.8549, 0.95294, 0.35294, 0.1,
    0, 0, 0, 0, 1, 1, 0.46667, 0.50196, 0.72549, 0.82549, 0.6, 1, 0.99216, 0.84706,
    0.9, 0.8, 1, 0.6, 0.6, 0.1, 0.6), g = c(0.862, 0.99608, 0.50196, 0.34118, 0.50196,
    0.75294, 0.50196, 0.50196, 0.84314, 0, 0.3, 0.72157, 0.50196, 0.4, 0, 0.59608,
    0.20392, 0.16863, 0.6, 0, 0, 0, 0, 1, 0.90196, 0.46667, 0.50196, 0.015686,
    0.43922, 0, 0, 0.38824, 0.078431, 0.2, 0, 0.6, 0.4, 0.8, 0.8, 0.5), b = c(0.988,
    0.78824, 0, 0, 0, 0.75294, 0.50196, 0, 0, 0.71765, 0.9, 0.99216, 1, 1, 0,
    0.36078, 0.086275, 0.027451, 0, 0, 0, 0, 0, 0.31765, 0.18824, 0.23529, 0,
    0.78431, 0.99608, 0, 0, 0.30196, 0.15294, 0.2, 0, 0, 0, 0, 0, 0)) %>%
mutate(color = rgb(r, g, b, 1)) -> c.map

# Load RData files and save markers and grids
load_marx <- function(path) {
  # Model name
  mod <- path %>% stringr::str_extract('cd.[0-9]+')
  cat('\nLoading markers and grids [', mod, ']', sep = '')
  # Load .RData file
  load(path)
  # Take the markers dataframe and ...
  cat('\nCropping and filtering markers [', mod, ']', sep = '')
  marx <- get(mod)$marx %>%
  # Crop tsteps before marker type change
    slice(if(is.na(which(type != dplyr::lag(type))[1])) seq_len(n())
          else seq_len(which(type != dplyr::lag(type))[1])-1) %>%
  # Depth from crust surface
  mutate('z' = z-18) %>%
  # Add rock type color map
  left_join(c.map, by = 'type') %>%
  select(-c(r, g, b))
  # Extract time from marker dataframe
  time <- unique(marx$time)
  # Add time metadata field to grids
  grid <- purrr::map2(get(mod)$grid, time, ~{attr(.x, 'time') <- .y; .x})
  # Save
  cat('\nSaving markers and grids [', mod, ']', sep = '')
  # Save markers and grids separately
  assign(paste0(mod, '.grid'), grid, envir = .GlobalEnv)
  assign(paste0(mod, '.marx'), marx, envir = .GlobalEnv)
}

# Save available features
c(
  'tsteps',
  'under.three.kbar',
  'under.ten.kbar',
  'above.thirty.kbar',
  'under.one.hundred.c',
  'under.five.hundred.c',
  'above.seven.hundred.c',
  'max.P',
  'med.P',
  'mean.P',
  'iqr.P',
  'max.T',
  'med.T',
  'mean.T',
  'iqr.T',
  'up.dx',
  'down.dx',
  'runup.dx',
  'rundown.dx',
  'sum.dx',
  'sumup.dx',
  'sumdown.dx',
  'up.dP',
  'down.dP',
  'runup.dP',
  'rundown.dP',
  'sum.dP',
  'sumup.dP',
  'sumdown.dP',
  'up.dT',
  'down.dT',
  'runup.dT',
  'rundown.dT',
  'sum.dT',
  'sumup.dT',
  'sumdown.dT') -> features

# Computes marker features
marx_ft <- function(df, features = 'all') {
  if(features == 'all') {
    cat('\nComputing features:',
        'tsteps',
        'under.three.kbar',
        'under.ten.kbar',
        'above.thirty.kbar',
        'under.one.hundred.c',
        'under.five.hundred.c',
        'above.seven.hundred.c',
        'max.P',
        'med.P',
        'iqr.P',
        'max.T',
        'med.T',
        'iqr.T',
        'up.dx',
        'down.dx',
        'runup.dx',
        'rundown.dx',
        'sum.dx',
        'sumup.dx',
        'sumdown.dx',
        'up.dP',
        'down.dP',
        'runup.dP',
        'rundown.dP',
        'sum.dP',
        'sumup.dP',
        'sumdown.dP',
        'up.dT',
        'down.dT',
        'runup.dT',
        'rundown.dT',
        'sum.dT',
        'sumup.dT',
        'sumdown.dT',
        sep = '\n')
  } else {
    cat('\nComputing features:', features, sep = '\n')
  }
  # Compute features
  df %>%
  summarise(
    tsteps = n(),
    under.three.kbar = sum(P < 3e3),
    under.ten.kbar = sum(P < 1e4),
    above.thirty.kbar = sum(P > 3e4),
    under.one.hundred.c = sum(T < 373),
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
    up.dx = sum(diff(x) > 0),
    down.dx = sum(diff(x) < 0),
    runup.dx = {rn <- rle(diff(x) > 0); rn$lengths[which(rn$values == TRUE)] %>% max()},
    rundown.dx = {rn <- rle(diff(x) > 0); rn$lengths[which(rn$values == FALSE)] %>% max()},
    sum.dP = sum(diff(P)),
    sumup.dP = sum(diff(P)[which(diff(P) > 0)]),
    sumdown.dP = sum(diff(P)[which(diff(P) < 0)]),
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
  # Save
  if(features == 'all' | length(features) == 0) {
    return(d)
  } else {
    return(d %>% select(features))
  }
}

# Markers motion movie
marx_motion_mov <- function(df, name, class = FALSE) {
  cat('\nAnimating marker motions [', name, ']', sep = '')
  # Color by type or class
  if(class) {
  # Draw plot
  df %>%
  mutate('time' = round(time/1e6, 3),
         'type' = as.factor(type),
         'class' = as.factor(class)) %>%
  ggplot() +
  geom_point(aes(x = x/1000, y = z/1000, color = class, group = id), size = 0.3) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = 'Distance [km]',
    y = 'Depth [km]',
    color = 'Rock Type') +
  scale_y_reverse() +
  coord_fixed() +
  scale_color_brewer(palette = 'Paired') +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_motion_class.mp4')
  cat('\nSaving movie to figs/anim/', name, '_motion_class.mp4', sep = '')
  } else {
  # Draw plot
  df %>%
  mutate('time' = round(time/1e6, 3),
         'type' = as.factor(type)) %>%
  ggplot() +
  geom_point(aes(x = x/1000, y = z/1000, color = type, group = id), size = 0.3) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = 'Distance [km]',
    y = 'Depth [km]',
    color = 'Rock Type') +
  scale_y_reverse() +
  coord_fixed() +
  scale_color_manual(
    values = unique(c.map$color)[as.factor(c.map$type) %>% unique() %>% order()],
    breaks = unique(c.map$type)[as.factor(c.map$type) %>% unique() %>% order()],
    na.value = 'white') +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_motion.mp4')
  cat('\nSaving movie to figs/anim/', name, '_motion.mp4', sep = '')
  }
  # Save video
  cat('\nSaving movie to figs/anim/', name, '_motion.mp4', sep = '')
  animate(
    p,
    fps = 30,
    duration = 20,
    width = 8,
    height = 4,
    units = 'in',
    res = 300,
    renderer = av_renderer(file = fname)
  )
}

# Markers boxplot movie
marx_boxplot_mov <- function(df, name) {
  cat('\nAnimating marker boxplots [', name, ']', sep = '')
  # Draw plot
  df %>%
  mutate(
    'time' = round(time/1e6, 3),
    'x' = x/1000,
    'z' = z/1000,
    'P' = P/1e4,
    'T' = T-273
  ) %>%
  rename(
    'P [kbar]' = P,
    'T [C]' = T,
    'x [km]' = x,
    'z [km]' = z
  ) %>%
  pivot_longer(-c(id, tstep, type, time, color), names_to = 'param') %>%
  group_by(param) %>%
  ggplot() +
  geom_boxplot(aes(x = value, y = param, group = param)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = NULL,
    y = NULL) +
  facet_wrap(~param, scales = 'free') +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  # Save video
  cat('\nSaving movie to figs/anim/', name, '_boxplot.mp4', sep = '')
  animate(
    p,
    fps = 30,
    duration = 20,
    width = 4,
    height = 4,
    units = 'in',
    res = 300,
    renderer = av_renderer(file = paste0('figs/anim/', name, '_boxplot.mp4'))
  )
}

# Markers PT movie
marx_PT_mov <- function(df, name, class = FALSE) {
  cat('\nPlotting marker PT paths [', name, ']', sep = '')
  # Color by type or class
  if(class) {
  # Draw plot
  df %>%
  mutate('time' = round(time/1e6, 3),
         'type' = as.factor(type),
         'class' = as.factor(class)) %>%
  ggplot() +
  geom_point(aes(x = T-273, P/1e4, color = class, group = id), size = 0.3) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = 'Temperature [C]',
    y = 'Pressure [GPa]',
    color = 'Rock Type') +
  coord_cartesian(
    xlim = c(0, 1500),
    ylim = c(0, 4)) +
  scale_color_brewer(palette = 'Paired') +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_PT_class.mp4')
  cat('\nSaving movie to figs/anim/', name, '_PT_class.mp4', sep = '')
  } else {
  # Draw plot
  df %>%
  mutate('time' = round(time/1e6, 3),
         'type' = as.factor(type)) %>%
  ggplot() +
  geom_point(aes(x = T-273, P/1e4, color = type, group = id), size = 0.3) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = 'Temperature [C]',
    y = 'Pressure [GPa]',
    color = 'Rock Type') +
  coord_cartesian(
    xlim = c(0, 1500),
    ylim = c(0, 4)) +
  scale_color_manual(
    values = unique(c.map$color)[as.factor(c.map$type) %>% unique() %>% order()],
    breaks = unique(c.map$type)[as.factor(c.map$type) %>% unique() %>% order()],
    na.value = 'white') +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_PT.mp4')
  cat('\nSaving movie to figs/anim/', name, '_PT.mp4', sep = '')
  }
  # Save video
  animate(
    p,
    fps = 30,
    duration = 20,
    width = 7,
    height = 7,
    units = 'in',
    res = 300,
    renderer = av_renderer(file = fname)
  )
}

# Markers features plot
marx_features_plot <- function(df, name) {
  cat('\nPlotting marker features [', name, ']', sep = '')
  # Draw plot
  df %>%
  pivot_longer(-id, names_to = 'param') %>%
  group_by(param) %>%
  ggplot() +
  geom_histogram(aes(x = value, group = param), bins = 20) +
  facet_wrap(~param, scale = 'free') +
  labs(
    title = paste0('[', name, '] Marker Features'),
    x = NULL,
    y = NULL) +
  theme_classic() -> p
  # Save video
  cat('\nSaving plot to figs/anim/', name, '_features.png', sep = '')
  ggsave(
    plot = p,
    filename = paste0('figs/features/', name, '_features.png'),
    device = 'png',
    type = 'cairo',
    width = 7,
    height = 7,
    units = 'in',
    dpi = 300
  )
}

draw_grid <- function(
  nodes = NULL,
  time,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  arrows = FALSE,
  leg.pos = 'right',
  base.size = 11,
  p.type = c(
    'stress',
    'strain',
    'density',
    'temperature',
    'viscosity',
    'stream'),
  v.pal = 'magma',
  v.direction = 1,
  transparent = TRUE) {
  n <- nodes
  if (p.type == 'temperature') {
    if (!is.null(n)) {
      n %>% ggplot() +
      geom_contour_fill(
        aes(x = x/1000, y = z/1000, z = tk - 273),
        size = 0.1,
        color = NA,
        breaks = c(0, seq(100, 1900, 200))) +
      geom_text_contour(
        aes(x = x/1000, y = z/1000, z = tk - 273),
        stroke = 0.2,
        size = 3,
        breaks = c(0, seq(100, 1900, 200))) +
      labs(
        x = 'km',
        y = 'km',
        fill = bquote(degree * C),
        title = paste0('Temperature  ', time, ' Ma')) +
      coord_equal(expand = F) +
      scale_y_reverse(limits = c(box[2], box[1])) +
      scale_x_continuous(limits = c(box[3], box[4])) +
      scale_fill_viridis_c(option = v.pal, direction = v.direction, na.value = 'transparent') +
      theme_minimal(base_size = base.size) +
      theme(legend.position = leg.pos, axis.text = element_text(color = 'black')) -> p
      if(arrows == TRUE){
        p +
        geom_arrow(
          aes(x = x/1000, y = z/1000, dx = vx, dy = (-vz)),
          skip = 5,
          alpha = 0.3,
          show.legend = F) -> p
      }
    }
  } else if (p.type == 'stress') {
    if (!is.null(n)) {
      n %>% ggplot() +
      geom_contour_fill(
        aes(x = x/1000, y = z/1000, z = log10(sii)),
        size = 0.1,
        color = NA,
        breaks = c(0, seq(100, 1900, 200))) +
      geom_text_contour(
        aes(x = x/1000, y = z/1000, z = tk - 273),
        stroke = 0.2,
        size = 3,
        breaks = c(0, seq(100, 1900, 200))) +
      labs(
        x = 'km',
        y = 'km',
        fill = bquote(log(Pa)),
        title = paste0('Log stress  ', time, ' Ma')) +
      coord_equal(expand = F) +
      scale_y_reverse(limits = c(box[2], box[1])) +
      scale_x_continuous(limits = c(box[3], box[4])) +
      scale_fill_viridis_c(option = v.pal, direction = v.direction, na.value = 'transparent') +
      theme_minimal(base_size = base.size) +
      theme(legend.position = leg.pos, axis.text = element_text(color = 'black')) -> p
      if(arrows == TRUE){
        p +
        geom_arrow(
          aes(x = x/1000, y = z/1000, dx = vx, dy = (-vz)),
          skip = 5,
          alpha = 0.3,
          show.legend = F) -> p
      }
    }
  } else if (p.type == 'strain') {
    if (!is.null(n)) {
      n %>% ggplot() +
      geom_contour_fill(
        aes(x = x/1000, y = z/1000, z = log10(eii)),
        size = 0.1,
        color = NA,
        breaks = c(0, seq(100, 1900, 200))) +
      geom_text_contour(
        aes(x = x/1000, y = z/1000, z = tk - 273),
        stroke = 0.2,
        size = 3,
        breaks = c(0, seq(100, 1900, 200))) +
      labs(
        x = 'km',
        y = 'km',
        fill = bquote(log(s^-1)),
        title = paste0('Log strain rate  ', time, ' Ma')) +
      coord_equal(expand = F) +
      scale_y_reverse(limits = c(box[2], box[1])) +
      scale_x_continuous(limits = c(box[3], box[4])) +
      scale_fill_viridis_c(option = v.pal, direction = v.direction, na.value = 'transparent') +
      theme_minimal(base_size = base.size) +
      theme(legend.position = leg.pos, axis.text = element_text(color = 'black')) -> p
      if(arrows == TRUE){
        p +
        geom_arrow(
          aes(x = x/1000, y = z/1000, dx = vx, dy = (-vz)),
          skip = 5,
          alpha = 0.3,
          show.legend = F) -> p
      }
    }
  } else if (p.type == 'density') {
    if (!is.null(n)) {
      n %>% ggplot() +
      geom_contour_fill(
        aes(x = x/1000, y = z/1000, z = ro),
        size = 0.1,
        color = NA,
        breaks = c(0, seq(100, 1900, 200))) +
      geom_text_contour(
        aes(x = x/1000, y = z/1000, z = tk - 273),
        stroke = 0.2,
        size = 3,
        breaks = c(0, seq(100, 1900, 200))) +
      labs(
        x = 'km',
        y = 'km',
        fill = bquote(kg~m^-3),
        title = paste0('Density  ', time, ' Ma')) +
      coord_equal(expand = F) +
      scale_y_reverse(limits = c(box[2], box[1])) +
      scale_x_continuous(limits = c(box[3], box[4])) +
      scale_fill_viridis_c(option = v.pal, direction = v.direction, na.value = 'transparent') +
      theme_minimal(base_size = base.size) +
      theme(legend.position = leg.pos, axis.text = element_text(color = 'black')) -> p
      if(arrows == TRUE){
        p +
        geom_arrow(
          aes(x = x/1000, y = z/1000, dx = vx, dy = (-vz)),
          skip = 5,
          alpha = 0.3,
          show.legend = F) -> p
      }
    }
  } else if (p.type == 'viscosity') {
    if (!is.null(n)) {
      n %>% ggplot() +
      geom_contour_fill(
        aes(x = x/1000, y = z/1000, z = log10(nu)),
        size = 0.1,
        color = NA,
        breaks = c(0, seq(100, 1900, 200))) +
      geom_text_contour(
        aes(x = x/1000, y = z/1000, z = tk - 273),
        stroke = 0.2,
        size = 3,
        breaks = c(0, seq(100, 1900, 200))) +
      labs(
        x = 'km',
        y = 'km',
        fill = bquote(log(Pa %.% s)),
        title = paste0('Log Viscosity  ', time, ' Ma')) +
      coord_equal(expand = F) +
      scale_y_reverse(limits = c(box[2], box[1])) +
      scale_x_continuous(limits = c(box[3], box[4])) +
      scale_fill_viridis_c(option = v.pal, direction = v.direction, na.value = 'transparent') +
      theme_minimal(base_size = base.size) +
      theme(legend.position = leg.pos, axis.text = element_text(color = 'black')) -> p
      if(arrows == TRUE){
        p +
        geom_arrow(
          aes(x = x/1000, y = z/1000, dx = vx, dy = (-vz)),
          skip = 5,
          alpha = 0.3,
          show.legend = F) -> p
      }
    }
  } else if (p.type == 'stream') {
    if (!is.null(n)) {
      n %>% ggplot() +
        geom_contour_fill(
          aes(x = x / 1000, y = z / 1000, z = tk - 273),
          size = 0.1,
          color = 'white',
          alpha = 0.5,
          breaks = c(0, seq(100, 1900, 200)),
          show.legend = F) +
        geom_streamline(
          aes(
            x = x / 1000,
            y = z / 1000,
            dx = vx,
            dy = (-vz),
            color = sqrt(..dx.. ^ 2 + ..dy.. ^ 2) * 31540000 * 100, alpha = ..step..),
          S = 5,
          dt = 31540000 * 1000 / 5,
          arrow = NULL,
          L = 10,
          res = 1,
          skip = 5,
          lineend = 'round',
          size = 0.3) +
        geom_contour(
          aes(x = x / 1000, y = z / 1000, z = tk - 273),
          size = 0.15,
          color = 'white',
          na.rm = T,
          breaks = c(0, seq(100, 1900, 200))) +
        geom_text_contour(
          aes(x = x / 1000, y = z / 1000, z = tk - 273),
          stroke = 0.2,
          size = 3,
          breaks = c(0, seq(100, 1900, 200))) +
        labs(
          x = 'km',
          y = 'km',
          color = bquote(cm ~ yr ^ -1),
          title = paste0('Streams  ', time, ' Ma')) +
        guides(alpha = F, fill = F) +
        coord_equal(expand = F) +
        scale_y_reverse(limits = c(box[2], box[1])) +
        scale_x_continuous(limits = c(box[3], box[4])) +
        scale_fill_gradient(low = 'grey90', high = 'grey10') +
        scale_color_viridis_c(option = v.pal, direction = v.direction, na.value = 'transparent') +
        theme_minimal(base_size = base.size) +
        theme(
          legend.position = leg.pos,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = 'black')
        ) -> p
    }
  }
  if (transparent == TRUE) {
    p +
    theme(
      plot.background = element_rect(fill = 'transparent', color = NA),
      panel.background = element_rect(fill = 'transparent', color = NA),
      legend.background = element_rect(fill = 'transparent',
      color = NA)
    ) -> p
  }
  return(p)
}

# b.ic ----
# Bayesian Information Criterion
b.ic <- function(
  df,
  features = 'all',
  G = NULL,
  init = NULL,
  scale = FALSE) {
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
    gmmData <- gmmData
  } else {
    gmmData <- gmmData %>%
      tibble() %>%
      scale() %>%
      replace_na(0)
  }
  if(!is.null(init)) {
    BIC <- mclustBIC(gmmData, G = G, initialization = list(hcPairs = hc(gmmData, use = init)))
  } else {
    BIC <- mclustBIC(gmmData, G = G)
  }
  print(summary(BIC))
  return(BIC)
}
# gmm ----
# Gaussian mixing model (GMM) classification
gmm <- function(df, features = 'all', G = NULL, bic = NULL) {
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
    m <- Mclust(gmmData, G = G, x = bic)
  } else {
    gmmData <- gmmData %>%
      tibble() %>%
      scale() %>%
      replace_na(0)
    m <- Mclust(gmmData, G = G, x = bic)
  }
  print(summary(m))
  df$class <- as.factor(m$classification)
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
  return(list(model = m))
}
