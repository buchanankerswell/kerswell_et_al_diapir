# Load packages
# Quiet loading
sshhh <- function(p){
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, character.only=TRUE)))}

# Package list
c('magrittr', 'tidyr', 'readr', 'purrr',
  'mclust', 'ggforce', 'dplyr', 'ggrepel',
  'patchwork', 'gridExtra', 'gganimate',
  'ggridges', 'progress', 'viridis', 'metR',
  'parallel', 'progress') -> p.list

#cat('Loading libraries:', p.list, sep = '\n')

# auto-load quietly
sapply(p.list, sshhh)

#cat('Loading functions\n')

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
      pr <- matrix(NA, znumz, xnumx) # Pressure [Pa]
      vx <- matrix(NA, znumz, xnumx) # Velocity [m/s]
      vz <- matrix(NA, znumz, xnumx) # Velocity [m/s]
      exx <- matrix(NA, znumz, xnumx) # Strain rate [1/s]
      ezz <- matrix(NA, znumz, xnumx) # Strain rate [1/s]
      exz <- matrix(NA, znumz, xnumx) # Strain rate [1/s]
      sxx <- matrix(NA, znumz, xnumx) # Stress [Pa]
      szz <- matrix(NA, znumz, xnumx) # Stress [Pa]
      sxz <- matrix(NA, znumz, xnumx) # Stress [Pa]
      ro <- matrix(NA, znumz, xnumx) # Density [kg/m^3]
      nu <- matrix(NA, znumz, xnumx) # Viscosity [Pa s]
      nd <- matrix(NA, znumz, xnumx) # ?
      mu <- matrix(NA, znumz, xnumx) # Standard viscosity for node [Pa s]
      ep <- matrix(NA, znumz, xnumx) # Surface trace
      et <- matrix(NA, znumz, xnumx) # Free array
      pr0 <- matrix(NA, znumz, xnumx) # Last cycle pressure [Pa]
      prb <- matrix(NA, znumz, xnumx) # ?
      dv <- matrix(NA, znumz, xnumx) # ?
      tk <- matrix(NA, znumz, xnumx) # Temperature [K]
      cp <- matrix(NA, znumz, xnumx) # Heat capacity [J/kg]
      kt <- matrix(NA, znumz, xnumx) # Thermal conductivity [Wt/m K]
      ht <- matrix(NA, znumz, xnumx) # Heat sources [W/m^3]
      eii <- matrix(1, znumz, xnumx)*1e-16 # Strain rate tensor [1/s]
      sii <- matrix(1, znumz, xnumx)*1e+4 # Stress tensor [Pa]
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
  mutate('z' = z-18000) %>%
  # Add rock type color map
  left_join(c.map, by = 'type') %>%
  select(-c(r, g, b))
  # Extract time from marker dataframe
  time <- unique(marx$time)
  tstep <- unique(marx$tstep)
  # Add time metadata field to grids
  grid <- purrr::map2(get(mod)$grid[tstep], time, ~{attr(.x, 'time') <- .y; .x})
  # Surface trace
  cat('\nDefining surface trace [', mod, ']', sep = '')
  purrr::map2_lgl(grid, time, ~{
    # Reshape surface data
    .x %>%
    select(z, x, ep) %>%
    pivot_wider(names_from = x, values_from = ep) %>%
    select(-1) %>%
    unlist(use.names = F) -> ep
    # Double gx vector by adding intermed interces
    approx(seq_along(unique(.x$x)), unique(.x$x), n = length(unique(.x$x))*2-1)$y -> gx
    # Make datframe
    tibble(x = gx, ep = ep[ep != 0][-1]) -> df
    # Find where trench interferes with convergence region
    if(df$x[which.min(df$ep)] >= 500000) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }) -> interf
  # Add attribute representing the timestep of interference
  if(all(interf == FALSE)) {
    cut <- length(interf)
  } else {
    cut <- min(which(interf == TRUE))
  }
  attr(marx, 'tcut') <- cut
  cat(
    '\nInterference with convergence region at [', mod, ']:',
    '\ntstep:',
    cut,
    '\ntime:',
    round(time[cut]/1e6, 1), 'Ma',
    sep = '')
  # Save
  cat('\nSaving markers and grids [', mod, ']\n', sep = '')
  # Save markers and grids separately
  assign(paste0(mod, '.grid'), grid, envir = .GlobalEnv)
  assign(paste0(mod, '.marx'), marx, envir = .GlobalEnv)
}

# Save available features
c(
  'tsteps',
  'above.three.kbar',
  'above.ten.kbar',
  'above.twenty.kbar',
  'above.thirty.kbar',
  'above.fourty.kbar',
  'above.fifty.kbar',
  'above.one.hundred.c',
  'above.two.hundred.c',
  'above.three.hundred.c',
  'above.four.hundred.c',
  'above.five.hundred.c',
  'above.seven.hundred.c',
  'above.eight.hundred.c',
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
  'sumdown.dT'
) -> features

# Computes marker features
marx_ft <- function(df, features = 'all') {
  if(features == 'all') {
    cat(
      '\nComputing features:',
      'tsteps',
      'above.three.kbar',
      'above.ten.kbar',
      'above.twenty.kbar',
      'above.thirty.kbar',
      'above.fourty.kbar',
      'above.fifty.kbar',
      'above.one.hundred.c',
      'above.two.hundred.c',
      'above.three.hundred.c',
      'above.four.hundred.c',
      'above.five.hundred.c',
      'above.seven.hundred.c',
      'above.eight.hundred.c',
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
      'sumdown.dT',
      sep = '\n'
    )
  } else {
    cat('\nComputing features:', features, sep = '\n')
  }
  # Compute features
  df %>%
  summarise(
    tsteps = n(),
    above.three.kbar = sum(P > 3e3),
    above.ten.kbar = sum(P > 1e4),
    above.twenty.kbar = sum(P > 2e4),
    above.thirty.kbar = sum(P > 3e4),
    above.fourty.kbar = sum(P > 4e4),
    above.one.hundred.c = sum(T > 373),
    above.two.hundred.c = sum(T > 473),
    above.three.hundred.c = sum(T > 573),
    above.four.hundred.c = sum(T > 673),
    above.five.hundred.c = sum(T > 773),
    above.seven.hundred.c = sum(T > 973),
    above.eight.hundred.c = sum(T > 1073),
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
    sum.dx = sum(diff(x)),
    sumup.dx = sum(diff(x)[which(diff(x) > 0)]),
    sumdown.dx = sum(diff(x)[which(diff(x) < 0)]),
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
    return(d %>% select(c(id, all_of(features))))
  }
}

# Gaussian mixture modelling (Scrucca et al., 2016) to classify recovered rocks
marx_classify <- function(marx.df, fts.df, k = 2) {
  # Crop markers > interference timecut
  tcut <- attr(marx.df, 'tcut')
  m <- marx.df %>% slice(1:tcut)
  # Save IDs
  ids <- fts.df$id
  # Fit Eigenvalue decomposition models and select best using BIC
  cat('\nSelecting best Gaussian mixture model by BIC\n')
  fts.df %>% ungroup() %>% select(-id) -> X
  # Try clustering
  mc.bic <- try(mclustBIC(X, G = seq_len(k)))
  # If clustering doesn't converge or throws error
  if(class(mc.bic) == 'try-error') {
    # Keep trying i times
    i <- 1 # Counter
    imax <- 20
    while(class(mc.bic) == 'try-error' && i <= imax) {
      if(i < imax) {
        mc.bic <- try(mclustBIC(X, G = seq_len(k)))
        i <- i + 1
      } else {
        stop('Clustering could not converge')
      }
    }
  }
  #print(summary(mc.bic))
  # GMM clustering using model picked by BIC
  mc <- Mclust(X, x = mc.bic, verbose = T)
  #print(summary(mc))
  # Make class id df
  class.df <- tibble(id = ids, class = mc$classification)
  # Join classification to markers data
  marx.class.df <- m %>% left_join(class.df, by = 'id') %>% replace(is.na(.), 0)
  # Max pressure for all markers
  mft <- marx.class.df %>%
  filter(class != 0) %>%
  summarise(maxP = max(P), sumdP = sum(diff(P)))
  # Get parameter centroids
  centroids <- t(mc$parameters$mean) %>% as_tibble() %>% mutate(class = 1:n(), .before = 'sum.dP')
  # Calculate lower bounds for features
  lowerB <- X %>% summarise(across(everything(), ~{median(.x) - 3/4*IQR(.x)}))
  # Classify group as recovered if maxP or sumdP centroid is below lower bound
  d <- centroids[centroids$sum.dP <= lowerB$sum.dP | centroids$max.P <= lowerB$max.P,]
  # Add recovered class
  recovered.df <- class.df %>%
  mutate(recovered = ifelse(class %in% d$class, TRUE, FALSE))
  # Join with marker data
  df <- marx.class.df %>%
  left_join(recovered.df, by = c('id', 'class')) %>%
  replace(is.na(.), FALSE)
  # Summarise maximum pressures
  #mP.class <- df %>% group_by(id, class) %>% summarise(maxP = max(P), sumdP = sum(diff(P)), .groups = 'keep')
  #mP <- df %>% group_by(id, recovered) %>% summarise(maxP = max(P), .groups = 'keep')
  # Change recovered to FALSE if recovered was TRUE and marker maxP was >= meadian + IQR
  #misclassed.marx <- which(mP$maxP >= median(mP$maxP) & mP$recovered == TRUE)
  #df$recovered[df$id %in% misclassed.marx] <- FALSE
  # Print results
  cat('\nRecovered classes:')
  df %>% slice(1) %>% group_by(class) %>% select(class, recovered) %>% table() %>% print()
  # Join classification (subducted or recovered) to markers data
  return(
    list(
      'marx' = df,
      'mc' = mc
    )
  )
}

# Calculate ratios, statistics, and P-T-x-z of markers
marx_stats <- function(df) {
  # Stats
  cat('\nCalculating statistics')
  df %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  summarise(
    recovered = sum(recovered),
    subducted = n()-recovered,
    ratio = recovered/subducted,
    perc = recovered/subducted * 100,
    max.P.rec = max(df$P[which(df$recovered == TRUE)]),
    med.P.rec = median(df$P[which(df$recovered == TRUE)]),
    iqr.P.rec = IQR(df$P[which(df$recovered == TRUE)]),
    mean.P.rec = mean(df$P[which(df$recovered == TRUE)]),
    sd.P.rec = sd(df$P[which(df$recovered == TRUE)]),
    min.P.rec = min(df$P[which(df$recovered == TRUE)]),
    max.T.rec = max(df$T[which(df$recovered == TRUE)]),
    med.T.rec = median(df$T[which(df$recovered == TRUE)]),
    iqr.T.rec = IQR(df$T[which(df$recovered == TRUE)]),
    mean.T.rec = mean(df$T[which(df$recovered == TRUE)]),
    sd.T.rec = sd(df$T[which(df$recovered == TRUE)]),
    min.T.rec = min(df$T[which(df$recovered == TRUE)]),
    max.P.sub = max(df$P[which(df$recovered == FALSE)]),
    med.P.sub = median(df$P[which(df$recovered == FALSE)]),
    iqr.P.sub = IQR(df$P[which(df$recovered == FALSE)]),
    mean.P.sub = mean(df$P[which(df$recovered == FALSE)]),
    sd.P.sub = sd(df$P[which(df$recovered == FALSE)]),
    min.P.sub = min(df$P[which(df$recovered == FALSE)]),
    max.T.sub = max(df$T[which(df$recovered == FALSE)]),
    med.T.sub = median(df$T[which(df$recovered == FALSE)]),
    iqr.T.sub = IQR(df$T[which(df$recovered == FALSE)]),
    mean.T.sub = mean(df$T[which(df$recovered == FALSE)]),
    sd.T.sub = sd(df$T[which(df$recovered == FALSE)]),
    min.T.sub = min(df$T[which(df$recovered == FALSE)])
  ) -> s
  cat('\nCalculating empirical probability distribution')
  cdf.P <- df %>%
  filter(recovered == TRUE) %>%
  summarise(maxP = max(P)) %>%
  arrange(maxP) %>%
  mutate(cdf = (row_number()-1)/n())
  cdf.T <- df %>%
  filter(recovered == TRUE) %>%
  summarise(maxT = max(T)) %>%
  arrange(maxT) %>%
  mutate(cdf = (row_number()-1)/n())
  return(list('stats' = s, 'cdfT' = cdf.T, 'cdfP' = cdf.P))
}

# Monte carlo sampling of marx_classify()
monte_carlo <- function(marx.df, fts.df, n = 100, k = 2) {
  # Progress bar
  pb <- progress_bar$new(
    format = 'Monte Carlo Sampling [:bar] :current/:total (:percent)',
    total = n, clear = FALSE, width = 80)
  purrr::map(seq_len(n), ~{
    pb$tick()
    marx_classify(marx.df, fts.df, k = k)$marx %>%
    marx_stats()
  })
}

# Markers motion movie
marx_motion_mov <- function(df, name, class = FALSE, recovered = FALSE) {
  cat('\nAnimating marker motions [', name, ']', sep = '')
  # Color by type, class, or recovered
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
    color = 'Classification') +
  scale_y_reverse() +
  coord_fixed() +
  scale_color_brewer(palette = 'Paired', breaks = as.character(seq_len(max(df$class)))) +
  theme_minimal() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_motion_class.mp4')
  } else if(recovered) {
  # Draw plot
  df %>%
  mutate('time' = round(time/1e6, 3),
         'type' = as.factor(type),
         'class' = as.factor(class)) %>%
  ggplot() +
  geom_point(aes(x = x/1000, y = z/1000, color = recovered, group = id), size = 0.3) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = 'Distance [km]',
    y = 'Depth [km]',
    color = 'Recovered') +
  scale_y_reverse() +
  coord_fixed() +
  scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'grey50')) +
  theme_minimal() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_motion_recovered.mp4')
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
  theme_minimal() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_motion.mp4')
  }
  # Save video
  cat('\nSaving movie to ', fname, sep = '')
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
marx_PT_mov <- function(df, name, class = FALSE, recovered = FALSE) {
  cat('\nAnimating marker PT paths [', name, ']', sep = '')
  # Color by type, class, or recovered
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
    color = 'Classification') +
  coord_cartesian(
    xlim = c(0, 1500),
    ylim = c(0, 4)) +
  scale_color_brewer(palette = 'Paired', breaks = as.character(seq_len(max(df$class)))) +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_PT_class.mp4')
  } else if(recovered) {
  # Draw plot
  df %>%
  mutate('time' = round(time/1e6, 3),
         'type' = as.factor(type),
         'class' = as.factor(class)) %>%
  ggplot() +
  geom_point(aes(x = T-273, P/1e4, color = recovered, group = id), size = 0.3) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  labs(
    title = paste0('[', name, '] Time: {round(frame_time, 2)} Ma'),
    x = 'Temperature [C]',
    y = 'Pressure [GPa]',
    color = 'Recovered') +
  coord_cartesian(
    xlim = c(0, 1500),
    ylim = c(0, 4)) +
  scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'grey50')) +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  ) +
  exit_fade() +
  transition_time(time) +
  ease_aes('linear') -> p
  fname <- paste0('figs/anim/', name, '_PT_recovered.mp4')
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
  }
  # Save video
  cat('\nSaving movie to ', fname, sep = '')
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

# Draw grids
draw_grid <- function(
  node,
  model,
  marx = NULL,
  class = c('type', 'class', 'recovered', 'none'),
  time,
  box = c(up = -18, down = 200, left = 0, right = 2000),
  arrows = FALSE,
  leg.pos = 'right',
  leg.dir = 'vertical',
  leg.dir.rec = 'vertical',
  base.size = 11,
  p.type = c('blank', 'temperature', 'viscosity', 'stream'),
  bk.alpha = 0.6,
  bk.col = rgb(1, 1, 1, 1),
  mk.alpha = 1,
  mk.size = 0.3,
  iso.alpha = 1,
  iso.col = 'white',
  iso.skip = 0,
  stm.alpha = 1,
  rec.col = 'white',
  sub.col = 'black',
  v.pal = 'magma',
  v.direction = 1,
  transparent = TRUE) {
  # Get time cutoff
  tcut <- time
  if(p.type == 'blank') {
    node %>%
    ggplot() +
    geom_contour(
      aes(x = x/1000, y = z/1000, z = tk - 273),
      color = 'black',
      breaks = c(0, seq(100, 1900, 200)),
      size = 0.3,
      alpha = iso.alpha) +
    geom_text_contour(
      aes(x = x/1000, y = z/1000, z = tk - 273),
      stroke = 0.2,
      size = 3,
      skip = iso.skip,
      breaks = c(0, seq(100, 1900, 200))) +
    labs(
      title = paste0('Markers [', model, '] ', tcut, ' Ma'),
      x = 'Distance [km]',
      y = 'Depth [km]') +
    guides(fill = guide_colorbar(direction = leg.dir, title.vjust = 1)) +
    scale_fill_viridis(option = v.pal, direction = v.direction) +
    scale_y_reverse() +
    coord_equal(
      expand = F,
      ylim = c(box[2], box[1]),
      xlim = c(box[3], box[4])) +
    theme_minimal(base_size = base.size) +
    theme(
      panel.grid = element_blank(),
      legend.position = leg.pos
    ) -> p
  } else if(p.type =='temperature') {
    node %>%
    ggplot() +
    geom_contour_fill(aes(x = x/1000, y = z/1000, z = tk - 273), alpha = bk.alpha) +
    geom_contour(
      aes(x = x/1000, y = z/1000, z = tk - 273),
      color = iso.col,
      breaks = c(0, seq(100, 1900, 200)),
      size = 0.3,
      alpha = iso.alpha) +
    geom_text_contour(
      aes(x = x/1000, y = z/1000, z = tk - 273),
      stroke = 0.2,
      size = 3,
      skip = iso.skip,
      breaks = c(0, seq(100, 1900, 200))) +
    labs(
      title = paste0('Temperature [', model, '] ', tcut, ' Ma'),
      x = 'Distance [km]',
      y = 'Depth [km]',
      fill = bquote(degree*C)) +
    guides(fill = guide_colorbar(direction = leg.dir, title.vjust = 1)) +
    scale_fill_viridis(option = v.pal, direction = v.direction) +
    scale_y_reverse() +
    coord_equal(
      expand = F,
      ylim = c(box[2], box[1]),
      xlim = c(box[3], box[4])) +
    theme_minimal(base_size = base.size) +
    theme(
      panel.grid = element_blank(),
      legend.position = leg.pos
    ) -> p
  } else if(p.type =='viscosity') {
    node %>%
    ggplot() +
    geom_contour_fill(aes(x = x/1000, y = z/1000, z = log10(nu)), alpha = bk.alpha) +
    geom_contour(
      aes(x = x/1000, y = z/1000, z = tk - 273),
      color = iso.col,
      breaks = c(0, seq(100, 1900, 200)),
      size = 0.3,
      alpha = iso.alpha) +
    geom_text_contour(
      aes(x = x/1000, y = z/1000, z = tk - 273),
      stroke = 0.2,
      size = 3,
      skip = iso.skip,
      breaks = c(0, seq(100, 1900, 200))) +
    labs(
      title = paste0('Log Viscosity [', model, '] ', tcut, ' Ma'),
      x = 'Distance [km]',
      y = 'Depth [km]',
      fill = 'Log[Pa s]') +
    guides(fill = guide_colorbar(direction = leg.dir, title.vjust = 1)) +
    scale_fill_viridis(option = v.pal, direction = v.direction) +
    scale_y_reverse() +
    coord_equal(
      expand = F,
      ylim = c(box[2], box[1]),
      xlim = c(box[3], box[4])) +
    theme_minimal(base_size = base.size) +
    theme(
      panel.grid = element_blank(),
      legend.position = leg.pos
    ) -> p
  } else if(p.type =='stream') {
    node %>%
    ggplot() +
    geom_streamline(
      aes(x = x / 1000, y = z / 1000, dx = vx, dy = (-vz),
          color = sqrt(..dx.. ^ 2 + ..dy.. ^ 2) * 31540000 * 100),
      alpha = stm.alpha,
      S = 5,
      dt = 31540000 * 1000 / 5,
      arrow = NULL,
      L = 10,
      res = 1,
      skip = 5,
      lineend = "round",
      size = 0.3) +
    geom_contour(
      aes(x = x / 1000, y = z / 1000, z = tk - 273),
      size = 0.3,
      color = iso.col,
      na.rm = T,
      breaks = c(0, seq(100, 1900, 200)),
      alpha = iso.alpha) +
    geom_text_contour(
      aes(x = x / 1000, y = z / 1000, z = tk - 273),
      stroke = 0.2,
      size = 3,
      skip = iso.skip,
      breaks = c(0, seq(100, 1900, 200))) +
    labs(
      x = 'Distance [km]',
      y = 'Depth [km]',
      fill = bquote(degree*C),
      color = bquote('['*cm~yr^-1*']'),
      title = paste0('Streams  [', model, '] ', tcut, ' Ma')) +
    guides(alpha = F, color = guide_colorbar(direction = leg.dir, title.vjust = 1), fill = F) +
    scale_y_reverse() +
    coord_equal(
      expand = F,
      ylim = c(box[2], box[1]),
      xlim = c(box[3], box[4])) +
    scale_color_viridis_c(option = v.pal, limits = c(0, 10), na.value = 'transparent') +
    theme_minimal(base_size = base.size) +
    theme(
      legend.position = leg.pos,
      panel.grid = element_blank(),
      panel.background = element_rect(color = NA, fill = bk.col)
    ) -> p
  }
  if(!is.null(marx)) {
    if(p.type != 'stream') {
      if(class == 'type') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000, color = as.factor(type)),
          shape = 15,
          alpha = mk.alpha,
          size = mk.size) +
        scale_color_manual(
          values = unique(c.map$color)[as.factor(c.map$type) %>% unique() %>% order()],
          breaks = unique(c.map$type)[as.factor(c.map$type) %>% unique() %>% order()],
          na.value = 'white') +
        labs(color = 'Rock type') +
        guides(color = guide_legend(direction = leg.dir.rec, override.aes = list(alpha = 1, size = 2)))
      } else if(class == 'class') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000, color = as.factor(class)),
          alpha = mk.alpha,
          shape = 15,
          size = mk.size) +
        labs(color = 'Class') +
        guides(color = guide_legend(direction = leg.dir.rec, nrow = 4, override.aes = list(alpha = 1, size = 2)))
      } else if(class == 'recovered') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000, color = as.factor(recovered)),
          shape = 15,
          alpha = mk.alpha,
          size = mk.size) +
        labs(color = 'Recovered') +
        scale_color_manual(values = c('TRUE' = rec.col, 'FALSE' = sub.col)) +
        guides(color = guide_legend(direction = leg.dir.rec, override.aes = list(alpha = 1, size = 2))) +
        theme(legend.key = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.5), color = NA))
      } else if(class == 'none') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000),
          shape = 15,
          alpha = mk.alpha,
          size = mk.size) +
        guides(color = guide_legend(direction = leg.dir.rec, nrow = 4, override.aes = list(alpha = 1, size = 2)))
      }
    } else {
      if(class == 'type') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000, shape = as.factor(type)),
          shape = 15,
          alpha = mk.alpha,
          size = mk.size) +
        labs(shape = 'Rock type') +
        guides(shape = guide_legend(direction = leg.dir.rec, override.aes = list(alpha = 1, size = 2)))
      } else if(class == 'class') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          alpha = mk.alpha,
          aes(x = x/1000, y = z/1000, shape = as.factor(class)),
          shape = 15,
          size = mk.size) +
        labs(shape = 'Class') +
        guides(shape = guide_legend(direction = leg.dir.rec, nrow = 4, override.aes = list(alpha = 1, size = 2)))
      } else if(class == 'recovered') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000, fill = as.factor(recovered)),
          shape = 22,
          color = 'transparent',
          size = mk.size) +
        labs(fill = 'Recovered') +
        scale_fill_manual(values = c('TRUE' = rec.col, 'FALSE' = sub.col)) +
        guides(fill = guide_legend(direction = leg.dir.rec, override.aes = list(alpha = 1, size = 2, color = NA))) +
        theme(legend.key = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.5), color = NA))
      } else if(class == 'none') {
        p <- p +
        geom_point(
          data = marx %>% filter(tstep == tcut),
          aes(x = x/1000, y = z/1000),
          alpha = mk.alpha,
          shape = 15,
          size = mk.size) +
        guides(color = guide_legend(direction = leg.dir.rec, nrow = 4, override.aes = list(alpha = 1, size = 2)))
      }
    }
  }
  if(arrows == TRUE){
    p <- p +
    geom_arrow(
      aes(x = x/1000, y = z/1000, dx = vx, dy = (-vz)),
      skip = 5,
      alpha = 0.3,
      show.legend = F)
  }
  if (transparent == TRUE) {
    p <- p +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent",
      color = NA))
  }
  return(p)
}
