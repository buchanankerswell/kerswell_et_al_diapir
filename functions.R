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
  'ggridges', 'gifski', 'progress',
  'parallel') -> p.list

cat('Loading libraries:', p.list, sep = '\n')

# auto-load quietly
sapply(p.list, sshhh)

cat('Loading functions\n')

# Read binary (.prn) files and trace markers
trace_marx <- function(
  prn.paths,
  marx.est = 2000000,
  area = c(500000, 1260000, 17500, 28500),
  markers = TRUE,
  grid = FALSE){
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
  fnames <- fpaths %>% stringr::str_extract('cd[a-z][0-9]+_[0-9]+')

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

  for(path in fpaths){
    # Filename
    f <- path
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
          stringr::str_extract(f, 'cd.[0-9]{2}.[0-9]+'),
          '] [:bar] :percent in: :elapsed'),
        total = xnumx*znumz,
        clear = FALSE,
        width = 60)
      # Read nodes information
      for(i in seq_len(xnumx)){
        for(j in seq_len(znumz)){
          vbuf <- readBin(f.prn, 'numeric', 1, 4)
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
          stringr::str_extract(f, 'cd.[0-9]{2}.[0-9]+'),
          '] [:bar] :percent in: :elapsed'),
        total = xnumx*znumz,
        clear = FALSE,
        width = 60)
      # Read nodes information
      for(i in seq_len(xnumx)){
        for(j in seq_len(znumz)){
          vbuf <- readBin(f.prn, 'numeric', 1, 4)
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
        format = 'Calc. Stress & Strain [:bar] :percent in: :elapsed',
        total = (xnumx*znumz)-2440,
        clear = FALSE,
        width = 60)
      for (i in seq_len(xnumx-2)){
        for (j in seq_len(znumz-2)){
          eii[j+1,i+1]=(exz[j+1,i+1]^2+((exx[j+1,i+1]+exx[j+2,i+1]+exx[j+1,i+2]+exx[j+2,i+2])/4)^2)^0.5;
          sii[j+1,i+1]=(sxz[j+1,i+1]^2+((sxx[j+1,i+1]+sxx[j+2,i+1]+sxx[j+1,i+2]+sxx[j+2,i+2])/4)^2)^0.5;
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
      pb.marx.find <- progress_bar$new(
        format = 'Reading Markers [:bar] :percent in: :elapsed',
        total = marknum,
        clear = FALSE,
        width = 60)
      pb.marx <- progress_bar$new(
        format = 'Update Marker PT [:bar] :percent in: :elapsed',
        total = mfind,
        clear = FALSE,
        width = 60)
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
            mpb[mfind,tstep] <- 1e-5*((1-dxm)*(1-dzm)*pr[j,i]+(dxm)*(1-dzm)*pr[j,i+1]+(1-dxm)*(dzm)*pr[j+1,i]+(dxm)*(dzm)*pr[j+1,i+1])
          }
          pb.marx.find$tick()
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
          mpb[mi,tstep] <- 1e-5*((1-dxm)*(1-dzm)*pr[j,i]+(dxm)*(1-dzm)*pr[j,i+1]+(1-dxm)*(dzm)*pr[j+1,i]+(dxm)*(dzm)*pr[j+1,i+1])
        }
        # .prn (tstep) counter
        tstep <- tstep+1
        pb.marx$tick()
      }
      # Close connection
      close(f.prn)
      if(mfind == 0){
        break
      }
    }
    # save grid
    rownames(pr) <- gz; colnames(pr) <- gx
    rownames(vx) <- gz; colnames(vx) <- gx
    rownames(vz) <- gz; colnames(vz) <- gx
    rownames(exx) <- gz; colnames(exx) <- gx
    rownames(ezz) <- gz; colnames(ezz) <- gx
    rownames(exz) <- gz; colnames(exz) <- gx
    rownames(sxx) <- gz; colnames(sxx) <- gx
    rownames(sxz) <- gz; colnames(sxz) <- gx
    rownames(ro) <- gz; colnames(ro) <- gx
    rownames(nu) <- gz; colnames(nu) <- gx
    rownames(nd) <- gz; colnames(nd) <- gx
    rownames(mu) <- gz; colnames(mu) <- gx
    rownames(ep) <- gz; colnames(ep) <- gx
    rownames(et) <- gz; colnames(et) <- gx
    rownames(pr0) <- gz; colnames(pr0) <- gx
    rownames(prb) <- gz; colnames(prb) <- gx
    rownames(dv) <- gz; colnames(dv) <- gx
    rownames(tk) <- gz; colnames(tk) <- gx
    rownames(cp) <- gz; colnames(cp) <- gx
    rownames(kt) <- gz; colnames(kt) <- gx
    rownames(ht) <- gz; colnames(ht) <- gx
    assign(paste0('grid.', stringr::str_extract(f, 'cd.[0-9]+.[0-9]+')), list(
      pr = pr %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'pr') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      vx = vx %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'vx') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      vz = vz %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'vz') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      exx = exx %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'exx') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      ezz = ezz %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'ezz') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      exz = exz %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'exz') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      sxx = sxx %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'sxx') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      szz = szz %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'szz') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      sxz = sxz %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'sxz') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      ro = ro %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'ro') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      nu = nu %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'nu') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      nd = nd %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'nd') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      mu = mu %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'mu') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      ep = ep %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'ep') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      et = et %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'et') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      pr0 = pr0 %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'pr0') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      prb = prb %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'prb') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      dv = dv %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'dv') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      tk = tk %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'tk') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      cp = cp %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'cp') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      kt = kt %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'kt') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer()),
      ht = ht %>% as_tibble(rownames = 'z', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-z, names_to = 'x', values_to = 'ht') %>% mutate('x' = x %>% as.integer(), 'z' = z %>% as.integer())
    ) %>% purrr::reduce(dplyr::left_join))
    grids[[g.count]] <- get(paste0('grid.', stringr::str_extract(f, 'cd.[0-9]+.[0-9]+')))
    g.count <- g.count + 1
  }
  # save markers
  assign(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')), list(
    m.ti = mti %>% as_tibble(rownames = 'id', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-id, names_to = 'tstep', values_to = 'time') %>% mutate('tstep' = tstep %>% substr(4,5) %>% as.integer()),
    m.xx = mxx %>% as_tibble(rownames = 'id', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-id, names_to = 'tstep', values_to = 'x') %>% mutate('tstep' = tstep %>% substr(4,5) %>% as.integer()),
    m.zz = mzz %>% as_tibble(rownames = 'id', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-id, names_to = 'tstep', values_to = 'z') %>% mutate('tstep' = tstep %>% substr(4,5) %>% as.integer()),
    m.tk = mtk %>% as_tibble(rownames = 'id', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-id, names_to = 'tstep', values_to = 'T') %>% mutate('tstep' = tstep %>% substr(4,5) %>% as.integer()),
    m.pb = mpb %>% as_tibble(rownames = 'id', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-id, names_to = 'tstep', values_to = 'P') %>% mutate('tstep' = tstep %>% substr(4,5) %>% as.integer()),
    m.ty = mty %>% as_tibble(rownames = 'id', .name_repair = ~ vctrs::vec_as_names(..., repair = 'unique', quiet = T)) %>% tidyr::pivot_longer(-id, names_to = 'tstep', values_to = 'type') %>% mutate('tstep' = tstep %>% substr(4,5) %>% as.integer())
  ) %>% purrr::reduce(left_join, by = c('id', 'tstep')) %>% group_by(id) %>% tidyr::drop_na())
  if(markers == TRUE & grid == TRUE) {
    # Print markers
    print(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+'))))
    print(table(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')))$type))
    return(list(grid = grids %>% purrr::set_names(fnames),
                marx = get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')))))
  } else if(markers == TRUE & grid == FALSE) {
    # Print markers
    print(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+'))))
    print(table(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+')))$type))
    return(get(paste0('marx.', f %>% stringr::str_extract('cd.[1-9]+'))))
  } else if(markers == FALSE & grid == TRUE) {
    return(grids)
  }
}

# Computes features (columns) for ML clustering
marx_ft <- function(df) {
  d <- df %>%
   add_count(tstep) %>%
   summarise(
     max.t = max(time),
     # max.P = max(P),
     # med.P = median(P),
     # iqr.P = IQR(P),
     # up.dP = sum(diff(P) > 0),
     # down.dP = sum(diff(P) < 0),
     runup.dP = {rn <- rle(diff(P) > 0); rn$lengths[which(rn$values == TRUE)] %>% max()},
     # rundown.dP = {rn <- rle(diff(P) > 0); rn$lengths[which(rn$values == FALSE)] %>% max()},
     sum.dP = sum(diff(P)),
     # sumup.dP = sum(diff(P)[which(diff(P) > 0)]),
     # sumdown.dP = sum(diff(P)[which(diff(P) < 0)]),
     # max.T = max(T),
     # med.T = median(T),
     # iqr.T = IQR(T),
     # up.dT = sum(diff(T) > 0),
     # down.dT = sum(diff(T) < 0),
     runup.dT = {rn <- rle(diff(T) > 0); rn$lengths[which(rn$values == TRUE)] %>% max()},
     # rundown.dT = {rn <- rle(diff(T) > 0); rn$lengths[which(rn$values == FALSE)] %>% max()},
     sum.dT = sum(diff(T)),
     # sumup.dT = sum(diff(T)[which(diff(T) > 0)]),
     # sumdown.dT = sum(diff(T)[which(diff(T) < 0)]),
     .groups = 'keep'
   )
  return(d)
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