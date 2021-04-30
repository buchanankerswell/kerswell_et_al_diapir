source('functions.R')

# Trace markers
prns <- list.files('data/I2VIS', pattern = '.prn', full.names = T)
mods <- prns %>% stringr::str_extract('cd[a-z]{1}[0-9]{2,3}')

# Compile filepaths
tibble(
  model = mods,
  path = prns) %>%
group_by(model) %>%
mutate(tstep = row_number()-1) -> files

# Trace markers
fun <- function(model) {
  fpaths <- files[files$model == model,]$path
  trace_marx(
    prn.paths = fpaths,
    marx.est = 500000,
    area = c(500000, 1260000, 17500, 28500),
    markers = T,
    grid = T)
}

models <- unique(mods)

marx <- parallel::mclapply(models, fun)

# Save
save(marx, file = 'data/marx.RData')