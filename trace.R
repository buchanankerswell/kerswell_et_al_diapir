source('functions.R')

# File paths (prn binaries)
prns <- list.files('../I2VIS', pattern = '.prn', full.names = T)
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
    marx.est = 20000,
    area = c(500000, 1260000, 17500, 28500),
    markers = T,
    grid = T)
}

models <- unique(mods)

parallel::mclapply(models, fun, mc.cores = 32) %>%
purrr::set_names(models) -> marx

# Save
purrr::walk2(marx, models, ~{
  assign(.y, .x)
  save(list = .y, file = paste0('data/', .y, '_marx.RData'))
})
