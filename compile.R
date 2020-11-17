load('data/marx.RData')
source('functions.R')

mod.names <- ls()[grepl('marx.*cd', ls())]
marx.lst <- mod.names %>% purrr::map(get) %>% set_names(mod.names)

load('data/coupling.RData')
d.coupling <- reducedAllModelsDataframe %>%
  as_tibble() %>%
  group_by(z1100)

marx <- tibble(
  model = mod.names,
  marx = marx.lst,
  marx.ft = purrr::map(marx.lst, marx_ft)
)

# save(marx, file = 'data/comp.RData')