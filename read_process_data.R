load('data.coupling.RData')
couplingList <- ls()
source('functions.R')

# Reduce kerswell_et_al_coupling dataset
data.coupling <- allModelsDataframe %>%
  select(model, plate.age, conv.velocity, tparam) %>%
  left_join(reducedAllModelsDataframe) %>%
  distinct() %>%
  tibble()
# Read in elevator data
data.elevator <- vroom('./dataset/elevator.txt')
# Filepaths to data
files <- tibble(
  paths = dir_ls('./dataset', regexp = '*marks.txt'),
  fnames = substring(paths, regexpr("/cd", paths) + 1, regexpr("_marks", paths) - 1)
)
# files <- files %>% slice_sample(n = 4)
# Read raw data
marx <- tibble(model = as.factor(files$fnames),
               filepath = files$paths) %>%
  left_join(data.coupling) %>%
  left_join(data.elevator) %>%
  mutate(data.raw = purrr::map(filepath, ~ data.get(.x)))
names(marx$data.raw) <- files$fnames
# Tidy data and add features for ML classification
marx <- marx %>%
  mutate(data.tidy = data.raw %>% 
           purrr::map(~ data.tidy(.x)),
         mark.ft = data.tidy %>%
           purrr::map(~ add.ft(.x)))
# Clean up environment
rm(list = couplingList)
rm(couplingList, data.coupling, data.elevator, files, bic.col, c.pal)
rm(list = lsf.str())
# Save data ----
save.image(file = 'marx.RData')
