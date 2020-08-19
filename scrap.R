load('marx.RData')
source('functions.R')
# Testing new workflow functions
# Remaining challenges --------------------------------------------
#  1. Visualize 1D (2D?) densities for summary variables
# 2. Sample unconstrained models and vis clusters (1D? 2D?) and BIC
# 3. Try dimension reduction (DR)
# 4. Fine-tune EM initialization
# General GMM model

c <- purrr::map_dfr(
  mods,
  .id = 'model',
  .f = function(mod) {
    if (mod.type == 'dr') {
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
