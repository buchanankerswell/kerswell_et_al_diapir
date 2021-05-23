source('functions.R')

load_marx('data/cdo46_marx.RData')

grid <- cdo46.grid[[10]]

grid %>% ggplot() +
geom_contour_fill(
  aes(x = x/1000, y = z/1000, z = tk - 273),
  size = 0.1,
  color = NA, breaks = c(0, seq(100, 1900, 200)))

gep <- grid[which(grid$ep != 0),] %>% select(z, x, ep)

grid %>%
mutate(n = row_number()) %>%
ggplot() +
geom_point(aes(x = x, y = z, color = tk)) +
scale_y_reverse() +
coord_equal()

grid$x[which(grid$x %%10 != 0)]

grid %>%
select(z, x, ep) %>%
summary()
ggplot() +
geom_point(aes(x = x, y = ep)) +
coord_cartesian(ylim = c(0, 100))

tibble(
  time = seq_len(length(cdo46.grid)),
  trench = cdo46.grid %>% purrr::map_dbl(~ .x$x[which.min(.x$et)])
) %>%
ggplot() +
geom_point(aes(x = time, y = trench))
