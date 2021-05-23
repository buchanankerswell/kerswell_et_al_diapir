# Load functions and libraries
source('functions.R')
# Load marker features
load('data/marx_features.RData')
load_marx('data/cdf78_marx.RData')
load_marx('data/cdd94_marx.RData')
load_marx('data/cdi62_marx.RData')
load_marx('data/cdo46_marx.RData')

j <- classify_recovered(cdd94.marx, marx.features$cdd94)
j <- classify_recovered(cdf78.marx, marx.features$cdf78)
j <- classify_recovered(cdi62.marx, marx.features$cdi62)
j <- classify_recovered(cdo46.marx, marx.features$cdo46)
j$marx$class[which(j$marx$recovered)[1]]
plot(j$mc, what = 'BIC')
plot(j$mc, what = 'classification')
plot(j$mc, what = 'uncertainty')

j$marx %>%
filter(tstep %in% seq(78, 79, 1)) %>%
ggplot() +
geom_point(aes(x = x, y = z, color = as.factor(class))) +
coord_fixed() +
scale_y_reverse()

j$marx %>%
filter(tstep %in% seq(78, 79, 1)) %>%
ggplot() +
geom_point(aes(x = T, y = P, color = as.factor(class)))

# Animate
marx_PT_mov(j$marx, 'cdf78', class = T)
marx_motion_mov(j$marx, 'cdf78', class = T)
