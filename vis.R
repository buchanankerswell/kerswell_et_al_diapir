load('data/comp.RData')

marx
marx$marx$marx.cdf78 %>% filter(tstep >= 10 & tstep <= 30) %>% ggplot() + geom_point(aes(x = x, y = z, color = as.factor(type), group = tstep), size = 0.01) + scale_y_reverse() + facet_wrap(~tstep)
