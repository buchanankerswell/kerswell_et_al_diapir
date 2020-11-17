# rm(list = ls())

dirs.list <- list.dirs('/Volumes/My Passport/numerical_models/common_depth_2017_fixed_geotherms/dataset', recursive = F)
dirs <- dirs.list[grepl('*km', dirs.list)]
    

# Close connection
# close(f.prn)