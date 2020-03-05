#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

## Purpose:
# read the .csv files generated for all scenarios
# report the following:
#    - nb and % of grid cell where sign. detected
#    - for all models and gcm combinations



## Libraries
#library("ncdf4")
library("readr")
suppressMessages(library("dplyr"))
library("ggplot2")

## read files
for (m in c("cwatm_", "matsiro_", "clm45_", "lpjml_", "pcr-globwb_", "watergap2_", "h08_")) { # add more as they come
  for (gcms in c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_")) {
    if ((m == "clm45_" & gcms %in% c("ipsl-cm5a-lr_", "miroc5_")) | (m == "h08_" & gcms == "miroc5_")) {
      next
    }
    if (gcms == "hadgem2-es_") {
      period1 <- "2012-2041"
      period2 <- "2147-2176"
    } else if (gcms == "ipsl-cm5a-lr_") {
      period1 <- "2015-2044"
      period2 <- "2150-2179"
    } else if (gcms == "miroc5_") {
      period1 <- "2035-2064"
      period2 <- "2243-2272"
    }
    for (id in c("max", "min", "p05", "p95")) {
      temp <- read_csv(paste0("/data01/julien/projects/extreme_trans_stab/OUT/", m, gcms, "rcp26_",
                              period1, "_rcp26_", period2, "_", id, ".csv"), col_types = cols(to_save = col_double()))
      temp$id     <- seq(1,259200)
      temp$indice <- id
      temp$model  <- m
      temp$gcm    <- gcms

      if (gcms == "hadgem2-es_" & m == "cwatm_" & id == "max") {
        agg_data <- temp
      } else {
        agg_data <- rbind(agg_data, temp)
      }
    }
  }
}

rm(temp)


#### aggregate
non_agg <- 
agg_data %>% group_by(model, gcm, indice) %>% summarise(n = n(),
                                                        n_pos = sum(to_save >= 0.05, na.rm = TRUE),
                                                        n_na = sum(is.na(to_save)))
non_agg$n_neg    <- (non_agg$n - non_agg$n_na) - non_agg$n_pos
non_agg$cen_neg  <- non_agg$n_neg / (non_agg$n - non_agg$n_na) * 100
non_agg$cen_pos  <- non_agg$n_pos / (non_agg$n - non_agg$n_na) * 100
non_agg

write_csv(non_agg, "./../OUT/trans_vs_stab_non_aggregated.csv")

### aggregate
non_agg %>% group_by(indice) %>% summarise(ave_cen_neg  = mean(cen_neg, na.rm = TRUE),
                                           ave_cen_pos = mean(cen_pos, na.rm = TRUE),
                                           sd_cen_neg  = sd(cen_neg, na.rm = TRUE),
                                           sd_cen_pos  = sd(cen_pos, na.rm = TRUE))

# model based
agg_model <- 
non_agg %>% group_by(indice, model) %>% summarise(ave_cen_neg  = mean(cen_neg, na.rm = TRUE),
                                           ave_cen_pos = mean(cen_pos, na.rm = TRUE),
                                           sd_cen_neg  = sd(cen_neg, na.rm = TRUE),
                                           sd_cen_pos  = sd(cen_pos, na.rm = TRUE))

write_csv(agg_model, "./../OUT/trans_vs_stab_only_model.csv")

#### spatial distribution
## two things to look at:
## 1) aggrement within the ensemble
## 2) mean or median change of the ensemble between transient and stabiliazed
grid_cell <-
agg_data %>% group_by(id, indice) %>% summarise(n = n(),
                                                n_sup = sum(to_save >= 0.05, na.rm = TRUE),
                                                n_na = sum(is.na(to_save)),
                                                n_neg = sum(to_save < 0.05, na.rm = TRUE))
grid_cell
## test
a <- grid_cell %>% filter(indice == "max")

a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360)
a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)

g <- ggplot(a %>% filter(n_na != 18))
g + geom_tile(aes(x = lon, y = lat, fill = n_neg))

# save independaently
write_csv(a, "./../OUT/grid_cells_max.csv")
a <- grid_cell %>% filter(indice == "min")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/grid_cells_min.csv")

a <- grid_cell %>% filter(indice == "p05")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/grid_cells_p05.csv")

a <- grid_cell %>% filter(indice == "p95")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/grid_cells_p95.csv")
