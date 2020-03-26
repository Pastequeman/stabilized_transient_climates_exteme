#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

## Purpose:
# This is evaluating the 95th of high flow only
#                        05th of low flow
# read the .csv files generated for all scenarios
# report the following:
#    - nb and % of grid cell where sign. detected
#    - for all models and gcm combinations

## Libraries
library("readr")
suppressMessages(library("dplyr"))

## read files
for (m in c("cwatm_", "matsiro_", "clm45_", "lpjml_", "pcr-globwb_", "watergap2_", "h08_")) { # add more as they come
  for (gcms in c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_")) {
    if ((m == "clm45_" & gcms %in% c("ipsl-cm5a-lr_", "miroc5_"))) {
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
    for (id in c("max","p05", "min", "p95")) { #  
      temp <- read_csv(paste0("/data01/julien/projects/extreme_trans_stab/OUT/high_value_", m, gcms, "rcp26_",
                              period1, "_rcp26_", period2, "_", id, ".csv"), col_types = cols(to_save = col_double(),
                                                                                              low = col_double(),
                                                                                              hgh = col_double()))
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
rm(temp) ; rm(period1) ; rm(period2) ; rm(id) ; rm(m) ; rm(gcms)

#### aggregate (all globe)
non_agg <- 
agg_data %>% group_by(model, gcm, indice) %>% summarise(n    = n(),
                                                      n_na  = sum(is.na(to_save)),
                                                      n_pos = sum(to_save < low | to_save > hgh, na.rm = TRUE))

non_agg$n_neg   <- (non_agg$n - non_agg$n_na) - non_agg$n_pos        # ok
non_agg$cen_neg <- non_agg$n_neg / (non_agg$n - non_agg$n_na) * 100  # ok 
non_agg$cen_pos <- non_agg$n_pos / (non_agg$n - non_agg$n_na) * 100  # ok

non_agg

write_csv(non_agg, "./../OUT/process/trans_vs_stab_non_aggregated_tail.csv")

### aggregate (all)
to_save <- non_agg %>% group_by(indice) %>% summarise(ave_cen_neg  = mean(cen_neg, na.rm = TRUE),
                                                    ave_cen_pos = mean(cen_pos, na.rm = TRUE),
                                                    sd_cen_neg  = sd(cen_neg, na.rm = TRUE),
                                                    sd_cen_pos  = sd(cen_pos, na.rm = TRUE))
to_save
write_csv(to_save, "./../OUT/process/condenced_table_tail.csv") ; rm(to_save)
# model based
agg_model <- 
non_agg %>% group_by(indice, model) %>% summarise(ave_cen_neg  = mean(cen_neg, na.rm = TRUE),
                                                  ave_cen_pos = mean(cen_pos, na.rm = TRUE),
                                                  sd_cen_neg  = sd(cen_neg, na.rm = TRUE),
                                                  sd_cen_pos  = sd(cen_pos, na.rm = TRUE))
write_csv(agg_model, "./../OUT/process/trans_vs_stab_only_model_tail.csv")

# GCM only
agg_model <- 
non_agg %>% group_by(indice, gcm) %>% summarise(ave_cen_neg  = mean(cen_neg, na.rm = TRUE),
                                                ave_cen_pos = mean(cen_pos, na.rm = TRUE),
                                                sd_cen_neg  = sd(cen_neg, na.rm = TRUE),
                                                sd_cen_pos  = sd(cen_pos, na.rm = TRUE))
write_csv(agg_model, "./../OUT/process/trans_vs_stab_only_gcm_tail.csv")

#### spatial distribution
## two things to look at:
## 1) aggrement within the ensemble
## 2) mean or median change of the ensemble between transient and stabiliazed
agg_data$n_pos <- ifelse(agg_data$to_save < agg_data$low | agg_data$to_save > agg_data$hgh, 1, 0)
agg_data$n_neg <- ifelse(agg_data$to_save >= agg_data$low & agg_data$to_save <= agg_data$hgh, 1, 0) 
agg_data$n_na  <- ifelse(is.na(agg_data$to_save), 1, 0)

grid_cell <- agg_data %>% group_by(id, indice) %>% summarise(n = n(), pos = sum(n_pos, na.rm = TRUE), neg = sum(n_neg, na.rm = TRUE), na = sum(n_na, na.rm = TRUE))
## test
a <- grid_cell %>% filter(indice == "max")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/process/grid_cells_max_tail.csv")

a <- grid_cell %>% filter(indice == "min")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/process/grid_cells_min_tail.csv")

a <- grid_cell %>% filter(indice == "p05")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/process/grid_cells_p05_tail.csv")

a <- grid_cell %>% filter(indice == "p95")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360) ; a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/process/grid_cells_p95_tail.csv")


## additional analysis to distinguish patterns
# gcm
grid_cell <-
  agg_data %>% group_by(id, indice, gcm) %>% summarise(n      = n(),
                                                       n_sup = sum(n_pos, na.rm = TRUE),
                                                       n_na = sum(is.na(to_save)),
                                                       n_neg = sum(n_neg, na.rm = TRUE))
write_csv(grid_cell, "./../OUT/process/grid_cells_gcms_all_ind_tail.csv")

# models
grid_cell <-
  agg_data %>% group_by(id, indice, model) %>% summarise(n      = n(),
                                                       n_sup = sum(n_pos, na.rm = TRUE),
                                                       n_na = sum(is.na(to_save)),
                                                       n_neg = sum(n_neg, na.rm = TRUE))
write_csv(grid_cell, "./../OUT/process/grid_cells_models_all_ind_tail.csv")
