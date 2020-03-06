#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

## Purpose:
# load the gembel distro
# compute stab - trans
# and see to wha extend it is different

## Libraries
library("readr")
suppressMessages(library("dplyr"))

## Read data
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
    for (id in c("max", "min", "p05", "p95")) {
      temp <- read_csv(paste0("/data01/julien/projects/extreme_trans_stab/OUT/gumbel/", m,
                              gcms, "rcp26_",period1, "_", id, ".csv"), col_types = cols(alpha = col_double(),
                                                                                espilon = col_double(),
                                                                                y100 = col_double()))
      temp$id     <- seq(1,259200)
      temp$indice <- id
      temp$model  <- m
      temp$gcm    <- gcms
      temp$climate <- "transient"
      if (gcms == "hadgem2-es_" & m == "cwatm_" & id == "max") {
        tran_data <- temp
      } else {
        tran_data <- rbind(tran_data, temp)
      }
      
      temp <- read_csv(paste0("/data01/julien/projects/extreme_trans_stab/OUT/gumbel/", m,
                              gcms, "rcp26_", period2, "_", id, ".csv"), col_types = cols(alpha = col_double(),
                                                                                espilon = col_double(),
                                                                                y100 = col_double()))
      temp$id     <- seq(1,259200)
      temp$indice <- id
      temp$model  <- m
      temp$gcm    <- gcms
      temp$climate <- "stabilized"
      if (gcms == "hadgem2-es_" & m == "cwatm_" & id == "max") {
        stab_data <- temp
      } else {
        stab_data <- rbind(stab_data, temp)
      }
    }
  }
}
# clean
rm(temp) ; rm(gcms) ; rm(m) ; rm(period1) ; rm(period2) ; rm(id)


## make diff file
diff_trans_minus_stab <- tibble(id = stab_data$id,
                                indice = stab_data$indice,
                                model = stab_data$model,
                                gcm = stab_data$gcm,
                                diff = tran_data$y100 - stab_data$y100)
# take median and mean for all indices
aggregated_diff <- 
diff_trans_minus_stab %>% group_by(id, indice) %>% summarise(ave_y100 = mean(diff, na.rm = TRUE),
                                                            med_y100 = median(diff, na.rm = TRUE) )
## save
a <- aggregated_diff %>% filter(indice == "max")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360)
a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/agg_100y_diff_max.csv")

a <- aggregated_diff %>% filter(indice == "min")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360)
a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/agg_100y_diff_min.csv")

a <- aggregated_diff %>% filter(indice == "p95")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360)
a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/agg_100y_diff_p95.csv")

a <- aggregated_diff %>% filter(indice == "p05")
a$lon <- rep(seq(-179.75, by = 0.5, length.out = 720), times = 360)
a$lat <- rep(seq(89.75, by = -0.5, length.out = 360), each = 720)
write_csv(a, "./../OUT/agg_100y_diff_p05.csv")
