#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

## Purpose:
# Assess the influence of natural variability on the detection of significance between 
# extreme indice for two periods

## Libraries
library("readr")
suppressMessages(library("dplyr"))

## Read files. Since there are a lot of combinations it is safer to split the files > I made 8 regions
sections <- c(1, 32400, 64800, 97200, 129600, 162000, 194400, 226800, 259200)

for (sec in 1:8) {
  for (m in c("cwatm_", "matsiro_", "clm45_", "lpjml_", "pcr-globwb_", "watergap2_", "h08_")) { # add more as they come
    for (gcms in c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_")) {
      if ((m == "clm45_" & gcms %in% c("ipsl-cm5a-lr_", "miroc5_"))) {
        next
      }
      # now the different pi sets
      for (i in c("pi12", "pi13", "pi14", "pi15")) { # !! Add remaining later
        if (i == "pi12") {
          if (m == "matsuro_" & gcms == "hadgem2-es_") { next }
          period1 <- "1661-1690"
          period2 <- "1691-1720"
        } else if (pi == "pi13") {
          period1 <- "1661-1690"
          period2 <- "1721-1750"
        } else if (pi == "pi14") {
          period1 <- "1661-1690"
          period2 <- "1751-1780"
        } else if (pi == "pi15") {
          period1 <- "1661-1690"
          period2 <- "1781-1810"
        } else if (pi == "pi16") {
          period1 <- "1661-1690"
          period2 <- "1811-1840"        
        }
        # get the proper files
        for (id in c("max", "min", "p05", "p95")) {
          temp <- read_csv(paste0("/data01/julien/projects/extreme_trans_stab/OUT/", m, gcms, "picontrol_",
                                  period1, "_picontrol_", period2, "_", id, ".csv"), col_types = cols(to_save = col_double()))
          temp$id     <- seq(1,259200)
          temp$indice <- id
          temp$model  <- m
          temp$gcm    <- gcms
          temp <- temp[sections[sec]:sections[(sec+1)], ] ## keep only subset to save memmory


          if (gcms == "hadgem2-es_" & m == "cwatm_" & id == "max") {
            agg_data <- temp
          } else {
            agg_data <- rbind(agg_data, temp)
          }
        } # indice loop
      } # pi period loop
    } # gcm loop
  } # models loop
  # done for a given region > process the agg_data file to get the upper and lower bounds
  rm(temp)
  # analyis start, here no aggregation yet
  non_agg <- 
    agg_data %>% group_by(model, gcm, indice) %>% summarise(n    = n(),
                                                          n_pos = sum(to_save >= 0.05, na.rm = TRUE),
                                                          n_na  = sum(is.na(to_save)))
  non_agg$n_neg    <- (non_agg$n - non_agg$n_na) - non_agg$n_pos
  non_agg$cen_neg  <- non_agg$n_neg / (non_agg$n - non_agg$n_na) * 100
  non_agg$cen_pos  <- non_agg$n_pos / (non_agg$n - non_agg$n_na) * 100


  # save and clean
  write_csv(agg_model, paste0("./../OUT/picontrol_nonaggr_reg_", sec, ".csv"))
  rm(non_agg)
}  # regions


