#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## Purpose:
# read the .csv files generated for all scenarios
# report the following:
#    - nb and % of grid cell where sign. detected
#    - for all models and gcm combinations
#    - for main table
## history:
# created: ?
# implemented Wilk criteria + clean-up: 08-09-2020

## Libraries
library("readr")
suppressMessages(library("dplyr"))
#
# land area
ar <- file("./../IN/grdara.hlf", open = "rb")
t <- readBin(ar, what = "numeric", n = 259200, size = 4, endian = "big")
close(ar) ; rm(ar)

## read main files
for (m in c("cwatm_", "matsiro_", "clm45_",
            "lpjml_", "pcr-globwb_", "watergap2_", "h08_")) {
  print(m) ## DEBUG
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
      temp <- read_csv(paste0("/data01/julien/projects/extreme_trans_stab/OUT/",
                              m, gcms, "rcp26_",
                              period1, "_rcp26_", period2, "_", id, ".csv"),
                       col_types = cols(to_save = col_double(),
                                        stat    = col_double(),
                                        low     = col_double(),
                                        hgh     = col_double())
                       )
      # Wilks
      all_stat <- temp$to_save[!is.na(temp$to_save)]
      criteria <- sort(all_stat)[length(which(sort(all_stat) <
                                                seq(1, length(all_stat)) /
                                                length(all_stat) * 0.1))]
      temp$sign1  <- ifelse(temp$to_save < criteria, 1, 0)
      #
      # other usefull info
      temp$area   <- t
      temp$id     <- seq(1, 259200)
      temp$indice <- id
      temp$model  <- m
      temp$gcm    <- gcms
       #
      if (gcms == "hadgem2-es_" & m == "cwatm_" & id == "max") {
        agg_data <- temp
      } else {
        agg_data <- rbind(agg_data, temp)
      }
    }
  }
}
rm(temp) ; rm(period1) ; rm(period2) ; rm(id) ; rm(m) ; rm(gcms) ; rm(all_stat)
rm(criteria)

# add informations to the dataset
# Koppen Geiger
koppen <- read_table2("./../IN/koppen.txt")
koppen$Cls[koppen$Lat == 20.75 & koppen$Lon == 17.25] <- "BWh"
# remove "BWh" and "EF"
temp <- tibble(lon = rep(seq(-179.75, 179.75, by = 0.5), times = 360),
               lat = rep(seq(89.75, -89.75, by = -0.5), each = 720),
               L   = seq(1, 259200, 1))
koppen <- koppen %>% left_join(temp, by = c("Lon" = "lon", "Lat"  = "lat"))
rm(temp)

# Remove the sea
agg_data <- agg_data %>% filter(!is.na(sign1))

# All stats for table 3
all_var <-
agg_data %>%
group_by(indice, model, gcm) %>%
summarise(n        = n(),
          tot_area = sum(area, na.rm = TRUE),
          wilk     = 100 * sum(sign1 * area, na.rm = TRUE) / tot_area)
#
## GWM based
gwm_based <-
all_var %>%
  group_by(indice, model) %>%
  summarise(ave_sig = mean(wilk, na.rm = TRUE),
            sd_sig  = sd(wilk, na.rm = TRUE))

## GCM based
gcm_based <-
all_var %>%
  group_by(indice, gcm) %>%
  summarise(ave_sig = mean(wilk, na.rm = TRUE),
            sd_sig  = sd(wilk, na.rm = TRUE))

## indice based
ind_based <-
all_var %>%
  group_by(indice) %>%
  summarise(ave_sig = mean(wilk, na.rm = TRUE),
            sd_sig  = sd(wilk, na.rm = TRUE))

## Save
write_csv(gwm_based, "./../OUT/process/gwm_based.csv")
write_csv(gcm_based, "./../OUT/process/gcm_based.csv")
write_csv(ind_based, "./../OUT/process/ind_based.csv")
write_csv(agg_data, "./../OUT/process/all_main_data.csv")
