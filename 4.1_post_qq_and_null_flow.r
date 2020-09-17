#!/usr/bin/env Rscript
args <-  commandArgs(trailingOnly = TRUE)

## History:
# ??/08/2020: first version
# 10/09/2020: clean-up and save dataset
# 14/09/2020: fix spelling and bug fix in output

## Purpose:
# read files indicating the number of low flow equal
# to zero for different periods

## Libraries
library("readr")
suppressMessages(library("dplyr"))

#### Jobs ####
# read
for (model in c("cwatm_", "matsiro_", "clm45_", "lpjml_",
                "pcr-globwb_", "watergap2_", "h08_")) {
  temp_min <- read_csv(paste0("./../OUT/count_0_",
                              model, "min.csv"),
                       col_types = cols(lon  = col_double(),
                                        lat  = col_double(),
                                        data = col_double()))
  temp_p95 <- read_csv(paste0("./../OUT/count_0_",
                              model, "p95.csv"),
                       col_types = cols(lon  = col_double(),
                                        lat  = col_double(),
                                        data = col_double()))
  temp_min$model <- model ; temp_p95$model <- model
  temp_min$ind   <- "min" ; temp_p95$ind   <- "p95"
  #
  temp_max <- read_csv(paste0("./../OUT/qq_", model, "max.csv"),
                       col_types = cols(lon  = col_double(),
                                        lat  = col_double(),
                                        data = col_double()))
  temp_p05 <- read_csv(paste0("./../OUT/qq_", model, "p05.csv"),
                       col_types = cols(lon  = col_double(),
                                        lat  = col_double(),
                                        data = col_double()))
  temp_max$model <- model ; temp_p05$model <- model
  temp_max$ind   <- "max" ; temp_p05$ind   <- "p05"
  #
  if (model == "cwatm_") {
    all_data_low <- rbind(temp_min, temp_p95)
    all_data_hgh <- rbind(temp_max, temp_p05)
  } else {
    temp <- rbind(temp_min, temp_p95)
    all_data_low <- rbind(all_data_low, temp)
    temp <- rbind(temp_max, temp_p05)
    all_data_hgh <- rbind(all_data_hgh, temp)
  }
}
# clean
rm(temp) ; rm(temp_min) ; rm(temp_p95) ; rm(model)
rm(temp_max) ; rm(temp_p05)

# create a file that will indicate the number of model to remove
# for every grid-cell due to flow always being equal to zero
all_data_low$mul <- ifelse(all_data_low$model == "clm45_",
                           1, 3)
all_data_low$pos <- ifelse(all_data_low$data == 180 &
                           all_data_low$model != "clm45_" |
                           all_data_low$data == 30 &
                           all_data_low$model == "clm45_",
                           all_data_low$mul, 0)

adjust_low_flow <-
  all_data_low %>%
  group_by(lon, lat, ind) %>%
  summarise(n = sum(mul, na.rm = TRUE),
            m = sum(pos, na.rm = TRUE))

# save
write_csv(adjust_low_flow, "./../OUT/adjust_low_flow.csv")

## Similarly, extract locations where QQ plot
# are equal to 1
all_data_hgh$mul <- ifelse(all_data_hgh$model == "clm45_",
                           1, 3)
all_data_hgh$pos <- ifelse(all_data_hgh$data >= 1.75 |
                             all_data_hgh$data <= 0.25,
                           all_data_hgh$mul, 0)

adjust_high_flow <-
  all_data_hgh %>%
  group_by(lon, lat, ind) %>%
  summarise(n = sum(mul, na.rm = TRUE),
            m = sum(pos, na.rm = TRUE))

# save
write_csv(adjust_high_flow, "./../OUT/adjust_high_flow.csv")
