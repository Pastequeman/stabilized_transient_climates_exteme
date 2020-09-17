#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## History
# > ??/08/2020: first version
# > 10/09/2020: clean-up

## Purpose:
# Report the number of 0 in the low flow time series
# and assess Q-Q plot !
# pull all GCMs values together

## User inputs
# none!

## Libraries
library("readr")
suppressMessages(library("dplyr"))

## Jobs
for (model in c("cwatm_", "matsiro_", "clm45_", "lpjml_",
                "pcr-globwb_", "watergap2_", "h08_")) {
  print(model) ## DEBUG
  counter <- 0
  if (model == "clm45_") {
    d_count <- 30
  } else {
    d_count <- 90
  }
  # Save values here
  # high flow
  trans_max <- matrix(data = 0, nrow = d_count, ncol = 259200)
  equib_max <- matrix(data = 0, nrow = d_count, ncol = 259200)
  trans_p05 <- matrix(data = 0, nrow = d_count, ncol = 259200)
  equib_p05 <- matrix(data = 0, nrow = d_count, ncol = 259200)
  # low flows
  trans_min <- matrix(data = 0, nrow = d_count, ncol = 259200)
  equib_min <- matrix(data = 0, nrow = d_count, ncol = 259200)
  trans_p95 <- matrix(data = 0, nrow = d_count, ncol = 259200)
  equib_p95 <- matrix(data = 0, nrow = d_count, ncol = 259200)
  #
  for (gcm in c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_")) {
    print(gcm) ## DEBUG
    if (model == "clm45_" & gcm != "hadgem2-es_") {
      next
    } else if (gcm == "hadgem2-es_") {
      period1 <- "2012-2041"
      period2 <- "2147-2176"
    } else if (gcm == "ipsl-cm5a-lr_") {
      period1 <- "2015-2044"
      period2 <- "2150-2179"
    } else if (gcm == "miroc5_") {
      period1 <- "2035-2064"
      period2 <- "2243-2272"
    }
    # The raw indices
    ## High flow: Max
    file_trans_max <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period1, "_", "max.bin"),
      open = "rb")
    file_equib_max <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period2, "_", "max.bin"),
      open = "rb")
    # High flow: p05
    file_trans_p05 <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period1, "_", "p05.bin"),
      open = "rb")
    file_equib_p05 <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period2, "_", "p05.bin"),
      open = "rb")
    ## Low flow: Min
    file_trans_min <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period1, "_", "min.bin"),
      open = "rb")
    file_equib_min <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period2, "_", "min.bin"),
      open = "rb")
    # Low flow: p95
    file_trans_p95 <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period1, "_", "p95.bin"),
      open = "rb")
    file_equib_p95 <- file(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             model, gcm, "rcp26_", period2, "_", "p95.bin"),
      open = "rb")
# all the data
    for (l in seq(1, 30)) {
      counter <- counter + 1
      # high flow
      trans_max[counter, ] <- readBin(file_trans_max, "numeric",
                                      n = 259200, size = 4, endian = "little")
      equib_max[counter, ] <- readBin(file_equib_max, "numeric",
                                      n = 259200, size = 4, endian = "little")
      trans_p05[counter, ] <- readBin(file_trans_p05, "numeric",
                                      n = 259200, size = 4, endian = "little")
      equib_p05[counter, ] <- readBin(file_equib_p05, "numeric",
                                      n = 259200, size = 4, endian = "little")
      # low flow
      trans_min[counter, ] <- readBin(file_trans_min, "numeric",
                                      n = 259200, size = 4, endian = "little")
      equib_min[counter, ] <- readBin(file_equib_min, "numeric",
                                      n = 259200, size = 4, endian = "little")
      trans_p95[counter, ] <- readBin(file_trans_p95, "numeric",
                                      n = 259200, size = 4, endian = "little")
      equib_p95[counter, ] <- readBin(file_equib_p95, "numeric",
                                      n = 259200, size = 4, endian = "little")
    }
    close(file_trans_max) ; close(file_equib_max)
    close(file_equib_p05) ; close(file_trans_p05)
    close(file_trans_min) ; close(file_equib_min)
    close(file_equib_p95) ; close(file_trans_p95)
  } # done with gcms

  ### Process the data
  ### high flow only
  # Q-Q
  quant_trans_max <- apply(trans_max, 2,
                           function(x) quantile(x, seq(0, 1, by = 0.05),
                                                na.rm = TRUE, type = 8))
  quant_equib_max <- apply(equib_max, 2,
                           function(x) quantile(x, seq(0, 1, by = 0.05),
                                                na.rm = TRUE, type = 8))
  #
  quant_trans_p05 <- apply(trans_p05, 2,
                           function(x) quantile(x, seq(0, 1, by = 0.05),
                                                na.rm = TRUE, type = 8))
  quant_equib_p05 <- apply(equib_p05, 2,
                           function(x) quantile(x, seq(0, 1, by = 0.05),
                                                na.rm = TRUE, type = 8))

  #
  qq_max <- sapply(seq_len(ncol(quant_equib_max)),
                   function(j) {
                     if (sum(is.na(quant_trans_max[, j])) > 4 |
                           sum(is.na(quant_equib_max[, j])) > 4) {
                       NA
                     } else {
                       summary(lm(quant_trans_max[, j] ~ 0 +
                                    quant_equib_max[, j]))$coefficients[1]
                     }
                   })
  #
  qq_p05 <- sapply(seq_len(ncol(quant_equib_p05)),
                   function(j) {
                     if (sum(is.na(quant_trans_p05[, j])) > 4 |
                           sum(is.na(quant_equib_p05[, j])) > 4) {
                       NA
                     } else {
                       summary(lm(quant_trans_p05[, j] ~ 0 +
                                    quant_equib_p05[, j]))$coefficients[1]
                     }
                   })

  ## low flow only
  # min flow
  count_0_min <- sapply(seq_len(ncol(trans_min)),
                        function(j) {
                          if (sum(is.na(trans_min[, j])) > 4 |
                                sum(is.na(equib_min[, j])) > 4) {
                            NA
                          } else {
                            sum(trans_min[, j] == 0) +
                              sum(equib_min[, j] == 0)
                          }
                        })
  #
  # p95 flow
  count_0_p95 <- sapply(seq_len(ncol(trans_p05)),
                        function(j) {
                          if (sum(is.na(trans_p95[, j])) > 4 |
                                sum(is.na(equib_p95[, j])) > 4) {
                            NA
                          } else {
                            sum(trans_p95[, j] == 0) +
                              sum(equib_p95[, j] == 0)
                          }
                        })

  # done with data processing
  # reformat as data-frame
  # low flow: min
  temp <- tibble(lon = rep(seq(-179.75, 179.75, 0.5), times = 360),
                 lat = rep(seq(89.75, -89.75, -0.5), each = 720),
                 data = count_0_min)
  # remove the sea
  count_0_min <- temp[!is.na(temp$data), ]
  # low flow: p95
  temp <- tibble(lon = rep(seq(-179.75, 179.75, 0.5), times = 360),
                 lat = rep(seq(89.75, -89.75, -0.5), each = 720),
                 data = count_0_p95)
  # remove the sea
  count_0_p95 <- temp[!is.na(temp$data), ]
  #
  # high flow: max
  temp <- tibble(lon = rep(seq(-179.75, 179.75, 0.5), times = 360),
                 lat = rep(seq(89.75, -89.75, -0.5), each = 720),
                 data = qq_max)
  # remove the sea
  qq_max <- temp[!is.na(temp$data), ]
  # high flow: p05
  temp <- tibble(lon = rep(seq(-179.75, 179.75, 0.5), times = 360),
                 lat = rep(seq(89.75, -89.75, -0.5), each = 720),
                 data = qq_p05)
  # remove the sea
  qq_p05 <- temp[!is.na(temp$data), ]
  # save the dataset
  ## Low flow data
  write_csv(count_0_min, paste0("./../OUT/count_0_", model, "min.csv"))
  write_csv(count_0_p95, paste0("./../OUT/count_0_", model, "p95.csv"))
  ## High flow data
  write_csv(qq_max, paste0("./../OUT/qq_", model, "max.csv"))
  write_csv(qq_p05, paste0("./../OUT/qq_", model, "p05.csv"))
} # models
