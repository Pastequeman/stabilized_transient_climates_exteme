#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## Purpose
# extract daily discharge
# and extreme values (high and low)
# from specific locations

## History:
# ??/??/2020: created
# 09/09/2020: big clean-up

## Libraries:
library("readr")
suppressMessages(library("dplyr"))

#### Code #####
# Amazon, Mississipi, Niger, and Murray
basins_id <- c(1, 3, 11, 19)

# find the river mouth of these basins
f_riv <- file("./../IN/riv_num_/rivnum.WATCH.hlf", open = "rb")
riv_number <- readBin(f_riv, what = "numeric", n = 259200,
                      size = 4, endian = "big")
close(f_riv)

# seq (river sequence)
f_riv <- file("./../IN/riv_seq_/rivseq.WATCH.hlf", open = "rb")
riv_seq <- readBin(f_riv, what = "numeric", n = 259200,
                   size = 4, endian = "big")
close(f_riv)

#nxl (L of the next cell)
f_riv <- file("./../IN/riv_nxl_/rivnxl.WATCH.hlf", open = "rb")
riv_nxl <- readBin(f_riv, what = "double", n = 259200,
                   size = 4, endian = "big")
close(f_riv)

## Merge and clean
riv <- tibble(L   = seq(1, 259200),
              lon = rep(seq(-179.75, 179.75, 0.5), times = 360),
              lat = rep(seq(89.75, -89.75, -0.5), each = 720),
              num = riv_number,
              seq = riv_seq,
              nxl = riv_nxl)
rm(riv_number) ; rm(riv_nxl) ; rm(riv_seq) ; rm(f_riv)

cell_id <- c(0, 0, 0, 0)
cc <- 0
for (i in basins_id) {
  cc <- cc + 1
  cell_id[cc] <-
  riv[riv$num == i &
        riv$seq == max(riv[riv$num == i, ]$seq), ]$L
}
# clean
rm(cc) ; rm(riv)

## Main jobs
for (L in cell_id) {
  for (model in c("cwatm_", "matsiro_", "clm45_",
                  "lpjml_", "pcr-globwb_", "watergap2_", "h08_")) {
    for (gcm in c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_")) {
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
      # extreme for this cell
      # min
      charac_min <- read_csv(paste0("../OUT/", model, gcm,
                                    "rcp26_", period1,
                                    "_rcp26_", period2, "_", "min.csv"),
                             col_types = cols(to_save = col_double(),
                                              stat    = col_double(),
                                              low     = col_double(),
                                              hgh     = col_double())
                             )
      # max
      charac_max <- read_csv(paste0("../OUT/", model, gcm,
                                    "rcp26_", period1,
                                    "_rcp26_", period2, "_", "max.csv"),
                             col_types = cols(to_save = col_double(),
                                              stat    = col_double(),
                                              low     = col_double(),
                                              hgh     = col_double())
                             )
      # prep data.frame to save
      to_save <- tibble(model = model,
                        gcm = gcm,
                        L = L,
                        trans_min = rep(0, 30),
                        equ_min   = rep(0, 30),
                        prop_min  = ifelse(charac_min[L, ]$stat <
                                             charac_min[L, ]$low |
                                             charac_min[L, ]$stat >
                                             charac_min[L, ]$hgh, 1, 0),
                        trans_max = rep(0, 30),
                        equ_max   = rep(0, 30),
                        prop_max  = ifelse(charac_max[L, ]$stat <
                                             charac_max[L, ]$low |
                                             charac_max[L, ]$stat >
                                             charac_max[L, ]$hgh, 1, 0))
      #
      ## Max
      f <- paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
                  model, gcm, "rcp26_", period1, "_", "max.bin")
      file_trans_max <- file(f, open = "rb")
      f <- paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
                  model, gcm, "rcp26_", period2, "_", "max.bin")
      file_equ_max   <- file(f, open = "rb")
      ## Min
      f <- paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
                  model, gcm, "rcp26_", period1, "_", "min.bin")
      file_trans_min <- file(f, open = "rb")
      f <- paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
                  model, gcm, "rcp26_", period2, "_", "min.bin")
      file_equ_min   <- file(f, open = "rb")
      ## clean
      rm(f)
      ## Save
      for (l in seq(1, 30)) {
        t1 <- readBin(file_trans_max, "numeric", n = 259200,
                      size = 4, endian = "little")
        t2 <- readBin(file_equ_max, "numeric", n = 259200,
                      size = 4, endian = "little")
        t3 <- readBin(file_trans_min, "numeric", n = 259200,
                      size = 4, endian = "little")
        t4 <- readBin(file_equ_min, "numeric", n = 259200,
                      size = 4, endian = "little")
        to_save[l, 4] <- t3[L]
        to_save[l, 5] <- t4[L]
        to_save[l, 7] <- t1[L]
        to_save[l, 8] <- t2[L]
      }
      close(file_trans_max) ; close(file_equ_max)
      close(file_trans_min) ; close(file_equ_min)
      # append
      write_csv(to_save, "../OUT/distribution_v2.csv", append = TRUE)
    }
  }
}
