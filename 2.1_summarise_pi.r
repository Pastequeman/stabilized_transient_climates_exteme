#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## Purpose
# retrieve the stat for the pi-control
# construct 95th interval conf
## history
# 15/06/2020: added support for bootstrap
# 08/09/2020: implemented Wilk and clean-up
## Libraries
library("readr")
library("tibble")

## input data
# land area
ar <- file("./../IN/grdara.hlf", open = "rb")
t <- readBin(ar, what = "numeric", n = 259200, size = 4, endian = "big")
close(ar) ; rm(ar)

## global data
for (models in c("h08", "watergap2", "lpjml", "matsiro",
                 "cwatm", "pcr-globwb", "clm45")) {
  for (gcms in c("hadgem2-es", "miroc5", "ipsl-cm5a-lr")) {
    if (gcms %in% c("miroc5", "ipsl-cm5a-lr") & models == "clm45") {
      next
    }
    # list all the file concerned
    file <- list.files("./../OUT/", pattern = paste0("^", models, "_",
                                                     gcms, "_picontrol"))
    file_boot <- list.files("./../OUT/pi_process",
                            pattern = paste0("^", models,
                                             "_", gcms, "_picontrol"))
    #
    for (ind in c("max", "p05", "min", "p95")) {
      sub_file <- file[grep(pattern = ind, file)]
      sub_file_boot <- file_boot[grep(pattern = ind, file_boot)]
      temp_matrix <- matrix(data = 0, nrow = 259200,
                            ncol = (length(sub_file)
                              + 100 * length(sub_file_boot)))
      #
      # Step 1: Default files
      for (i in seq_len(length(sub_file))) {
        temp <- read_csv(paste0("./../OUT/", sub_file[i]),
                         col_types = cols(to_save = col_double(),
                                          stat    = col_double(),
                                          low     = col_double(),
                                          hgh     = col_double()))
        temp_matrix[, i] <- temp$to_save
      }
      #
      # Step 2: bootstrap data
      for (i in seq_len(length(sub_file_boot))) {
        file_to_r <- file(paste0("./../OUT/pi_process/",
                                 sub_file_boot[i]), open = "rb")
        for (j in 1:100) {
          temp_matrix[, (length(sub_file) + i * 100 - 100 + j)] <-
            readBin(file_to_r, size = 4, what = "numeric",
                    n  = 259200, endian = "little")
        }
        close(file_to_r)
      }
      #
      # done filling the matrix for that combination of model, gcm, and indice
      low <- apply(temp_matrix[, seq_len(length(sub_file))], 1,
                   function(x) quantile(x, 0.025, na.rm = TRUE))
      hgh <- apply(temp_matrix[, seq_len(length(sub_file))], 1,
                   function(x) quantile(x, 0.975, na.rm = TRUE))
      #
      temp <- tibble(L    = seq(1, 259200),
                     low  = low,
                     hgh  = hgh,
                     area = t)
      ## Wilk 2016
      # recompute the FDR for every col of the matrix
      # then evaluate global significance
      # report finaly the 95 percentage
      # total land surface:
      tot_lnd_s <- sum(temp[!is.na(temp_matrix[, i]), ]$area)
      for (i in seq_len(length(sub_file))) {
        #
        temp_2 <- temp_matrix[, i]
        temp_2 <- temp_2[!is.na(temp_2) & temp_2 != 0] # remove the sea
        temp_2 <- sort(temp_2, decreasing = FALSE)
        # The cells where this relation is fullfill
        cell <- which(temp_2 <= (seq(1, length(temp_2)) / length(temp_2) * 0.1))
        if (length(cell) == 0) {
          next
        }
        # Pick the max p value belonging for such cell
        sig <- (100 * sum(temp[temp_matrix[, i] < (temp_2[cell[length(cell)]]) &
                               !is.na(temp_matrix[, i]), ]$area) / tot_lnd_s)
        tt <-
        tibble(model    = models,
               gcm      = gcms,
               it       = i,
               indice   = ind,
               sig_area = sig)
        write_csv(tt, "./../OUT/pi_control_summary.csv", append = TRUE)
        rm(tt)
      }
      #
      ## clean and save
      rm(low) ; rm(hgh) ; rm(temp_matrix)
      # save
      write_csv(temp, paste0("./../OUT/pi_process/range_pi_stat_",
                             models, "_", gcms, "_", ind, ".csv"))
    }
  }
}
