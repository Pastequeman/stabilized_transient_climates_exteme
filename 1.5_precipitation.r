#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## Purpose:
# read the pcp files for a given period

## user inputs
if (length(args) != 3) {
  stop("Wrong nb of argument:\n
           GCM   -> HadGEM2-ES / IPSL-CM5A-LR / MIROC5
           PERIOD -> YYYY YYYY")
} else {
  gcm    <- args[1]
  period <- as.integer(c(args[2], args[3]))
}


## Libraries
library("tibble")
library("readr")
suppressMessages(library("dplyr"))
library("ncdf4")

## mask
file_mask <- nc_open(
  "/data01/julien/models/GW_ISIMIP/ISIMIP2b_landseamask_generic.nc4")
seamask <- ncvar_get(file_mask, "LSM")
seamask <- ifelse(is.na(seamask), FALSE, TRUE)
nc_close(file_mask)


## initialization
seq_start <- seq(2011, 2291, 10)
i <- 0

# global ave. of these
output    <- tibble(years    = seq(period[1], period[2]),
                    mean_pcp = double(30),
                    cdd      = double(30),
                    cwd      = double(30),
                    d10      = double(30),
                    d20      = double(30),
                    wdp      = double(30))

## Main loop
for (years in seq(period[1], period[2], 1)) {
  i <- i + 1
  if (years <= 2010) {
    y_start <- 2006
    y_end   <- 2010
    yrs_str <- "20060101-20101231"
  } else {
    # find closest years + lower which end in xxx1
    y_start <- years
    while (! y_start %in% seq_start) {
      y_start <- y_start - 1
    }
    y_end <- y_start + 9
    yrs_str <- paste0(y_start, "0101-",
                      y_end, "1231")
  }

  # the file to open
  file_to_open <-
  nc_open(paste0("/data01/julien/models/GW_ISIMIP/pcp/rcp26/pr_day_",
                 gcm, "_rcp26_r1i1p1_EWEMBI_",
                 yrs_str, ".nc4"))
  # where to start
  if (years == y_start) {
    file_st <- 1
  } else {
    file_st <-
        seq(
            as.Date(paste0(y_start, "/1/1")),
            as.Date(paste0((years - 1), "/12/31")),
            "day")
  }
  file_stop <-
      seq(as.Date(paste0(years, "/1/1")),
          as.Date(paste0(years, "/12/31")),
          "day")
  # get the data
  pcp <-
  ncvar_get(file_to_open, "pr",
            start = c(1, 1, (length(file_st) + 1)),
            count = c(-1, -1, length(file_stop)))
  nc_close(file_to_open)

  # remove the sea and adjust unit to mm/day
  # also new matrix with no sea
  lnd_mat <- array(data = 0, dim = c(67420, dim(pcp)[3]))
  for (j in seq(1, dim(pcp)[3])) {
    lnd_mat[, j] <- matrix(pcp[, , j][as.vector(seamask)]
                             * 86400,
                             nrow = 67420,
                             byrow = TRUE)
  }

  # transform, annual ave. in all cell
  annual_ave <- apply(lnd_mat, 1, mean)
  # consecutive dry day
  cdd <- apply(lnd_mat, 1, function(x) {
    max(rle(x < 1)$length[rle(x <= 1)$values])
  })
  cdd <- ifelse(!is.finite(cdd), 0, cdd)
  # concecutive wet day everywhere
  ccwet <- apply(lnd_mat, 1, function(x) {
    max(rle(x >= 1)$length[rle(x >= 1)$values])
  })
  ccwet <- ifelse(!is.finite(ccwet), 0, ccwet)
  # annual count pcp >= 10
  annual_10 <- apply(lnd_mat, 1, function(x) {
    sum(x > 10)
  })
  # annual count pcp >= 20
  annual_20 <- apply(lnd_mat, 1, function(x) {
    sum(x > 20)
  })
  # annual wet-day pcp
  annual_wdp <- apply(lnd_mat, 1, function(x) {
    sum(x[x >= 1])
  })
  # save lnd dis and take average
  if (years == period[1]) {
    temp_ave   <- annual_ave
    temp_cdd   <- cdd
    temp_ccwet <- ccwet
    temp_10    <- annual_10
    temp_20    <- annual_20
    temp_wdp   <- annual_wdp
  } else {
    temp_ave   <- (temp_ave + annual_ave) / 2
    temp_cdd   <- (temp_cdd + cdd) / 2
    temp_ccwet <- (temp_ccwet + ccwet) / 2
    temp_10    <- (temp_10 + annual_10) / 2
    temp_20    <- (temp_20 + annual_20) / 2
    temp_wdp   <- (temp_wdp + annual_wdp) / 2
  }

  # save
  output[i, 2:7] <- c(mean(annual_ave, na.rm = TRUE),
                      mean(cdd, na.rm = TRUE),
                      mean(ccwet, na.rm = TRUE),
                      mean(annual_10, na.rm = TRUE),
                      mean(annual_20, na.rm = TRUE),
                      mean(annual_wdp, na.rm = TRUE))
}

# save the final outputs
write_csv(output,
          paste0("../OUT/pcp_", gcm, "_",
                 period[1], "-", period[2],
                 ".csv"))

# write the temp files
write_csv(enframe(temp_ave),
          paste0("../OUT/pcp_ave", gcm, "_",
                 period[1], "-", period[2],
                 ".csv"))
write_csv(enframe(temp_cdd),
          paste0("../OUT/pcp_cdd", gcm, "_",
                 period[1], "-", period[2],
                 ".csv"))
write_csv(enframe(temp_ccwet),
          paste0("../OUT/pcp_ccwet", gcm, "_",
                 period[1], "-", period[2],
                 ".csv"))
write_csv(enframe(temp_10),
          paste0(
            "../OUT/pcp_10", gcm, "_",
            period[1], "-", period[2],
            ".csv"))
write_csv(enframe(temp_20),
          paste0(
            "../OUT/pcp_20", gcm, "_",
            period[1], "-", period[2],
            ".csv"))
write_csv(enframe(temp_wdp),
          paste0("../OUT/pcp_wdp", gcm, "_",
                 period[1], "-", period[2],
                 ".csv"))
