#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## Purpose:
# for two 30-year period, compute extreme indices
# for discharge and check whereas they are different
# at each grid cells
# use directly the netcdf files
#
# 1) Read required files for first period
# 2) Compute extreme indices for first period
# 3) Read file for second period
# 4) Compute extreme indices for second period
# 5) perform t test > write out
# 6) also perform gumbel fitting and extract parameters

## CC inputs
if (length(args) != 9) {
  stop("Wrong nb of argument:\n
           MODEL -> CLM45 / cwatm / h08 / watergap2 /
                    pcr-globwb / lpjml / matsiro
           GCM   -> hadgem2-es / ipsl-cm5a-lr / miroc5
           EXP1  -> hist / picontrol / rcp26 / rcp60
           EXP2  -> hist / picontrol / rcp26 / rcp60
           PERIOD1 -> YYYY YYYY
           PERIOD2  -> YYYY YYYY
           PERMUTATION -> yes/no")
} else {
  model       <- args[1]
  gcm         <- args[2]
  experiments <- c(args[3], args[4])
  period_1    <- as.integer(c(args[5], args[6]))
  period_2    <- as.integer(c(args[7], args[8]))
  permutation <- args[9]
}

## here set random period for pi-control
if (experiments[1] == "picontrol" |
      experiment[2] == "picontrol") {
  if (model == "matsiro" & gcm == "hadgem2-es") {
    # because some data are not available
    a <- seq(1671, 1860)
  } else {
    a <- seq(1661, 1860)
  }
  st1 <- sample(a, 1)
  while (st1 >= 1832) {
    st1 <- sample(a, 1)
  }
  st2 <- sample(a, 1)
  while ((st2 >= 1832) | (st2 >= st1 & st2 <= (st1 + 29)) |
           ((st2 + 29) <= (st1 + 29) & (st2 + 29) <=
              (st1 + 29))) {
                st2 <- sample(a, 1)
              }
  period_1 <- c(st1, (st1 + 29))
  period_2 <- c(st2, (st2 + 29))
  # clean
  rm(a) ; rm(st1) ; rm(st2)
  print(c(period_1, period_2)) ; print(model)
}



## Libraries
library("ncdf4")
library("readr")
suppressMessages(library("dplyr"))


## The matrix storing storing the extreme indices
outputs <- list(max_dis_v1 = matrix(data = 0, nrow = 30, ncol = 259200),
                p05_dis_v1 = matrix(data = 0, nrow = 30, ncol = 259200),
                min_dis_v1 = matrix(data = 0, nrow = 30, ncol = 259200),
                p95_dis_v1 = matrix(data = 0, nrow = 30, ncol = 259200),
                max_dis_v2 = matrix(data = 0, nrow = 30, ncol = 259200),
                p05_dis_v2 = matrix(data = 0, nrow = 30, ncol = 259200),
                min_dis_v2 = matrix(data = 0, nrow = 30, ncol = 259200),
                p95_dis_v2 = matrix(data = 0, nrow = 30, ncol = 259200)
                )
#
leap_years <- c(1664, 1668, 1672, 1676, 1680, 1684, 1688, 1692, 1696, 1704,
                1708, 1712, 1716, 1720, 1724, 1728, 1732, 1736, 1740, 1744,
                1748, 1752, 1756, 1760, 1764, 1768, 1772, 1776, 1780, 1784,
                1788, 1792, 1796, 1804, 1808, 1812, 1816, 1820, 1824, 1828,
                1832, 1836, 1840, 1844, 1848, 1852, 1856, 1860, 1864, 1868,
                1872, 1876, 1880, 1884, 1888, 1892, 1896, 1904, 1908, 1912,
                1916, 1920, 1924, 1928, 1932, 1936, 1940, 1944, 1948, 1952,
                1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992,
                1996, 2000, 2004, 2008, 2012, 2016, 2020, 2024, 2028, 2032,
                2036, 2040, 2044, 2048, 2052, 2056, 2060, 2064, 2068, 2072,
                2076, 2080, 2084, 2088, 2092, 2096, 2104, 2108, 2112, 2116,
                2120, 2124, 2128, 2132, 2136, 2140, 2144, 2148, 2152, 2156,
                2160, 2164, 2168, 2172, 2176, 2180, 2184, 2188, 2192, 2196,
                2204, 2208, 2212, 2216, 2220, 2224, 2228, 2232, 2236, 2240,
                2244, 2248, 2252, 2256, 2260, 2264, 2268, 2272, 2276, 2280,
                2284, 2288, 2292, 2296)
out_var <- c("max", "p05", "min", "p95", "max", "p05", "min", "p95")
## read the files
period <- 0
for (EXP in experiments) {
  period <- period + 1

  # if saved extreme indice already exist, skip
  if (period == 1) {
    # if the binary file exist use it as it is faster
    if (file.exists(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             casefold(model), "_", gcm, "_",
             experiments[1], "_", as.character(period_1[1]),
             "-", as.character(period_1[2]), "_", "max.bin"))) {
      # load data and next
      for (i in c("max", "p05", "min", "p95")) {
        f_write <- file(
          paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
                 casefold(model), "_", gcm, "_",
                 experiments[1], "_", as.character(period_1[1]), "-",
                 as.character(period_1[2]), "_", i, ".bin"),
          open = "rb")
        # loop through the 30-years
        for (j in 1:30) {
          outputs[[which(i == c("max", "p05", "min", "p95"))]][j, ] <-
            readBin(f_write, size = 4, what = "numeric", n = 259200,
                    endian = "little")
        }
        close(f_write)
        outputs[[which(i == c("max", "p05", "min", "p95"))]] <-
          apply(outputs[[which(i == c("max", "p05", "min", "p95"))]], 2,
                function(x) {
                  ifelse(is.finite(x), x, NA)
                }
                )
      }
      rm(i) ; rm(j) ; rm(f_write)
      next
    }
  } else {
    # this is for the second periods, similarly
    # check if the bin file exist as it is faster to read
    if (file.exists(
      paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
             casefold(model), "_", gcm, "_", experiments[2], "_",
             as.character(period_2[1]), "-", as.character(period_2[2]),
             "_", "max.bin"))) {
      # load data and next
      for (i in c("max", "p05", "min", "p95")) {
        f_write <- file(
          paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
                 casefold(model), "_", gcm, "_",
                 experiments[2], "_", as.character(period_2[1]),
                 "-", as.character(period_2[2]), "_", i, ".bin"),
          open = "rb")
        for (j in 1:30) {
          outputs[[(4 + which(i == c("max", "p05", "min", "p95")))]][j, ] <-
            readBin(f_write, size = 4, what = "numeric", n = 259200,
                    endian = "little")
        }
        close(f_write)
        outputs[[(4 + which(i == c("max", "p05", "min", "p95")))]] <-
          apply(outputs[[(4 + which(i == c("max", "p05", "min", "p95")))]], 2,
                function(x) {
                  ifelse(is.finite(x), x, NA)
                }
                )
      }
      rm(i) ; rm(j) ; rm(f_write)
      next
    }
  }

  cc    <- 0
  pp    <- get(paste0("period_", period))
  temp  <- matrix(data = 0, nrow = 30, ncol = 259200)
  if (EXP == "picontrol") {
    soc <- "1860soc_co2"
  } else {
    soc <- "2005soc_co2"
  }
  for (years in seq(pp[1], pp[2])) {

    cc <- cc + 1
    days <- seq(as.Date(paste(years, "/1/1", sep = "")),
                as.Date(paste(years, "/12/31", sep = "")), "day")

    ## Derive the folder from experiment and year being processed
    yrs_seq <- seq(1661, 2301, 10)
    if (years <= 1860) {
      folder <- "pre-industrial"
    } else if (years >= 1861 & years <= 2005) {
      folder <- "historical"
    } else if (years >= 2006 & years <= 2099) {
      folder <- "future"
    } else {
      folder <- "future_extended"
    }

    if (years >= 2006 & years <= 2010) {
      # Special case 1
      per_str <- "2006_2010"
      ini <- 2006
    } else if (years >= 2000 & years <= 2005) {
      # Special case 2
      per_str <- "2001_2005"
      ini <- 2001
    } else if (years >= 2091 & years <= 2099) {
      per_str <- "2091_2099"
      ini <- 2091

    } else if (years >= 2291 & years <= 2299) {
      per_str <- "2291_2299"
      ini <- 2291
    } else {
      # Other
      if (years %in% seq(1661, 2291, by = 10)) {
        per_str <- paste0(years, "_", years + 9)
        ini <- years
      } else {
        per_str <- paste0(yrs_seq[yrs_seq >= years][1] - 10,
                          "_", yrs_seq[yrs_seq >= years][1] - 1)
        ini <- yrs_seq[yrs_seq >= years][1] - 10
      }
    }

    # finaly read the file
    file_1 <- nc_open(
      paste0("/data01/julien/models/GW_ISIMIP/", model, "/", gcm, "/",
             folder, "/", casefold(model), "_", gcm, "_ewembi_",
             experiments[period], "_", soc, "_dis_global_daily_",
             per_str, ".nc4"))

    # get the correct year
    if (years == ini) {
      start_y <- 1
    } else {
      start_y <- seq(as.Date(paste0(ini, "/1/1")),
                     as.Date(paste((years - 1), "/12/31", sep = "")),
                     "day")
      if (model %in% c("CLM45", "lpjml")) {
        start_y <- start_y[!(format(start_y, "%m") == "02" &
                               format(start_y, "%d") == "29")]
      }
    }

    # start:
    if (model %in% c("CLM45", "lpjml") & years %in% leap_years) {
      adj <- -1
    } else {
      adj <- 0
    }
    discharge <- ncvar_get(file_1, "dis",
                           start = c(1, 1, length(start_y)),
                           count = c(-1, -1, (length(days) + adj)))

    nc_close(file_1)
    # 4 indices (the order of the data is correct! no need to transpose here)
    if (period == 1) {
      outputs$max_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      max(x, na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
      outputs$p05_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      quantile(x, 0.95,
                                                               na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
      outputs$min_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      min(x, na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
      outputs$p95_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      quantile(x, 0.05,
                                                               na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
    } else {
      outputs$max_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      max(x, na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
      outputs$p05_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      quantile(x, 0.95,
                                                               na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
      outputs$min_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      min(x, na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
      outputs$p95_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2),
                                                  function(x) {
                                                    if (sum(is.na(x)) > 300) {
                                                      NA
                                                    } else {
                                                      quantile(x, 0.05,
                                                               na.rm = TRUE)
                                                    }
                                                  }
                                                  ))
    }
    # carrefull ! the matrix are 720x360 (swapped)
  } # year loop
  # save indice to binary files in order to save time in the future
  if (period == 1) {
    for (i in c("max", "p05", "min", "p95")) {
      f_write <- file(
        paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
               casefold(model), "_", gcm, "_",
               experiments[1], "_", as.character(period_1[1]),
               "-", as.character(period_1[2]), "_", i, ".bin"),
        open = "wb")
      writeBin(as.vector(t(outputs[[which(i == c("max", "p05",
                                                 "min", "p95"))]])),
               f_write, size = 4, endian = "little")
      close(f_write)
    }
  } else { # the period 2
    for (i in c("max", "p05", "min", "p95")) {
      f_write <- file(
        paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/",
               casefold(model), "_", gcm, "_",
               experiments[2], "_", as.character(period_2[1]),
               "-", as.character(period_2[2]), "_", i, ".bin"),
        open = "wb")
      writeBin(as.vector(t(outputs[[(4 + which(i == c("max", "p05",
                                                      "min", "p95")))]])),
               f_write, size = 4, endian = "little")
      close(f_write)
    }
  }
}  # experiment loop

##### part 2: inference and fitting #####
## the 30 extrema data for both sets are loaded
for (i in 1:4) {
  to_save <- sapply(seq_len(ncol(outputs[[i]])),
                    function(j) {
                      if (sum(is.na(outputs[[i]][, j])) > 25 |
                            sum(is.na(outputs[[(i + 4)]][, j])) > 25) {
                        NA
                      } else {
                        unname(ks.test(outputs[[i]][, j],
                                       outputs[[(i + 4)]][, j])$p.value)
                      }
                    }
                    )
  to_save2 <- sapply(seq_len(ncol(outputs[[i]])),
                     function(j) {
                       if (sum(is.na(outputs[[i]][, j])) > 25 |
                              sum(is.na(outputs[[(i + 4)]][, j])) > 25) {
                         NA
                       } else {
                         unname(ks.test(outputs[[i]][, j],
                                        outputs[[(i + 4)]][, j])$statistic)
                       }
                     }
                     )
  # pool all data together, randomly split 1000 time and get the t interval

  if (permutation == "yes") {
    pooled <- rbind(outputs[[i]], outputs[[i + 4]])
    estimates <- matrix(nrow = 100, ncol = 259200)
    for (k in 1:1000) {
      picked <- sample(seq(1:60), 30)
      # save it
      estimates[k, ] <- sapply(seq_len(ncol(pooled)),
                               function(j) {
                                 if (sum(is.na(pooled[picked, j])) > 25 |
                                       sum(is.na(pooled[-picked, j])) > 25) {
                                   NA
                                 } else {
                                   unname(ks.test(pooled[picked, j],
                                                  pooled[-picked, j])$statistic)
                                 }
                               }
                               )
    }
    to_sa$to_save <- as_vector(to_save)
    to_sa$stat <- to_save2
    to_sa$low <- apply(estimates, 2,
                       function(x) quantile(x, 0.025, na.rm = TRUE))
    to_sa$hgh <- apply(estimates, 2,
                       function(x) quantile(x, 0.975, na.rm = TRUE))
  } else {
    to_sa <- as.data.frame(to_save)
    to_sa$stat <- to_save2
  }
  # save the 100 stat ! (estimate)
  file_1 <- file(
    paste0("/data01/julien/projects/extreme_trans_stab/OUT/pi_process/",
           casefold(model), "_", gcm, "_",
           experiments[1], "_", as.character(period_1[1]), "-",
           as.character(period_1[2]), "_",
           experiments[2], "_", as.character(period_2[1]), "-",
           as.character(period_2[2]), "_",
           out_var[i], ".bin"), open = "wb")
  writeBin(as.vector(t(estimates)), File_1, size = 4, endian = "little")
  close(file_1) ; rm(file_1)
  # save
  write_csv(to_sa,
            paste0("/data01/julien/projects/extreme_trans_stab/OUT/",
                   casefold(model), "_", gcm, "_",
                   experiments[1], "_", as.character(period_1[1]),
                   "-", as.character(period_1[2]), "_",
                   experiments[2], "_", as.character(period_2[1]),
                   "-", as.character(period_2[2]), "_",
                   out_var[i], ".csv"))
} # indice loop
rm(to_save) ; rm(pooled) ; rm(estimates) ; rm(to_save2)


#### Gumbel fitting
#if (experiments[1] == "rcp26") {
  for (i in 1:8) {
    if (i %in% c(1, 2, 5, 6)) {
      outputs[[i]] <- apply(outputs[[i]], 2,
                            function(x) sort(x, na.last = TRUE))
    } else {
      outputs[[i]] <- outputs[[i]] * -1
      outputs[[i]] <- apply(outputs[[i]], 2,
                            function(x) sort(x, decreasing = TRUE,
                                             na.last = TRUE))
    }
    # fit L-moment
    m1 <- apply(outputs[[i]], 2,
                function(x) if (sum(is.na(x)) >= 25) {
                  NA
                } else {
                  mean(x, na.rm = TRUE)
                }
                )
    m2 <- apply(outputs[[i]], 2,
                function(x) if (sum(is.na(x)) >= 25) {
                  NA
                } else {
                  1 / 30 * sum((seq(1, 30, 1) - 1) / (30 - 1) * x)
                }
                )
    l1 <- m1
    l2 <- 2 * m2 - m1
    #
    alpha   <- l2 / log(2)
    epsilon <- l1 - alpha * 0.57721
    y100_ext   <- epsilon - alpha * log(-log(1 - 1 / 100))
    # save
    if (i < 5) {
      write_csv(data.frame(alpha = alpha, espilon = epsilon, y100 = Y100_ext),
                paste0("/data01/julien/projects/extreme_trans_stab/OUT/gumbel/",
                       casefold(model), "_", gcm, "_",
                       experiments[1], "_", as.character(period_1[1]),
                       "-", as.character(period_1[2]), "_",
                     out_var[i], ".csv"))
    } else {
      write_csv(data.frame(alpha = alpha, espilon = epsilon, y100 = Y100_ext),
                paste0("/data01/julien/projects/extreme_trans_stab/OUT/gumbel/",
                       casefold(model), "_", gcm, "_",
                       experiments[2], "_", as.character(period_2[1]),
                       "-", as.character(period_2[2]), "_",
                       out_var[i], ".csv"))
    }
  }
#}
