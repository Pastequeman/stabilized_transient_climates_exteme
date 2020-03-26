#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

## Purpose:
# for two 30-year period, compute extreme indices
# Thene instead of checking all series
# check 95th percentile (influence on gumbel fitting)
# at each grid cells
# use directly the netcdf files
#
# 1) Read required files for first period
# 2) Compute extreme indices for first period
# 3) Read file for second period
# 4) Compute extreme indices for second period
# 5) perform t test > write out

## CC inputs
if (length(args) != 9) {
  stop("Wrong nb of argument:\n
           MODEL -> CLM45 / cwatm
           GCM   -> hadgem2-es ipsl-cm5a-lr miroc5
           EXP1  -> hist / picontrol / rcp26 / rcp60
           EXP2 -> hist / picontrol / rcp26 / rcp60
           period1 -> YYYY YYYY
           period2  -> YYYY YYYY
           PERMUTATION -> yes/no")
} else {
  MODEL       <- args[1]                           # Other options:
  GCM         <- args[2]                           # possible choice: hadgem2-es, 
  EXPERIMENTS <- c(args[3], args[4])               # The 2 experiments: hist, picontrol, rcp26, rcp85
  PERIOD_1    <- as.integer(c(args[5], args[6]))   # 30 year period
  PERIOD_2    <- as.integer(c(args[7], args[8]))   # 30 year period
  PERMUTATION <- args[9]
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
leap_years <- c(1664, 1668, 1672, 1676, 1680, 1684, 1688, 1692, 1696, 1704, 1708, 1712, 1716, 1720, 1724, 1728, 1732, 1736, 1740, 1744, 1748, 1752, 1756, 1760, 1764, 1768, 1772, 1776, 1780, 1784, 1788, 1792, 1796, 1804, 1808, 1812, 1816, 1820, 1824, 1828, 1832, 1836, 1840, 1844, 1848, 1852, 1856, 1860, 1864, 1868, 1872, 1876, 1880, 1884, 1888, 1892, 1896, 1904, 1908, 1912, 1916, 1920, 1924, 1928, 1932, 1936, 1940, 1944, 1948, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2024, 2028, 2032, 2036, 2040, 2044, 2048, 2052, 2056, 2060, 2064, 2068, 2072, 2076, 2080, 2084, 2088, 2092, 2096, 2104, 2108, 2112, 2116, 2120, 2124, 2128, 2132, 2136, 2140, 2144, 2148, 2152, 2156, 2160, 2164, 2168, 2172, 2176, 2180, 2184, 2188, 2192, 2196, 2204, 2208, 2212, 2216, 2220, 2224, 2228, 2232, 2236, 2240, 2244, 2248, 2252, 2256, 2260, 2264, 2268, 2272, 2276, 2280, 2284, 2288, 2292, 2296)

out_var <- c("max", "p05", "min", "p95", "max", "p05", "min", "p95")
## read the files
PERIOD <- 0

for (EXP in EXPERIMENTS) {
  PERIOD <- PERIOD + 1
  #count <- 1
  # if bin file exist just load it
  if (PERIOD == 1) {
    if (file.exists(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", casefold(MODEL), "_", GCM, "_",
                           EXPERIMENTS[1], "_", as.character(PERIOD_1[1]), "-", as.character(PERIOD_1[2]), "_", "max.bin"))) {
      # load data and next
      for (i in c("max", "p05", "min", "p95")) {
        f_write <- file(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", casefold(MODEL), "_", GCM, "_",
                               EXPERIMENTS[1], "_", as.character(PERIOD_1[1]), "-", as.character(PERIOD_1[2]), "_", i, ".bin"), open = "rb")
        for (j in 1:30) {
          outputs[[which(i == c("max", "p05", "min", "p95"))]][j,] <- readBin(f_write, size = 4, what = "numeric", n = 259200, endian = "little")
        }
        
        close(f_write)
        outputs[[which(i == c("max", "p05", "min", "p95"))]] <- apply(outputs[[which(i == c("max", "p05", "min", "p95"))]], 2, function(x) {ifelse(is.finite(x),x , NA)})
      }
      rm(i) ; rm(j) ; rm(f_write)
      next
    } 
  } else {
    if (file.exists(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", casefold(MODEL), "_", GCM, "_",
                           EXPERIMENTS[2], "_", as.character(PERIOD_2[1]), "-", as.character(PERIOD_2[2]), "_", "max.bin"))) {
      # load data and next
      for (i in c("max", "p05", "min", "p95")) {
        f_write <- file(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", casefold(MODEL), "_", GCM, "_",
                               EXPERIMENTS[2], "_", as.character(PERIOD_2[1]), "-", as.character(PERIOD_2[2]), "_", i, ".bin"), open = "rb")
        for (j in 1:30) {
          outputs[[(4+which(i == c("max", "p05", "min", "p95")))]][j,] <- readBin(f_write, size = 4, what = "numeric", n = 259200, endian = "little")
        }
        
        close(f_write)
        outputs[[(4+which(i == c("max", "p05", "min", "p95")))]] <- apply(outputs[[(4+which(i == c("max", "p05", "min", "p95")))]], 2, function(x) {ifelse(is.finite(x),x , NA)})
      }
      rm(i) ; rm(j) ; rm(f_write)
      next
    }    
  }
  
  cc    <- 0
  PP    <- get(paste0("PERIOD_", PERIOD))
  temp  <- matrix(data = 0, nrow = 30, ncol = 259200)
  if (EXP == "picontrol") {
    soc <- "1860soc_co2"
  } else {
    soc <- "2005soc_co2"
  }
  for (years in seq(PP[1], PP[2])) {
    #print(years) ## DEBUG
    cc <- cc + 1
    days <- seq(as.Date(paste(years, "/1/1", sep = "")), as.Date(paste(years, "/12/31", sep = "")), "day")

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
        per_str <- paste0(years, "_", years+9)
        ini <- years
      } else {
        per_str <- paste0(yrs_seq[yrs_seq >= years][1] - 10, "_", yrs_seq[yrs_seq >= years][1] - 1)
        ini <- yrs_seq[yrs_seq >= years][1] - 10
      }
    }

    # finaly read the file
    file_1 <- nc_open(paste0("/data01/julien/models/GW_ISIMIP/", MODEL, "/", GCM, "/", folder, "/" ,casefold(MODEL), "_", GCM, "_ewembi_",EXPERIMENTS[PERIOD], "_", soc, "_dis_global_daily_", per_str, ".nc4"))

    # get the correct year
    if (years == ini) {
      start_y <- 1
    } else {
      start_y <- seq(as.Date(paste0(ini, "/1/1")), as.Date(paste((years - 1), "/12/31", sep = "")), "day")
      if (MODEL %in% c("CLM45", "lpjml")) {
        start_y <- start_y[!(format(start_y,"%m") == "02" & format(start_y, "%d") == "29")]
      }
    }

# 
    # start :   ;
    if (MODEL %in% c("CLM45", "lpjml") & years %in% leap_years) {
      adj <- -1
    } else {
      adj <- 0
    }
    discharge <- ncvar_get(file_1, "dis", start = c(1, 1, length(start_y)), count = c(-1, -1, (length(days)+adj)))

    nc_close(file_1)
    # 4 indices (the order of the data is correct! no need to transpose here)
    if (PERIOD == 1) {
      
      outputs$max_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {max(x, na.rm = TRUE)}}))
      outputs$p05_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {quantile(x, 0.95 ,na.rm = TRUE)}}))
      outputs$min_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {min(x, na.rm = TRUE)}}))
      outputs$p95_dis_v1[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {quantile(x, 0.05, na.rm = TRUE)}}))
    } else {
      outputs$max_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {max(x, na.rm = TRUE)}}))
      outputs$p05_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {quantile(x, 0.95 ,na.rm = TRUE)}}))
      outputs$min_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {min(x, na.rm = TRUE)}}))
      outputs$p95_dis_v2[cc, ] <- as.vector(apply(discharge, c(1, 2), function(x) {if (sum(is.na(x)) > 300) {NA} else {quantile(x, 0.05, na.rm = TRUE)}}))
    }
    # carrefull ! the matrix are 720x360 (swapped)
  } # year loop
  # save indice to save time
  if (PERIOD == 1) {
    for (i in c("max", "p05", "min", "p95")) {
      f_write <- file(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", casefold(MODEL), "_", GCM, "_",
                             EXPERIMENTS[1], "_", as.character(PERIOD_1[1]), "-", as.character(PERIOD_1[2]), "_", i, ".bin"), open = "wb")
      writeBin(as.vector(t(outputs[[which(i == c("max", "p05", "min", "p95"))]])), f_write, size = 4, endian = "little")
      close(f_write)
    }
  } else {
    for (i in c("max", "p05", "min", "p95")) {
      f_write <- file(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", casefold(MODEL), "_", GCM, "_",
                             EXPERIMENTS[2], "_", as.character(PERIOD_2[1]), "-", as.character(PERIOD_2[2]), "_", i, ".bin"), open = "wb")
      writeBin(as.vector(t(outputs[[(4+which(i == c("max", "p05", "min", "p95")))]])), f_write, size = 4, endian = "little")
      close(f_write)
    }    
  }
  # end save bin files
}  # experiment loop
#rm(cc) ; rm(temp) ; rm(PP) ; rm(PERIOD)
#rm(years) ; rm(leap_years) ; rm(yrs_seq) ; rm(ini)
#rm(folder) ; rm(discharge) ; rm(days) ; rm(adj) ; rm(per_str) ; rm(soc)
##### part 2: inference and fitting #####    ## the 30 extrema data for both sets are loaded
for (i in 1:4) {
  if (i %in% c(1, 2)) {
    target <- 0.95
  } else {
    target <- 0.05
  }
  to_save <- sapply(1:ncol(outputs[[i]]), function(j) { if ( sum(is.na(outputs[[i]][,j])) > 25 | sum(!is.finite(outputs[[i]][,j])) > 25 | sum(is.na(outputs[[(i+4)]][,j])) > 25 | sum(!is.finite(outputs[[(i+4)]][,j])) > 25) {NA} else {quantile(outputs[[i]][,j], target, na.rm = TRUE) - quantile(outputs[[(i+4)]][,j], target, na.rm = TRUE)} })

  if (PERMUTATION == "yes") {
    # pool all data together, randomly split 1000 time and get the t interval
    pooled <- rbind(outputs[[i]], outputs[[i+4]])
    estimates <- matrix(nrow = 100, ncol = 259200)
    for (k in 1:100) {
      #print(k) ## DEBUG
      picked <- sample(seq(1:60), 30)
      #save it
      estimates[k,] <- sapply(1:ncol(pooled), function(j) { if ( sum(is.na(pooled[picked,j])) > 25 | sum(!is.finite(pooled[picked,j])) > 25 | sum(is.na(pooled[-picked,j])) > 25 | sum(!is.finite(pooled[-picked,j]))) {NA} else {quantile(pooled[picked,j], target, na.rm = TRUE) - quantile(pooled[-picked,j], target, na.rm = TRUE)}})  
    }
    to_save <- as.data.frame(to_save)
    to_save$low <- apply(estimates, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
    to_save$hgh <- apply(estimates, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  } else {
    to_save <- as.data.frame(to_save)
  }
  
  # save
  write_csv(to_save, paste0("/data01/julien/projects/extreme_trans_stab/OUT/high_value_", casefold(MODEL), "_", GCM, "_",
                            EXPERIMENTS[1], "_", as.character(PERIOD_1[1]), "-", as.character(PERIOD_1[2]), "_",
                            EXPERIMENTS[2], "_", as.character(PERIOD_2[1]), "-", as.character(PERIOD_2[2]), "_",
                            out_var[i], ".csv"))
} # all indices 
rm(to_save) ; rm(pooled) ; rm(estimates)
