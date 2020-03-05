#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

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
# 7) 

## CC inputs
if (length(args) != 8) {
  stop("Wrong nb of argument:\n
           MODEL -> CLM45 / cwatm
           GCM   -> hadgem2-es ipsl-cm5a-lr miroc5
           EXP1  -> hist / picontrol / rcp26 / rcp60
           EXP2 -> hist / picontrol / rcp26 / rcp60
           period1 -> YYYY YYYY
           period2  -> YYYY YYYY")
} else {
  MODEL       <- args[1]                           # Other options:
  GCM         <- args[2]                           # possible choice: hadgem2-es, 
  EXPERIMENTS <- c(args[3], args[4])               # The 2 experiments: hist, picontrol, rcp26, rcp85
  PERIOD_1    <- as.integer(c(args[5], args[6]))   # 30 year period
  PERIOD_2    <- as.integer(c(args[7], args[8]))   # 30 year period
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
out_var <- c("max", "p05", "min", "p95", "max", "p05", "min", "p95")
## read the files
PERIOD <- 0
for (EXP in EXPERIMENTS) {
  PERIOD <- PERIOD + 1
  #count <- 1
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
    }

# 
    # start :   ; 
    discharge <- ncvar_get(file_1, "dis", start = c(1, 1, length(start_y)), count = c(-1, -1, length(days)))

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
}  # experiment loop
rm(cc) ; rm(temp) ; rm(PP) ; rm(PERIOD)

##### part 2: inference and fitting #####    ## the 30 extrema data for both sets are loaded
for (i in 1:4) {
  to_save <- sapply(1:ncol(outputs[[i]]), function(j) { if ( sum(is.na(outputs[[i]][,j])) > 25 | sum(is.na(outputs[[(i+4)]][,j])) > 25) {NA} else {ks.test(outputs[[i]][,j], outputs[[(i+4)]][,j])$p.value} })
  # save
  write_csv(as.data.frame(to_save), paste0("/data01/julien/projects/extreme_trans_stab/OUT/", casefold(MODEL), "_", GCM, "_",
                 EXPERIMENTS[1], "_", as.character(PERIOD_1[1]), "-", as.character(PERIOD_1[2]), "_",
                 EXPERIMENTS[2], "_", as.character(PERIOD_2[1]), "-", as.character(PERIOD_2[2]), "_",
                 out_var[i], ".csv"))
}
rm(to_save)


#### Gumbel fitting
for (i in 1:8) {
  outputs[[i]] <- apply(outputs[[i]], 2, function(x) sort(x, na.last = TRUE))
  # fit L-moment
  M1 <- apply(outputs[[i]], 2, function(x) if (sum(is.na(x)) >= 25) {NA} else {mean(x, na.rm = TRUE)})
  M2 <- apply(outputs[[i]], 2, function(x) if (sum(is.na(x)) >= 25) {NA} else {1 / 30 * sum((seq(1, 30, 1) - 1) / (30 - 1) * x)})
  L1 <- M1
  L2 <- 2 * M2 - M1
  
  alpha   <- L2 / log(2)
  epsilon <- L1 - alpha * 0.57721
  Y100_ext   <- epsilon - alpha * log(-log(1 - 1/100))
  # save
  if (i < 5) {
    write_csv(data.frame(alpha = alpha,espilon = epsilon,y100 = Y100_ext),
              paste0("/data01/julien/projects/extreme_trans_stab/OUT/gumbel/", casefold(MODEL), "_", GCM, "_",
                     EXPERIMENTS[1], "_", as.character(PERIOD_1[1]), "-", as.character(PERIOD_1[2]), "_",
                     out_var[i], ".csv"))      
  } else {
      write_csv(data.frame(alpha = alpha,espilon = epsilon,y100 = Y100_ext),
            paste0("/data01/julien/projects/extreme_trans_stab/OUT/gumbel/", casefold(MODEL), "_", GCM, "_",
                   EXPERIMENTS[2], "_", as.character(PERIOD_2[1]), "-", as.character(PERIOD_2[2]), "_",
                   out_var[i], ".csv"))  
  }

}
