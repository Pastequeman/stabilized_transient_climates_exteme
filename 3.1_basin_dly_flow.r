#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


## Purpose:
# Extract daily flow for specific locations


# user input
if (length(args) != 4) {
  stop("Wrong nb of argument:\n
           MODEL -> CLM45 / cwatm / h08 / watergap2 /
                    pcr-globwb / lpjml / matsiro\n
           gcm   -> hadgem2-es / ipsl-cm5a-lr / miroc5\n
           EXP  -> hist / picontrol / rcp26 / rcp60\n
           status -> transient / stabilized")
} else {
  model       <- args[1]
  gcm         <- args[2]
  experiment  <- args[3]
  status      <- args[4]
}

if (status == "transient") {
  if (gcm == "hadgem2-es") {
   period <- c(2012, 2041)
  } else if (gcm == "ipsl-cm5a-lr") {
    period <- c(2015, 2044)
  } else {
    period <- c(2035, 2064)
  }
} else {
  if (gcm == "hadgem2-es") {
    period <- c(2147, 2176)
  } else if (gcm == "ipsl-cm5a-lr") {
    period <- c(2150, 2179)
  } else {
    period <- c(2243, 2272)
  }
}

# libraries required
library("ncdf4")
library("readr")
suppressMessages(library("dplyr"))

# list of leap year since some models ignore them
leap_years <- c(1664, 1668, 1672, 1676, 1680, 1684,
                1688, 1692, 1696, 1704, 1708, 1712,
                1716, 1720, 1724, 1728, 1732, 1736,
                1740, 1744, 1748, 1752, 1756, 1760,
                1764, 1768, 1772, 1776, 1780, 1784,
                1788, 1792, 1796, 1804, 1808, 1812,
                1816, 1820, 1824, 1828, 1832, 1836,
                1840, 1844, 1848, 1852, 1856, 1860,
                1864, 1868, 1872, 1876, 1880, 1884,
                1888, 1892, 1896, 1904, 1908, 1912,
                1916, 1920, 1924, 1928, 1932, 1936,
                1940, 1944, 1948, 1952, 1956, 1960,
                1964, 1968, 1972, 1976, 1980, 1984,
                1988, 1992, 1996, 2000, 2004, 2008,
                2012, 2016, 2020, 2024, 2028, 2032,
                2036, 2040, 2044, 2048, 2052, 2056,
                2060, 2064, 2068, 2072, 2076, 2080,
                2084, 2088, 2092, 2096, 2104, 2108,
                2112, 2116, 2120, 2124, 2128, 2132,
                2136, 2140, 2144, 2148, 2152, 2156,
                2160, 2164, 2168, 2172, 2176, 2180,
                2184, 2188, 2192, 2196, 2204, 2208,
                2212, 2216, 2220, 2224, 2228, 2232,
                2236, 2240, 2244, 2248, 2252, 2256,
                2260, 2264, 2268, 2272, 2276, 2280,
                2284, 2288, 2292, 2296)

if (experiment == "picontrol") {
    soc <- "1860soc_co2"
} else {
    soc <- "2005soc_co2"
}

for (years in seq(period[1], period[2])) {

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
  file_flow <- nc_open(paste0("/data01/julien/models/GW_ISIMIP/",
                              model, "/", gcm, "/", folder, "/",
                              casefold(model), "_", gcm, "_ewembi_",
                              experiment, "_", soc, "_dis_global_daily_",
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
      days <- days[!(format(days, "%m") == "02" &
                       format(days, "%d") == "29")]
    }
  }
   # start :   ;
  if (model %in% c("CLM45", "lpjml") & years %in% leap_years) {
    adj <- -1
  } else {
    adj <- 0
  }
   # one year of discharge, globally
  discharge <- ncvar_get(file_flow, "dis",
                         start = c(1, 1, length(start_y)),
                         count = c(-1, -1, (length(days))))
  nc_close(file_flow)

  basins_id <- 129140 #bas$id
  format_0.5 <- tibble(L   = seq(1, 259200),
                       lon = rep(seq(1, by = 1, length.out = 720), times = 360),
                       lat = rep(seq(1, by = 1, length.out = 360), each = 720))
   #
  for (b in 1:length(basins_id)) {
    # dim: 365 days
    temp <- discharge[format_0.5[basins_id[b], ]$lon,
                      format_0.5[basins_id[b], ]$lat, ]
    to_save <- tibble(id     = basins_id[b],
                      years  = years,
                      date   = days,
                      status = status,
                      flow   = temp)
    write_csv(to_save,
              paste0("../OUT/daily_flow/", basins_id[b],
                     "_", model, "_", gcm, "_", experiment,
                     "_daily_flow.csv"), append = TRUE)
  }
} # year loop
