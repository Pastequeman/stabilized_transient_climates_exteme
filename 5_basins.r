#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

## Purpose:


## Libraries:
library("readr")
suppressMessages(library("dplyr"))

## user input
#none

## read data
all_data <- read_csv("./../OUT/process/grid_cells_gcms_all_ind.csv")

# rivers characteristics
# number
f_riv <- file("./../IN/riv_num_/rivnum.WATCH.hlf", open = "rb")
riv_number <- readBin(f_riv, what = "numeric", n = 259200, size = 4, endian = "big")
close(f_riv)

# seq (river sequence)
f_riv <- file("./../IN/riv_seq_/rivseq.WATCH.hlf", open = "rb")
riv_seq <- readBin(f_riv, what = "numeric", n = 259200, size = 4, endian = "big")
close(f_riv)

#nxl (distance to next L)
f_riv <- file("./../IN/riv_nxl_/rivnxl.WATCH.hlf", open = "rb")
riv_nxl <- readBin(f_riv, what = "numeric", n = 259200, size = 4, endian = "big")
close(f_riv)

## clean
rm(f_riv)
riv <- tibble(L   = seq(1, 259200),
              lon = rep(seq(-179.75, 179.75, 0.5), times = 360),
              lat = rep(seq(89.75, -89.75, -0.5), each = 720),
              num = riv_number,
              seq = riv_seq,
              nxl = riv_nxl)
rm(riv_number) ; rm(riv_nxl) ; rm(riv_seq)

## jobs
# flag cell where n_neg and n_neg2 >= 6
all_data$flag1 <- ifelse(all_data$n_neg == 6 & all_data$gcm != "hadgem2-es" |
                        all_data$n_neg == 7 & all_data$gcm == "hadgem2-es" , 1, 0)
all_data$flag2 <- ifelse(all_data$n_neg2 == 6 & all_data$gcm != "hadgem2-es" |
                        all_data$n_neg2 == 7 & all_data$gcm == "hadgem2-es" , 1, 0)

ind1 <- all_data %>% filter(flag1 == 1) %>% select(id, indice, gcm)
ind2 <- all_data %>% filter(flag2 == 1) %>% select(id, indice, gcm)

ind1 <- ind1 %>% left_join(riv, by = c("id" = "L"))
ind2 <- ind2 %>% left_join(riv, by = c("id" = "L"))

##
t1 <- ind1 %>% group_by(indice, gcm, num) %>% summarise(n = n())
t2 <- ind2 %>% group_by(indice, gcm, num) %>% summarise(n = n())

## keep only basin/river with at least 40 grid-cell marked
keep1 <- t1 %>% filter(n >= 30)
keep2 <- t2 %>% filter(n >= 30)
# example how to keep relevant data
#keep1 %>% ungroup() %>% filter(indice == "max") %>% select(n)
#keep2 %>% ungroup() %>% filter(indice == "max") %>% select(n)

# now identify the L of the slected river which is the most downstream?

final_loc1 <- riv[riv$num %in% keep1$num & riv$L == riv$nxl,]
final_loc2 <- riv[riv$num %in% keep2$num & riv$L == riv$nxl,]

## test: keep1
for (ind in c("max", "p05", "min", "p95")) {
  a <- ind1 %>% filter(num %in% keep1[keep1$indice == ind, ]$num, indice == ind)
  for (it in unique(a$num)) {
    print(c(it, max(a[a$num == it, ]$nxl)))
  }
}

## assess natural variability

## for these location, read the indice files for all models / GCM and 2 period
for (m in c("cwatm_", "matsiro_", "clm45_", "lpjml_", "pcr-globwb_", "watergap2_", "h08_")) { # add more as they come
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
      file_in1 <- file(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", m, gcms, "rcp26_",
                              period1, "_", id, ".bin"), open = "rb")
      file_in2 <- file(paste0("/data01/julien/projects/extreme_trans_stab/OUT/indice/", m, gcms, "rcp26_",
                              period2, "_", id, ".bin"), open = "rb")
      for (y in 1:30) {
        temp <- readBin(file_in1, what = "numeric", n = 259200, size = 4, endian = "little")
        temp <- temp - readBin(file_in2, what = "numeric", n = 259200, size = 4, endian = "little")
        
        final_loc1$flow  <- temp[final_loc1$L]
        final_loc1$year  <- y
        final_loc1$ind   <- id
        final_loc1$gcm   <- gcms
        final_loc1$model <- m
        if (gcms == "hadgem2-es_" & m == "cwatm_" & id == "max" & y == 1) {
          all_data <- final_loc1
        } else {
          all_data <- rbind(all_data, final_loc1)
        }
      }
      close(file_in1) ; close(file_in2)
    }
  }
}
# clean
rm(temp) ; rm(file_in1) ; rm(file_in2)

## process the data to keep gcm, L, indice...
gcm_data <- 
all_data %>% group_by(L, num, seq, ind, gcm, year) %>% summarise(ave_flow = mean(flow, na.rm = TRUE),
                                                                 low_flow = quantile(flow, 0.025, na.rm = TRUE),
                                                                 hgt_flow = quantile(flow, 0.025, na.rm = TRUE))

## test
g <- ggplot()
g + geom_line(data = gcm_data %>% filter(L == 156673, ind == "min"), aes(x = year, y = ave_flow, color = gcm)) + 
    geom_line(data = gcm_data %>% filter(L == 156673, ind == "p95"), aes(x = year, y = ave_flow, color = gcm))
