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
              nxl = riv_seq)
rm(riv_number) ; rm(riv_nxl) ; rm(riv_seq)

## jobs
# flag cell where n_neg and n_neg2 >= 6
all_data$flag1 <- ifelse(all_data$n_neg >= 6, 1, 0)
all_data$flag2 <- ifelse(all_data$n_neg2 >= 6, 1, 0)

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
