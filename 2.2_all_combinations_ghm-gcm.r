#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

## Purpose:
# investigate all possible combinations of gcm/ghm
# and report consensus
# only <66%

## User inputs
# none

# Libraries
library("tibble")
library("readr")
suppressMessages(library("dplyr"))

# read initial file
all_main_data_filter <-
read_csv("./../IN/all_main_data.csv",
         col_types = cols(to_save = col_double(),
                          stat    = col_double(),
                          low     = col_double(),
                          hgh     = col_double(),
                          sign1   = col_integer(),
                          area    = col_double(),
                          id      = col_integer(),
                          indice  = col_character(),
                          model   = col_character(),
                          gcm     = col_character(),
                          Lat     = col_double(),
                          Lon     = col_double(),
                          Cls     = col_character()
                          )
         )

##### Jobs #####
# all possible permutation
perm <- gtools::permutations(8, 6, c(0, 1, 2, 3, 4, 5, 6, 7),
                             repeats.allowed = TRUE)
colnames(perm) <- c("cwatm", "h08", "lpjml", "matsiro",
                    "pcr-globwb", "watergap2")
model_to_include <- c("cwatm_", "h08_", "lpjml_", "matsiro_",
                      "pcr-globwb_", "watergap2_")
gcm_to_include <- c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_")

# store outputs here
robustness_combination <- tibble(cases   = integer(1),
                                 low_all = double(1),
                                 low_fil = double(1),
                                 hgh_all = double(1),
                                 hgh_fil = double(1),
                                 n_model = integer(1),
                                 n_ghm   = integer(1))

# function taht will report the gcms involve for all ghm
involved_gcm <- function(x) {
  if (x == 0) {
    return(NULL)
  } else if (x == 1) {
    return("hadgem2-es_")
  } else if (x == 2) {
    return("ipsl-cm5a-lr_")
  } else if (x == 3) {
    return("miroc5_")
  } else if (x == 4) {
    return(c("hadgem2-es_", "ipsl-cm5a-lr_"))
  } else if (x == 5) {
    return(c("hadgem2-es_", "miroc5_"))
  } else if (x == 6) {
    return(c("miroc5_", "ipsl-cm5a-lr_"))
  } else if (x == 7) {
    return(c("hadgem2-es_", "ipsl-cm5a-lr_", "miroc5_"))
  }
}

# return the nb of ghm
nb_ghm <- function(x) {
  if (x != 0) {
    return(1)
  } else {
    return(0)
  }
}

# main job and loop
for (i in 2:nrow(perm)) {
  print(i)
  robustness_combination[1, 1] <- i                 # iteration
  robustness_combination[1, 6] <- sum(perm[i, ])    # nb of models
  robustness_combination[1, 7] <- sum(sapply(perm[i, ],
                                             nb_ghm))    # nb of ghm

  # the subsample, according to  perm
  # low streamflow
   submodels_low <-
    all_main_data_filter %>%
    filter(indice == "min",
           model == "cwatm_"      & gcm %in% involved_gcm(perm[i, 1]) |
           model == "h08_"        & gcm %in% involved_gcm(perm[i, 2]) |
           model == "lpjml_"      & gcm %in% involved_gcm(perm[i, 3]) |
           model == "matsiro_"    & gcm %in% involved_gcm(perm[i, 4]) |
           model == "pcr-globwb_" & gcm %in% involved_gcm(perm[i, 5]) |
           model == "watergap2_"  & gcm %in% involved_gcm(perm[i, 6])
           ) %>%
    group_by(id, Cls) %>%
    summarise(n           = n(),
              tot_signi   = sum(sign1, na.rm = TRUE),
              m_area      = mean(area),
              .groups     = "keep")

  # high streamflow
  submodels_hgh <-
    all_main_data_filter %>%
    filter(indice == "max",
           model == "cwatm_"      & gcm %in% involved_gcm(perm[i, 1]) |
           model == "h08_"        & gcm %in% involved_gcm(perm[i, 2]) |
           model == "lpjml_"      & gcm %in% involved_gcm(perm[i, 3]) |
           model == "matsiro_"    & gcm %in% involved_gcm(perm[i, 4]) |
           model == "pcr-globwb_" & gcm %in% involved_gcm(perm[i, 5]) |
           model == "watergap2_"  & gcm %in% involved_gcm(perm[i, 6])
           ) %>%
    group_by(id, Cls) %>%
    summarise(n           = n(),
              tot_signi   = sum(sign1, na.rm = TRUE),
              m_area      = mean(area),
              .groups     = "keep")

  ## evaluate consensus and save it
  # 1
  robustness_combination[1, 2] <-
    submodels_low %>%
    group_by(id, Cls, n, tot_signi) %>%
    tally(name = "tot") %>%
    ungroup() %>%
    mutate(tot_grid = sum(tot, na.rm = TRUE)) %>%
    filter(tot_signi >= 66 / 100 * n) %>%
    summarise(n    = n(),
              frac = sum(n) / mean(tot_grid) * 100) %>%
    select(frac) %>%
    unname()

  # 2
  robustness_combination[1, 3] <-
    submodels_low %>%
    filter(!Cls %in% c("BWh", "BWk")) %>%
    group_by(id, Cls, n, tot_signi) %>%
    tally(name = "tot") %>%
    ungroup() %>%
    mutate(tot_grid = sum(tot, na.rm = TRUE)) %>%
    filter(tot_signi >= 66 / 100 * n) %>%
    summarise(n    = n(),
              frac = sum(n) / mean(tot_grid) * 100) %>%
    select(frac) %>%
    unname()

  # 3
  robustness_combination[1, 4] <-
    submodels_hgh %>%
    group_by(id, Cls, n, tot_signi) %>%
    tally(name = "tot") %>%
    ungroup() %>%
    mutate(tot_grid = sum(tot, na.rm = TRUE)) %>%
    filter(tot_signi >= 66 / 100 * n) %>%
    summarise(n    = n(),
              frac = sum(n) / mean(tot_grid) * 100) %>%
    select(frac) %>%
    unname()

  # 4
  robustness_combination[1, 5] <-
    submodels_hgh %>%
    filter(!Cls %in% c("BWh", "BWk")) %>%
    group_by(id, Cls, n, tot_signi) %>%
    tally(name = "tot") %>%
    ungroup() %>%
    mutate(tot_grid = sum(tot, na.rm = TRUE)) %>%
    filter(tot_signi >= 66 / 100 * n) %>%
    summarise(n    = n(),
              frac = sum(n) / mean(tot_grid) * 100) %>%
    select(frac) %>%
    unname()

  # write
  write_csv(robustness_combination, "./../OUT/all_combinations.csv",
            append = TRUE)
}
