## Dependancy
library(tidyverse)
library(ggplot2)
library(ggdensity)
library(scales)
library(patchwork)
library(ggpubr)
library(readr)
library(PRROC)
source("utilities/FONCTIONS.R")

# DATA
## OBS
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)
survey_obs <- read.csv("data/survey_obs.csv", header = TRUE, row.names = 1)
traits <- read.csv("data/traits.csv", header = TRUE, row.names = 1)

## MEM
# Load precomputed community indices from MEM and Proxy
CM_mem <- read.csv("data/community_indices/CM_mem.csv", header = TRUE, row.names = 1)
colnames(CM_mem) <- paste0(colnames(CM_mem), "_mem")
CM_mem$ID <- rownames(CM_mem)

CSTD_mem <- read.csv("data/community_indices/CSTD_mem.csv", header = TRUE, row.names = 1)
colnames(CSTD_mem) <- paste0(colnames(CSTD_mem), "_mem")
CSTD_mem$ID <- rownames(CSTD_mem)

SR_mem <- read.csv("data/community_indices/SR_mem.csv", header = TRUE, row.names = 1)
colnames(SR_mem) <- paste0(colnames(SR_mem), "_mem")
SR_mem$ID <- rownames(SR_mem)


## SSDM
sdm_raw_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_raw_", i, ".csv"), header = TRUE, row.names = 1))
sdm_raw <- do.call(rbind, sdm_raw_list)
sdm_raw <- sdm_raw[rownames(survey_obs), ]

## Degradations of SSDM :
degraded_predictions <- function(df, w = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df_noisy <- t(apply(df, 1, function(probs) {
    shuffled <- sample(probs)
    (1 - w) * probs + w * shuffled
  }))
  df_noisy <- as.data.frame(df_noisy)
  colnames(df_noisy) <- colnames(df)
  rownames(df_noisy) <- rownames(df)
  
  return(df_noisy)
}
#sdm_degrad_1 <- degraded_predictions(sdm_raw, w = 1)
#write.csv(sdm_degrad_1, file="data/sdm_degrad_1.csv")
sdm_degrad_1 <- readr::read_csv("data/sdm_degrad_1.csv") |>
  tibble::column_to_rownames(var = "...1")

#sdm_degrad_0.8 <- degraded_predictions(sdm_raw, w = 0.8)
#write.csv(sdm_degrad_0.8, file="data/sdm_degrad_0.8.csv")
sdm_degrad_0.8 <- readr::read_csv("data/sdm_degrad_0.8.csv") |>
  tibble::column_to_rownames(var = "...1")

#sdm_degrad_0.6 <- degraded_predictions(sdm_raw, w = 0.6)
#write.csv(sdm_degrad_0.6, file="data/sdm_degrad_0.6.csv")
sdm_degrad_0.6 <- readr::read_csv("data/sdm_degrad_0.6.csv") |>
  tibble::column_to_rownames(var = "...1")

#sdm_degrad_0.4 <- degraded_predictions(sdm_raw, w = 0.4)
#write.csv(sdm_degrad_0.4, file="data/sdm_degrad_0.4.csv")
sdm_degrad_0.4 <- readr::read_csv("data/sdm_degrad_0.4.csv") |>
  tibble::column_to_rownames(var = "...1")

#sdm_degrad_0.2 <- degraded_predictions(sdm_raw, w = 0.2)
#write.csv(sdm_degrad_0.2, file="data/sdm_degrad_0.2.csv")
sdm_degrad_0.2 <- readr::read_csv("data/sdm_degrad_0.2.csv") |>
  tibble::column_to_rownames(var = "...1")

#sdm_degrad_0 <- degraded_predictions(sdm_raw, w = 0)
#write.csv(sdm_degrad_0, file="data/sdm_degrad_0.csv")
sdm_degrad_0 <- readr::read_csv("data/sdm_degrad_0.csv") |>
  tibble::column_to_rownames(var = "...1")

sdm_degrad_1_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_degrad_1_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_degrad_1_opth <- do.call(rbind, sdm_degrad_1_opth_list)
sdm_degrad_1_opth <- sdm_degrad_1_opth[rownames(survey_obs), ]

sdm_degrad_0.8_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_degrad_0.8_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_degrad_0.8_opth <- do.call(rbind, sdm_degrad_0.8_opth_list)
sdm_degrad_0.8_opth <- sdm_degrad_0.8_opth[rownames(survey_obs), ]

sdm_degrad_0.6_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_degrad_0.6_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_degrad_0.6_opth <- do.call(rbind, sdm_degrad_0.6_opth_list)
sdm_degrad_0.6_opth <- sdm_degrad_0.6_opth[rownames(survey_obs), ]

sdm_degrad_0.4_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_degrad_0.4_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_degrad_0.4_opth <- do.call(rbind, sdm_degrad_0.4_opth_list)
sdm_degrad_0.4_opth <- sdm_degrad_0.4_opth[rownames(survey_obs), ]

sdm_degrad_0.2_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_degrad_0.2_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_degrad_0.2_opth <- do.call(rbind, sdm_degrad_0.2_opth_list)
sdm_degrad_0.2_opth <- sdm_degrad_0.2_opth[rownames(survey_obs), ]

sdm_degrad_0_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_degrad_0_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_degrad_0_opth <- do.call(rbind, sdm_degrad_0_opth_list)
sdm_degrad_0_opth <- sdm_degrad_0_opth[rownames(survey_obs), ]

## CAMELIA PREDICTIONS
### camelia mem
camelia_mem_raw <- load_camelia_predictions("mem", "raw")
camelia_mem_opth <- load_camelia_predictions("mem", "opth")

camelia_degrad_1_raw <- load_camelia_predictions("degrad_1mem", "raw")
camelia_degrad_1_opth <- load_camelia_predictions("degrad_1mem", "opth")

camelia_degrad_0.8_raw <- load_camelia_predictions("degrad_0.8mem", "raw")
camelia_degrad_0.8_opth <- load_camelia_predictions("degrad_0.8mem", "opth")

camelia_degrad_0.6_raw <- load_camelia_predictions("degrad_0.6mem", "raw")
camelia_degrad_0.6_opth <- load_camelia_predictions("degrad_0.6mem", "opth")

camelia_degrad_0.4_raw <- load_camelia_predictions("degrad_0.4mem", "raw")
camelia_degrad_0.4_opth <- load_camelia_predictions("degrad_0.4mem", "opth")

camelia_degrad_0.2_raw <- load_camelia_predictions("degrad_0.2mem", "raw")
camelia_degrad_0.2_opth <- load_camelia_predictions("degrad_0.2mem", "opth")

camelia_degrad_0_raw <- load_camelia_predictions("degrad_0mem", "raw")
camelia_degrad_0_opth <- load_camelia_predictions("degrad_0mem", "opth")

#write.csv(camelia_mem_raw, "results/camelia_mem_raw.csv")
#write.csv(camelia_mem_opth, "results/camelia_mem_opth.csv")

#write.csv(camelia_degrad_1_raw, "results/camelia_degrad_1_raw.csv")
#write.csv(camelia_degrad_1_opth, "results/camelia_degrad_1_opth.csv")

#write.csv(camelia_degrad_0.8_raw, "results/camelia_degrad_0.8_raw.csv")
#write.csv(camelia_degrad_0.8_opth, "results/camelia_degrad_0.8_opth.csv")

#write.csv(camelia_degrad_0.6_raw, "results/camelia_degrad_0.6_raw.csv")
#write.csv(camelia_degrad_0.6_opth, "results/camelia_degrad_0.6_opth.csv")

#write.csv(camelia_degrad_0.4_raw, "results/camelia_degrad_0.4_raw.csv")
#write.csv(camelia_degrad_0.4_opth, "results/camelia_degrad_0.4_opth.csv")

#write.csv(camelia_degrad_0.2_raw, "results/camelia_degrad_0.2_raw.csv")
#write.csv(camelia_degrad_0.2_opth, "results/camelia_degrad_0.2_opth.csv")

#write.csv(camelia_degrad_0_raw, "results/camelia_degrad_0_raw.csv")
#write.csv(camelia_degrad_0_opth, "results/camelia_degrad_0_opth.csv")

# COMMUNITY INDICES
# Compute indices for all datasets
indices_obs <- compute_indices(survey_obs, "obs")

indices_sdm_raw <- compute_indices(sdm_raw, "sdm_raw")
indices_sdm_opth <- compute_indices(sdm_opth, "sdm_opth")
indices_camelia_mem_raw <- compute_indices(camelia_mem_raw, "camelia_mem_raw")
indices_camelia_mem_opth <- compute_indices(camelia_mem_opth, "camelia_mem_opth")

indices_sdm_degrad_1_raw <- compute_indices(sdm_degrad_1, "sdm_degrad_1")
indices_sdm_degrad_1_opth <- compute_indices(sdm_degrad_1_opth, "sdm_degrad_1_opth")
indices_camelia_degrad_1_raw <- compute_indices(camelia_degrad_1_raw, "camelia_degrad_1_raw")
indices_camelia_degrad_1_opth <- compute_indices(camelia_degrad_1_opth, "camelia_degrad_1_opth")

indices_sdm_degrad_0.8_raw <- compute_indices(sdm_degrad_0.8, "sdm_degrad_0.8")
indices_sdm_degrad_0.8_opth <- compute_indices(sdm_degrad_0.8_opth, "sdm_degrad_0.8_opth")
indices_camelia_degrad_0.8_raw <- compute_indices(camelia_degrad_0.8_raw, "camelia_degrad_0.8_raw")
indices_camelia_degrad_0.8_opth <- compute_indices(camelia_degrad_0.8_opth, "camelia_degrad_0.8_opth")

indices_sdm_degrad_0.6_raw <- compute_indices(sdm_degrad_0.6, "sdm_degrad_0.6")
indices_sdm_degrad_0.6_opth <- compute_indices(sdm_degrad_0.6_opth, "sdm_degrad_0.6_opth")
indices_camelia_degrad_0.6_raw <- compute_indices(camelia_degrad_0.6_raw, "camelia_degrad_0.6_raw")
indices_camelia_degrad_0.6_opth <- compute_indices(camelia_degrad_0.6_opth, "camelia_degrad_0.6_opth")

indices_sdm_degrad_0.4_raw <- compute_indices(sdm_degrad_0.4, "sdm_degrad_0.4")
indices_sdm_degrad_0.4_opth <- compute_indices(sdm_degrad_0.4_opth, "sdm_degrad_0.4_opth")
indices_camelia_degrad_0.4_raw <- compute_indices(camelia_degrad_0.4_raw, "camelia_degrad_0.4_raw")
indices_camelia_degrad_0.4_opth <- compute_indices(camelia_degrad_0.4_opth, "camelia_degrad_0.4_opth")

indices_sdm_degrad_0.2_raw <- compute_indices(sdm_degrad_0.2, "sdm_degrad_0.2")
indices_sdm_degrad_0.2_opth <- compute_indices(sdm_degrad_0.2_opth, "sdm_degrad_0.2_opth")
indices_camelia_degrad_0.2_raw <- compute_indices(camelia_degrad_0.2_raw, "camelia_degrad_0.2_raw")
indices_camelia_degrad_0.2_opth <- compute_indices(camelia_degrad_0.2_opth, "camelia_degrad_0.2_opth")

indices_sdm_degrad_0_raw <- compute_indices(sdm_degrad_0, "sdm_degrad_0")
indices_sdm_degrad_0_opth <- compute_indices(sdm_degrad_0_opth, "sdm_degrad_0_opth")
indices_camelia_degrad_0_raw <- compute_indices(camelia_degrad_0_raw, "camelia_degrad_0_raw")
indices_camelia_degrad_0_opth <- compute_indices(camelia_degrad_0_opth, "camelia_degrad_0_opth")

# Merge CM indices for all methods
idx <- "CM"
CM_tot_raw <- left_join(indices_obs[[idx]], indices_sdm_degrad_1_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.8_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.6_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.4_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.2_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_1_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.8_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.6_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.4_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.2_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0_raw[[idx]], by = "ID")

CM_tot_opth <- left_join(indices_obs[[idx]], indices_sdm_degrad_1_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.8_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.6_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.4_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.2_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_1_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.8_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.6_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.4_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.2_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0_opth[[idx]], by = "ID") %>%
  left_join(CM_mem, by = "ID")

# Merge CSTD indices for all methods
idx <- "CSTD"
CSTD_tot_raw <- left_join(indices_obs[[idx]], indices_sdm_degrad_1_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.8_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.6_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.4_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.2_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_1_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.8_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.6_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.4_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.2_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0_raw[[idx]], by = "ID")

CSTD_tot_opth <- left_join(indices_obs[[idx]], indices_sdm_degrad_1_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.8_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.6_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.4_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.2_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_1_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.8_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.6_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.4_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.2_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0_opth[[idx]], by = "ID") %>%
  left_join(CSTD_mem, by = "ID")

# Merge SR indices for all methods
idx <- "SR"
SR_tot_raw <- left_join(indices_obs[[idx]], indices_sdm_degrad_1_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.8_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.6_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.4_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.2_raw[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_1_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.8_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.6_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.4_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.2_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0_raw[[idx]], by = "ID")

SR_tot_opth <- left_join(indices_obs[[idx]], indices_sdm_degrad_1_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.8_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.6_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.4_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0.2_opth[[idx]], by = "ID") %>%
  left_join(indices_sdm_degrad_0_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_1_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.8_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.6_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.4_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0.2_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_degrad_0_opth[[idx]], by = "ID") %>%
  left_join(SR_mem, by = "ID")

#### -- FIGURE S7 -- ####
# Create a dataframe containing observed and predicted values for different models
figure_raw <- data.frame("obs_CM"=rep(CM_tot_raw$PLH_obs,6),
                         "Pred_CM"=c(CM_tot_raw$PLH_sdm_degrad_0,
                                     CM_tot_raw$PLH_sdm_degrad_0.2,
                                     CM_tot_raw$PLH_sdm_degrad_0.4,
                                     CM_tot_raw$PLH_sdm_degrad_0.6,
                                     CM_tot_raw$PLH_sdm_degrad_0.8,
                                     CM_tot_raw$PLH_sdm_degrad_1),
                         'obs_CSTD'=rep(CSTD_tot_raw$PLH_obs,6),
                         "Pred_CSTD"=c(CSTD_tot_raw$PLH_sdm_degrad_0,
                                       CSTD_tot_raw$PLH_sdm_degrad_0.2,
                                       CSTD_tot_raw$PLH_sdm_degrad_0.4,
                                       CSTD_tot_raw$PLH_sdm_degrad_0.6,
                                       CSTD_tot_raw$PLH_sdm_degrad_0.8,
                                       CSTD_tot_raw$PLH_sdm_degrad_1),
                         'obs_SR' =rep(SR_tot_raw$SR_obs,6),
                         'Pred_SR' =c(SR_tot_raw$SR_sdm_degrad_0,
                                      SR_tot_raw$SR_sdm_degrad_0.2,
                                      SR_tot_raw$SR_sdm_degrad_0.4,
                                      SR_tot_raw$SR_sdm_degrad_0.6,
                                      SR_tot_raw$SR_sdm_degrad_0.8,
                                      SR_tot_raw$SR_sdm_degrad_1),
                         "Models"=c(rep("100% NOISE", 4463), 
                                    rep("80% NOISE", 4463),
                                    rep("60% NOISE", 4463),
                                    rep("40% NOISE", 4463),
                                    rep("20% NOISE", 4463),
                                    rep("0% NOISE", 4463))
)
figure_raw$Models <- factor(figure_raw$Models, levels = c("0% NOISE", "20% NOISE","40% NOISE", "60% NOISE", "80% NOISE", "100% NOISE"))

# Compute R² and RMSE for each model
cor_metrics <- data.frame(
  Models = factor(c("0% NOISE", "20% NOISE","40% NOISE", "60% NOISE", "80% NOISE", "100% NOISE"), levels = c("0% NOISE", "20% NOISE","40% NOISE", "60% NOISE", "80% NOISE", "100% NOISE")),
  R2_CM = sapply(list(CM_tot_raw$PLH_sdm_degrad_0, 
                      CM_tot_raw$PLH_sdm_degrad_0.2,
                      CM_tot_raw$PLH_sdm_degrad_0.4,
                      CM_tot_raw$PLH_sdm_degrad_0.6,
                      CM_tot_raw$PLH_sdm_degrad_0.8,
                      CM_tot_raw$PLH_sdm_degrad_1), function(x) cor(CM_tot_raw$PLH_obs, x)^2),
  RMSE_CM = sapply(list(CM_tot_raw$PLH_sdm_degrad_0, 
                        CM_tot_raw$PLH_sdm_degrad_0.2,
                        CM_tot_raw$PLH_sdm_degrad_0.4,
                        CM_tot_raw$PLH_sdm_degrad_0.6,
                        CM_tot_raw$PLH_sdm_degrad_0.8,
                        CM_tot_raw$PLH_sdm_degrad_1), function(x) rmse(CM_tot_raw$PLH_obs, x)),
  R2_CSTD = sapply(list(CSTD_tot_raw$PLH_sdm_degrad_0, 
                        CSTD_tot_raw$PLH_sdm_degrad_0.2,
                        CSTD_tot_raw$PLH_sdm_degrad_0.4,
                        CSTD_tot_raw$PLH_sdm_degrad_0.6,
                        CSTD_tot_raw$PLH_sdm_degrad_0.8,
                        CSTD_tot_raw$PLH_sdm_degrad_1), function(x) cor(CSTD_tot_raw$PLH_obs, x)^2),
  RMSE_CSTD = sapply(list(CSTD_tot_raw$PLH_sdm_degrad_0,
                          CSTD_tot_raw$PLH_sdm_degrad_0.2,
                          CSTD_tot_raw$PLH_sdm_degrad_0.4,
                          CSTD_tot_raw$PLH_sdm_degrad_0.6,
                          CSTD_tot_raw$PLH_sdm_degrad_0.8,
                          CSTD_tot_raw$PLH_sdm_degrad_1), function(x) rmse(CSTD_tot_raw$PLH_obs, x)),
  R2_SR = sapply(list(SR_tot_raw$SR_sdm_degrad_0, 
                      SR_tot_raw$SR_sdm_degrad_0.2,
                      SR_tot_raw$SR_sdm_degrad_0.4,
                      SR_tot_raw$SR_sdm_degrad_0.6,
                      SR_tot_raw$SR_sdm_degrad_0.8,
                      SR_tot_raw$SR_sdm_degrad_1), function(x) cor(SR_tot_raw$SR_obs, x)^2),
  RMSE_SR = sapply(list(SR_tot_raw$SR_sdm_degrad_0, 
                        SR_tot_raw$SR_sdm_degrad_0.2,
                        SR_tot_raw$SR_sdm_degrad_0.4,
                        SR_tot_raw$SR_sdm_degrad_0.6,
                        SR_tot_raw$SR_sdm_degrad_0.8,
                        SR_tot_raw$SR_sdm_degrad_1), function(x) rmse(SR_tot_raw$SR_obs, x))
)

plot_fig_S7 <- function(obs, pred, lab_text, cor_column, rmse_column) {
  ggplot(figure_raw, aes(y = obs, x = pred, color = Models, fill = Models)) +
    geom_hdr() +
    geom_point(shape = 21, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("darkblue", "#2A4FAF", "#558CD3", "#81BFE7", "#ACDEEF", "#D8F6FF")) +
    scale_fill_manual(values = c("darkblue", "#2A4FAF", "#558CD3", "#81BFE7", "#ACDEEF", "#D8F6FF")) +
    facet_wrap(vars(Models), nrow = 1) +
    ylab(paste0("Observed ",lab_text)) +
    xlab(paste0("Predicted ",lab_text)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      strip.text = element_text(size = 8)
    ) +
    geom_label(data = cor_metrics,
               aes(x = -Inf, y = Inf,
                   label = sprintf("R² = %.2f\nRMSE = %.2f", round(get(cor_column), 2), round(get(rmse_column), 2))),
               color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)
}

# Generate plots
CM_raw <- plot_fig_S7(figure_raw$obs_CM, figure_raw$Pred_CM, "CM plant height", "R2_CM", "RMSE_CM")
CSTD_raw <- plot_fig_S7(figure_raw$obs_CSTD, figure_raw$Pred_CSTD, "CSTD plant height", "R2_CSTD", "RMSE_CSTD")
SR_raw <- plot_fig_S7(figure_raw$obs_SR, figure_raw$Pred_SR, "Species Richness", "R2_SR", "RMSE_SR")

# Combine plots
combined_plot <- SR_raw / CM_raw / CSTD_raw +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("figure_S7_sdm_degraded_predictions.pdf", 
       plot = combined_plot, 
       device = "pdf", 
       width = 8.27,
       height = 11.69,
       dpi = 300)

ggsave("figure_S7_sdm_degraded_predictions.png", 
       plot = combined_plot, 
       width = 8.27,
       height = 11.69,
       dpi = 300)

#### -- FIGURE S8 -- ####
# Create a dataframe containing observed and predicted with camelia values for different models
figure_raw <- data.frame("obs_CM"=rep(CM_tot_raw$PLH_obs,6),
                         "Pred_CM"=c(CM_tot_raw$PLH_camelia_degrad_0_raw,
                                     CM_tot_raw$PLH_camelia_degrad_0.2_raw,
                                     CM_tot_raw$PLH_camelia_degrad_0.4_raw,
                                     CM_tot_raw$PLH_camelia_degrad_0.6_raw,
                                     CM_tot_raw$PLH_camelia_degrad_0.8_raw,
                                     CM_tot_raw$PLH_camelia_degrad_1_raw),
                         'obs_CSTD'=rep(CSTD_tot_raw$PLH_obs,6),
                         "Pred_CSTD"=c(CSTD_tot_raw$PLH_camelia_degrad_0_raw,
                                       CSTD_tot_raw$PLH_camelia_degrad_0.2_raw,
                                       CSTD_tot_raw$PLH_camelia_degrad_0.4_raw,
                                       CSTD_tot_raw$PLH_camelia_degrad_0.6_raw,
                                       CSTD_tot_raw$PLH_camelia_degrad_0.8_raw,
                                       CSTD_tot_raw$PLH_camelia_degrad_1_raw),
                         'obs_SR' =rep(SR_tot_raw$SR_obs,6),
                         'Pred_SR' =c(SR_tot_raw$SR_camelia_degrad_0_raw,
                                      SR_tot_raw$SR_camelia_degrad_0.2_raw,
                                      SR_tot_raw$SR_camelia_degrad_0.4_raw,
                                      SR_tot_raw$SR_camelia_degrad_0.6_raw,
                                      SR_tot_raw$SR_camelia_degrad_0.8_raw,
                                      SR_tot_raw$SR_camelia_degrad_1_raw),
                         "Models"=c(rep("0% NOISE", 4463), 
                                    rep("20% NOISE", 4463),
                                    rep("40% NOISE", 4463),
                                    rep("60% NOISE", 4463),
                                    rep("80% NOISE", 4463),
                                    rep("100% NOISE", 4463))
)
figure_raw$Models <- factor(figure_raw$Models, levels = c("0% NOISE", "20% NOISE","40% NOISE", "60% NOISE", "80% NOISE", "100% NOISE"))

# Compute R² and RMSE for each model
cor_metrics <- data.frame(
  Models = factor(c("0% NOISE", "20% NOISE","40% NOISE", "60% NOISE", "80% NOISE", "100% NOISE"), levels = c("0% NOISE", "20% NOISE","40% NOISE", "60% NOISE", "80% NOISE", "100% NOISE")),
  R2_CM = sapply(list(CM_tot_raw$PLH_camelia_degrad_0_raw, 
                      CM_tot_raw$PLH_camelia_degrad_0.2_raw,
                      CM_tot_raw$PLH_camelia_degrad_0.4_raw,
                      CM_tot_raw$PLH_camelia_degrad_0.6_raw,
                      CM_tot_raw$PLH_camelia_degrad_0.8_raw,
                      CM_tot_raw$PLH_camelia_degrad_1_raw), function(x) cor(CM_tot_raw$PLH_obs, x)^2),
  RMSE_CM = sapply(list(CM_tot_raw$PLH_camelia_degrad_0_raw, 
                        CM_tot_raw$PLH_camelia_degrad_0.2_raw,
                        CM_tot_raw$PLH_camelia_degrad_0.4_raw,
                        CM_tot_raw$PLH_camelia_degrad_0.6_raw,
                        CM_tot_raw$PLH_camelia_degrad_0.8_raw,
                        CM_tot_raw$PLH_camelia_degrad_1_raw), function(x) rmse(CM_tot_raw$PLH_obs, x)),
  R2_CSTD = sapply(list(CSTD_tot_raw$PLH_camelia_degrad_0_raw, 
                        CSTD_tot_raw$PLH_camelia_degrad_0.2_raw,
                        CSTD_tot_raw$PLH_camelia_degrad_0.4_raw,
                        CSTD_tot_raw$PLH_camelia_degrad_0.6_raw,
                        CSTD_tot_raw$PLH_camelia_degrad_0.8_raw,
                        CSTD_tot_raw$PLH_camelia_degrad_1_raw), function(x) cor(CSTD_tot_raw$PLH_obs, x)^2),
  RMSE_CSTD = sapply(list(CSTD_tot_raw$PLH_camelia_degrad_0_raw,
                          CSTD_tot_raw$PLH_camelia_degrad_0.2_raw,
                          CSTD_tot_raw$PLH_camelia_degrad_0.4_raw,
                          CSTD_tot_raw$PLH_camelia_degrad_0.6_raw,
                          CSTD_tot_raw$PLH_camelia_degrad_0.8_raw,
                          CSTD_tot_raw$PLH_camelia_degrad_1_raw), function(x) rmse(CSTD_tot_raw$PLH_obs, x)),
  R2_SR = sapply(list(SR_tot_raw$SR_camelia_degrad_0_raw, 
                      SR_tot_raw$SR_camelia_degrad_0.2_raw,
                      SR_tot_raw$SR_camelia_degrad_0.4_raw,
                      SR_tot_raw$SR_camelia_degrad_0.6_raw,
                      SR_tot_raw$SR_camelia_degrad_0.8_raw,
                      SR_tot_raw$SR_camelia_degrad_1_raw), function(x) cor(SR_tot_raw$SR_obs, x)^2),
  RMSE_SR = sapply(list(SR_tot_raw$SR_camelia_degrad_0_raw, 
                        SR_tot_raw$SR_camelia_degrad_0.2_raw,
                        SR_tot_raw$SR_camelia_degrad_0.4_raw,
                        SR_tot_raw$SR_camelia_degrad_0.6_raw,
                        SR_tot_raw$SR_camelia_degrad_0.8_raw,
                        SR_tot_raw$SR_camelia_degrad_1_raw), function(x) rmse(SR_tot_raw$SR_obs, x))
)

plot_fig_S8 <- function(obs, pred, lab_text, cor_column, rmse_column) {
  ggplot(figure_raw, aes(y = obs, x = pred, color = Models, fill = Models)) +
    geom_hdr() +
    geom_point(shape = 21, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("darkblue", "#2A4FAF", "#558CD3", "#81BFE7", "#ACDEEF", "#D8F6FF")) +
    scale_fill_manual(values = c("darkblue", "#2A4FAF", "#558CD3", "#81BFE7", "#ACDEEF", "#D8F6FF")) +
    facet_wrap(vars(Models), nrow = 1) +
    ylab(paste0("Observed ",lab_text)) +
    xlab(paste0("Predicted ",lab_text)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      strip.text = element_text(size = 8)
    ) +
    geom_label(data = cor_metrics,
               aes(x = -Inf, y = Inf,
                   label = sprintf("R² = %.2f\nRMSE = %.2f", round(get(cor_column), 2), round(get(rmse_column), 2))),
               color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)
}

# Generate plots
CM_raw <- plot_fig_S8(figure_raw$obs_CM, figure_raw$Pred_CM, "CM plant height", "R2_CM", "RMSE_CM")
CSTD_raw <- plot_fig_S8(figure_raw$obs_CSTD, figure_raw$Pred_CSTD, "CSTD plant height", "R2_CSTD", "RMSE_CSTD")
SR_raw <- plot_fig_S8(figure_raw$obs_SR, figure_raw$Pred_SR, "Species Richness", "R2_SR", "RMSE_SR")

# Combine plots
combined_plot <- SR_raw / CM_raw / CSTD_raw +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("figure_S8_camelia_degraded_predictions.pdf", 
       plot = combined_plot, 
       device = "pdf", 
       width = 8.27,
       height = 11.69,
       dpi = 300)

ggsave("figure_S8_camelia_degraded_predictions.png", 
       plot = combined_plot, 
       width = 8.27,
       height = 11.69,
       dpi = 300)

# Species predictions evaluation : 
## RAW
AUC_raw <- data.frame("Species"=rep(colnames(survey_obs),12), "AUC"=rep(0,831*12), "pr-AUC"=rep(0,831*12), "Approach"=c(rep("DEGRAD w=1",831), rep("DEGRAD w=0.8",831), rep("DEGRAD w=0.6",831), rep("DEGRAD w=0.4",831), rep("DEGRAD w=0.2",831), rep("DEGRAD w=0",831), rep("CAMELIA DEGRAD w=1",831), rep("CAMELIA DEGRAD w=0.8",831), rep("CAMELIA DEGRAD w=0.6",831), rep("CAMELIA DEGRAD w=0.4",831), rep("CAMELIA DEGRAD w=0.2",831), rep("CAMELIA DEGRAD w=0",831)))
for (i in 1:831){
  AUC_raw$AUC[i] = auc(survey_obs[,i], as.numeric(unlist(sdm_degrad_1[,i])))
  pr_curve <- pr.curve(as.numeric(unlist(sdm_degrad_1[,i]))[survey_obs[,i]==1],as.numeric(unlist(sdm_degrad_1[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[i]= pr_curve$auc.integral
  AUC_raw$AUC[831+i] = auc(survey_obs[,i], as.numeric(unlist(sdm_degrad_0.8[,i])))
  pr_curve <- pr.curve(as.numeric(unlist(sdm_degrad_0.8[,i]))[survey_obs[,i]==1],as.numeric(unlist(sdm_degrad_0.8[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*2+i] = auc(survey_obs[,i], as.numeric(unlist(sdm_degrad_0.6[,i])))
  pr_curve <- pr.curve(as.numeric(unlist(sdm_degrad_0.6[,i]))[survey_obs[,i]==1],as.numeric(unlist(sdm_degrad_0.6[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*2+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*3+i] = auc(survey_obs[,i], as.numeric(unlist(sdm_degrad_0.4[,i])))
  pr_curve <- pr.curve(as.numeric(unlist(sdm_degrad_0.4[,i]))[survey_obs[,i]==1],as.numeric(unlist(sdm_degrad_0.4[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*3+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*4+i] = auc(survey_obs[,i], as.numeric(unlist(sdm_degrad_0.2[,i])))
  pr_curve <- pr.curve(as.numeric(unlist(sdm_degrad_0.2[,i]))[survey_obs[,i]==1],as.numeric(unlist(sdm_degrad_0.2[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*4+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*5+i] = auc(survey_obs[,i], as.numeric(unlist(sdm_degrad_0[,i])))
  pr_curve <- pr.curve(as.numeric(unlist(sdm_degrad_0[,i]))[survey_obs[,i]==1],as.numeric(unlist(sdm_degrad_0[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*5+i]= pr_curve$auc.integral  
  AUC_raw$AUC[831*6+i]=auc(survey_obs[,i], camelia_degrad_1_raw[,i])
  pr_curve <- pr.curve(as.numeric(unlist(camelia_degrad_1_raw[,i]))[survey_obs[,i]==1],as.numeric(unlist(camelia_degrad_1_raw[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*6+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*7+i]=auc(survey_obs[,i], camelia_degrad_0.8_raw[,i])
  pr_curve <- pr.curve(as.numeric(unlist(camelia_degrad_0.8_raw[,i]))[survey_obs[,i]==1],as.numeric(unlist(camelia_degrad_0.8_raw[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*7+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*8+i]=auc(survey_obs[,i], camelia_degrad_0.6_raw[,i])
  pr_curve <- pr.curve(as.numeric(unlist(camelia_degrad_0.6_raw[,i]))[survey_obs[,i]==1],as.numeric(unlist(camelia_degrad_0.6_raw[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*8+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*9+i]=auc(survey_obs[,i], camelia_degrad_0.4_raw[,i])
  pr_curve <- pr.curve(as.numeric(unlist(camelia_degrad_0.4_raw[,i]))[survey_obs[,i]==1],as.numeric(unlist(camelia_degrad_0.4_raw[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*9+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*10+i]=auc(survey_obs[,i], camelia_degrad_0.2_raw[,i])
  pr_curve <- pr.curve(as.numeric(unlist(camelia_degrad_0.2_raw[,i]))[survey_obs[,i]==1],as.numeric(unlist(camelia_degrad_0.2_raw[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*10+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*11+i]=auc(survey_obs[,i], camelia_degrad_0_raw[,i])
  pr_curve <- pr.curve(as.numeric(unlist(camelia_degrad_0_raw[,i]))[survey_obs[,i]==1],as.numeric(unlist(camelia_degrad_0_raw[,i]))[survey_obs[,i]==0], curve=T)
  AUC_raw$pr.AUC[831*11+i]= pr_curve$auc.integral
  print(i)
}

## BINARIZED USING PER SPECIES TSS-OPTIMIZED THRESHOLD
eval_species <- data.frame("species"=colnames(survey_obs), 
                           "prevalence"=colSums(survey_obs),
                           "Specificite_SDM_degrad_1"=rep(0,831),
                           "Sensibilite_SDM_degrad_1"=rep(0,831),
                           "Precision_SDM_degrad_1"=rep(0,831),
                           "TSS_SDM_degrad_1"=rep(0,831),
                           "Specificite_SDM_degrad_0.8"=rep(0,831),
                           "Sensibilite_SDM_degrad_0.8"=rep(0,831),
                           "Precision_SDM_degrad_0.8"=rep(0,831),
                           "TSS_SDM_degrad_0.8"=rep(0,831),
                           "Specificite_SDM_degrad_0.6"=rep(0,831),
                           "Sensibilite_SDM_degrad_0.6"=rep(0,831),
                           "Precision_SDM_degrad_0.6"=rep(0,831),
                           "TSS_SDM_degrad_0.6"=rep(0,831),
                           "Specificite_SDM_degrad_0.4"=rep(0,831),
                           "Sensibilite_SDM_degrad_0.4"=rep(0,831),
                           "Precision_SDM_degrad_0.4"=rep(0,831),
                           "TSS_SDM_degrad_0.4"=rep(0,831),
                           "Specificite_SDM_degrad_0.2"=rep(0,831),
                           "Sensibilite_SDM_degrad_0.2"=rep(0,831),
                           "Precision_SDM_degrad_0.2"=rep(0,831),
                           "TSS_SDM_degrad_0.2"=rep(0,831),
                           "Specificite_SDM_degrad_0"=rep(0,831),
                           "Sensibilite_SDM_degrad_0"=rep(0,831),
                           "Precision_SDM_degrad_0"=rep(0,831),
                           "TSS_SDM_degrad_0"=rep(0,831),
                           "Specificite_CAMELIA_SDM_degrad_1"=rep(0,831),
                           "Sensibilite_CAMELIA_SDM_degrad_1"=rep(0,831),
                           "Precision_CAMELIA_SDM_degrad_1"=rep(0,831),
                           "TSS_CAMELIA_SDM_degrad_1"=rep(0,831),
                           "Specificite_CAMELIA_SDM_degrad_0.8"=rep(0,831),
                           "Sensibilite_CAMELIA_SDM_degrad_0.8"=rep(0,831),
                           "Precision_CAMELIA_SDM_degrad_0.8"=rep(0,831),
                           "TSS_CAMELIA_SDM_degrad_0.8"=rep(0,831),
                           "Specificite_CAMELIA_SDM_degrad_0.6"=rep(0,831),
                           "Sensibilite_CAMELIA_SDM_degrad_0.6"=rep(0,831),
                           "Precision_CAMELIA_SDM_degrad_0.6"=rep(0,831),
                           "TSS_CAMELIA_SDM_degrad_0.6"=rep(0,831),
                           "Specificite_CAMELIA_SDM_degrad_0.4"=rep(0,831),
                           "Sensibilite_CAMELIA_SDM_degrad_0.4"=rep(0,831),
                           "Precision_CAMELIA_SDM_degrad_0.4"=rep(0,831),
                           "TSS_CAMELIA_SDM_degrad_0.4"=rep(0,831),
                           "Specificite_CAMELIA_SDM_degrad_0.2"=rep(0,831),
                           "Sensibilite_CAMELIA_SDM_degrad_0.2"=rep(0,831),
                           "Precision_CAMELIA_SDM_degrad_0.2"=rep(0,831),
                           "TSS_CAMELIA_SDM_degrad_0.2"=rep(0,831),
                           "Specificite_CAMELIA_SDM_degrad_0"=rep(0,831),
                           "Sensibilite_CAMELIA_SDM_degrad_0"=rep(0,831),
                           "Precision_CAMELIA_SDM_degrad_0"=rep(0,831),
                           "TSS_CAMELIA_SDM_degrad_0"=rep(0,831))

eval_species_opth <- eval_species
for (i in 1:831){
  VP=length(which(which(sdm_degrad_1_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(sdm_degrad_1_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(sdm_degrad_1_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(sdm_degrad_1_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SDM_degrad_1[i]=VP/(VP+FN)
  eval_species_opth$Specificite_SDM_degrad_1[i]=VN/(VN+FP)
  eval_species_opth$Precision_SDM_degrad_1[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_SDM_degrad_1[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(sdm_degrad_0.8_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(sdm_degrad_0.8_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(sdm_degrad_0.8_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(sdm_degrad_0.8_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SDM_degrad_0.8[i]=VP/(VP+FN)
  eval_species_opth$Specificite_SDM_degrad_0.8[i]=VN/(VN+FP)
  eval_species_opth$Precision_SDM_degrad_0.8[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_SDM_degrad_0.8[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(sdm_degrad_0.6_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(sdm_degrad_0.6_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(sdm_degrad_0.6_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(sdm_degrad_0.6_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SDM_degrad_0.6[i]=VP/(VP+FN)
  eval_species_opth$Specificite_SDM_degrad_0.6[i]=VN/(VN+FP)
  eval_species_opth$Precision_SDM_degrad_0.6[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_SDM_degrad_0.6[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(sdm_degrad_0.4_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(sdm_degrad_0.4_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(sdm_degrad_0.4_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(sdm_degrad_0.4_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SDM_degrad_0.4[i]=VP/(VP+FN)
  eval_species_opth$Specificite_SDM_degrad_0.4[i]=VN/(VN+FP)
  eval_species_opth$Precision_SDM_degrad_0.4[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_SDM_degrad_0.4[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(sdm_degrad_0.2_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(sdm_degrad_0.2_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(sdm_degrad_0.2_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(sdm_degrad_0.2_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SDM_degrad_0.2[i]=VP/(VP+FN)
  eval_species_opth$Specificite_SDM_degrad_0.2[i]=VN/(VN+FP)
  eval_species_opth$Precision_SDM_degrad_0.2[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_SDM_degrad_0.2[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(sdm_degrad_0_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(sdm_degrad_0_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(sdm_degrad_0_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(sdm_degrad_0_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SDM_degrad_0[i]=VP/(VP+FN)
  eval_species_opth$Specificite_SDM_degrad_0[i]=VN/(VN+FP)
  eval_species_opth$Precision_SDM_degrad_0[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_SDM_degrad_0[i] = 2 * VP/(2*VP+FP+FN)
  
  
  VP=length(which(which(camelia_degrad_1_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(camelia_degrad_1_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(camelia_degrad_1_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(camelia_degrad_1_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_SDM_1[i]=VP/(VP+FN)
  eval_species_opth$Specificite_CAMELIA_SDM_1[i]=VN/(VN+FP)
  eval_species_opth$Precision_CAMELIA_SDM_1[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_CAMELIA_SDM_1[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(camelia_degrad_0.8_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(camelia_degrad_0.8_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(camelia_degrad_0.8_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(camelia_degrad_0.8_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_SDM_0.8[i]=VP/(VP+FN)
  eval_species_opth$Specificite_CAMELIA_SDM_0.8[i]=VN/(VN+FP)
  eval_species_opth$Precision_CAMELIA_SDM_0.8[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_CAMELIA_SDM_0.8[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(camelia_degrad_0.6_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(camelia_degrad_0.6_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(camelia_degrad_0.6_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(camelia_degrad_0.6_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_SDM_0.6[i]=VP/(VP+FN)
  eval_species_opth$Specificite_CAMELIA_SDM_0.6[i]=VN/(VN+FP)
  eval_species_opth$Precision_CAMELIA_SDM_0.6[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_CAMELIA_SDM_0.6[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(camelia_degrad_0.4_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(camelia_degrad_0.4_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(camelia_degrad_0.4_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(camelia_degrad_0.4_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_SDM_0.4[i]=VP/(VP+FN)
  eval_species_opth$Specificite_CAMELIA_SDM_0.4[i]=VN/(VN+FP)
  eval_species_opth$Precision_CAMELIA_SDM_0.4[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_CAMELIA_SDM_0.4[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(camelia_degrad_0.2_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(camelia_degrad_0.2_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(camelia_degrad_0.2_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(camelia_degrad_0.2_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_SDM_0.2[i]=VP/(VP+FN)
  eval_species_opth$Specificite_CAMELIA_SDM_0.2[i]=VN/(VN+FP)
  eval_species_opth$Precision_CAMELIA_SDM_0.2[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_CAMELIA_SDM_0.2[i] = 2 * VP/(2*VP+FP+FN)
  
  VP=length(which(which(camelia_degrad_0_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP=length(which(which(camelia_degrad_0_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN=length(which(which(camelia_degrad_0_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN=length(which(which(camelia_degrad_0_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_SDM_0[i]=VP/(VP+FN)
  eval_species_opth$Specificite_CAMELIA_SDM_0[i]=VN/(VN+FP)
  eval_species_opth$Precision_CAMELIA_SDM_0[i]=VP/(VP+FP)
  eval_species_opth$Sorensen_CAMELIA_SDM_0[i] = 2 * VP/(2*VP+FP+FN)
  print(i)
}
eval_species_opth$TSS_SDM_degrad_1 = eval_species_opth$Specificite_SDM_degrad_1+eval_species_opth$Sensibilite_SDM_degrad_1-1
eval_species_opth$TSS_SDM_degrad_0.8 = eval_species_opth$Specificite_SDM_degrad_0.8+eval_species_opth$Sensibilite_SDM_degrad_0.8-1
eval_species_opth$TSS_SDM_degrad_0.6 = eval_species_opth$Specificite_SDM_degrad_0.6+eval_species_opth$Sensibilite_SDM_degrad_0.6-1
eval_species_opth$TSS_SDM_degrad_0.4 = eval_species_opth$Specificite_SDM_degrad_0.4+eval_species_opth$Sensibilite_SDM_degrad_0.4-1
eval_species_opth$TSS_SDM_degrad_0.2 = eval_species_opth$Specificite_SDM_degrad_0.2+eval_species_opth$Sensibilite_SDM_degrad_0.2-1
eval_species_opth$TSS_SDM_degrad_0 = eval_species_opth$Specificite_SDM_degrad_0+eval_species_opth$Sensibilite_SDM_degrad_0-1
eval_species_opth$TSS_CAMELIA_SDM_1 = eval_species_opth$Specificite_CAMELIA_SDM_1+eval_species_opth$Sensibilite_CAMELIA_SDM_1-1
eval_species_opth$TSS_CAMELIA_SDM_0.8 = eval_species_opth$Specificite_CAMELIA_SDM_0.8+eval_species_opth$Sensibilite_CAMELIA_SDM_0.8-1
eval_species_opth$TSS_CAMELIA_SDM_0.6 = eval_species_opth$Specificite_CAMELIA_SDM_0.6+eval_species_opth$Sensibilite_CAMELIA_SDM_0.6-1
eval_species_opth$TSS_CAMELIA_SDM_0.4 = eval_species_opth$Specificite_CAMELIA_SDM_0.4+eval_species_opth$Sensibilite_CAMELIA_SDM_0.4-1
eval_species_opth$TSS_CAMELIA_SDM_0.2 = eval_species_opth$Specificite_CAMELIA_SDM_0.2+eval_species_opth$Sensibilite_CAMELIA_SDM_0.2-1
eval_species_opth$TSS_CAMELIA_SDM_0 = eval_species_opth$Specificite_CAMELIA_SDM_0+eval_species_opth$Sensibilite_CAMELIA_SDM_0-1

table_S4 <- data.frame(
  "Stacking processus" = c(rep("Probability", 2),rep('Binazised by TSS threshold', 4)),
  "Indices" = c("mean AUC", "mean pr-AUC", "Sensibilite", "Specificite", "Precision", "TSS"),
  
  "DEGRAD w=1" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "DEGRAD w=1"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "DEGRAD w=1"]),
    mean(eval_species_opth$Sensibilite_SDM_degrad_1),
    mean(eval_species_opth$Specificite_SDM_degrad_1),
    mean(eval_species_opth$Precision_SDM_degrad_1),
    mean(eval_species_opth$TSS_SDM_degrad_1)
  ),
  "DEGRAD w=0.8" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "DEGRAD w=0.8"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "DEGRAD w=0.8"]),
    mean(eval_species_opth$Sensibilite_SDM_degrad_0.8),
    mean(eval_species_opth$Specificite_SDM_degrad_0.8),
    mean(eval_species_opth$Precision_SDM_degrad_0.8),
    mean(eval_species_opth$TSS_SDM_degrad_0.8)
  ),
  "DEGRAD w=0.6" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "DEGRAD w=0.6"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "DEGRAD w=0.6"]),
    mean(eval_species_opth$Sensibilite_SDM_degrad_0.6),
    mean(eval_species_opth$Specificite_SDM_degrad_0.6),
    mean(eval_species_opth$Precision_SDM_degrad_0.6),
    mean(eval_species_opth$TSS_SDM_degrad_0.6)
  ),
  "DEGRAD w=0.4" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "DEGRAD w=0.4"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "DEGRAD w=0.4"]),
    mean(eval_species_opth$Sensibilite_SDM_degrad_0.4),
    mean(eval_species_opth$Specificite_SDM_degrad_0.4),
    mean(eval_species_opth$Precision_SDM_degrad_0.4),
    mean(eval_species_opth$TSS_SDM_degrad_0.4)
  ),
  "DEGRAD w=0.2" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "DEGRAD w=0.2"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "DEGRAD w=0.2"]),
    mean(eval_species_opth$Sensibilite_SDM_degrad_0.2),
    mean(eval_species_opth$Specificite_SDM_degrad_0.2),
    mean(eval_species_opth$Precision_SDM_degrad_0.2),
    mean(eval_species_opth$TSS_SDM_degrad_0.2)
  ),
  "DEGRAD w=0" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "DEGRAD w=0"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "DEGRAD w=0"]),
    mean(eval_species_opth$Sensibilite_SDM_degrad_0),
    mean(eval_species_opth$Specificite_SDM_degrad_0),
    mean(eval_species_opth$Precision_SDM_degrad_0),
    mean(eval_species_opth$TSS_SDM_degrad_0)
  ),
  
  "CAMELIA DEGRAD w=1" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=1"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=1"]),
    mean(eval_species_opth$Sensibilite_CAMELIA_SDM_1),
    mean(eval_species_opth$Specificite_CAMELIA_SDM_1),
    mean(eval_species_opth$Precision_CAMELIA_SDM_1),
    mean(eval_species_opth$TSS_CAMELIA_SDM_1)
  ),
  "CAMELIA DEGRAD w=0.8" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.8"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.8"]),
    mean(eval_species_opth$Sensibilite_CAMELIA_SDM_0.8),
    mean(eval_species_opth$Specificite_CAMELIA_SDM_0.8),
    mean(eval_species_opth$Precision_CAMELIA_SDM_0.8),
    mean(eval_species_opth$TSS_CAMELIA_SDM_0.8)
  ),
  "CAMELIA DEGRAD w=0.6" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.6"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.6"]),
    mean(eval_species_opth$Sensibilite_CAMELIA_SDM_0.6),
    mean(eval_species_opth$Specificite_CAMELIA_SDM_0.6),
    mean(eval_species_opth$Precision_CAMELIA_SDM_0.6),
    mean(eval_species_opth$TSS_CAMELIA_SDM_0.6)
  ),
  "CAMELIA DEGRAD w=0.4" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.4"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.4"]),
    mean(eval_species_opth$Sensibilite_CAMELIA_SDM_0.4),
    mean(eval_species_opth$Specificite_CAMELIA_SDM_0.4),
    mean(eval_species_opth$Precision_CAMELIA_SDM_0.4),
    mean(eval_species_opth$TSS_CAMELIA_SDM_0.4)
  ),
  "CAMELIA DEGRAD w=0.2" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.2"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0.2"]),
    mean(eval_species_opth$Sensibilite_CAMELIA_SDM_0.2),
    mean(eval_species_opth$Specificite_CAMELIA_SDM_0.2),
    mean(eval_species_opth$Precision_CAMELIA_SDM_0.2),
    mean(eval_species_opth$TSS_CAMELIA_SDM_0.2)
  ),
  "CAMELIA DEGRAD w=0" = c(
    mean(AUC_raw$AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0"]), 
    mean(AUC_raw$pr.AUC[AUC_raw$Approach == "CAMELIA DEGRAD w=0"]),
    mean(eval_species_opth$Sensibilite_CAMELIA_SDM_0),
    mean(eval_species_opth$Specificite_CAMELIA_SDM_0),
    mean(eval_species_opth$Precision_CAMELIA_SDM_0),
    mean(eval_species_opth$TSS_CAMELIA_SDM_0)
  )
)


write.csv(table_S4, file="table_S9_species_predictions.csv")
