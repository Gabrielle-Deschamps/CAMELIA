## Dependancy
library(tidyverse)
library(ggplot2)
library(ggdensity)
library(scales)
library(patchwork)
library(ggpubr)
library(readr)
source("utilities/FONCTIONS.R")

# DATA
## OBS
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)
survey_obs <- read.csv("data/survey_obs.csv", header = TRUE, row.names = 1)
traits <- read.csv("data/traits.csv", header = TRUE, row.names = 1)

## SSDM
sdm_raw_1 <- read.csv("results/sdm_raw_1.csv", header = TRUE, row.names = 1)
sdm_raw_2 <- read.csv("results/sdm_raw_2.csv", header = TRUE, row.names = 1)
sdm_raw_3 <- read.csv("results/sdm_raw_3.csv", header = TRUE, row.names = 1)
sdm_raw_4 <- read.csv("results/sdm_raw_4.csv", header = TRUE, row.names = 1)
sdm_raw <- rbind(sdm_raw_1, sdm_raw_2, sdm_raw_3, sdm_raw_4)
sdm_raw <- sdm_raw[rownames(survey_obs), ]

sdm_opth_1 <- read.csv("results/sdm_opth_1.csv", header = TRUE, row.names = 1)
sdm_opth_2 <- read.csv("results/sdm_opth_2.csv", header = TRUE, row.names = 1)
sdm_opth_3 <- read.csv("results/sdm_opth_3.csv", header = TRUE, row.names = 1)
sdm_opth_4 <- read.csv("results/sdm_opth_4.csv", header = TRUE, row.names = 1)
sdm_opth <- rbind(sdm_opth_1, sdm_opth_2, sdm_opth_3, sdm_opth_4)
sdm_opth <- sdm_opth[rownames(survey_obs), ]

sdm_prr <- prr(sdm_raw, SR_mem)

## CAMELIA PREDICTIONS
### camelia obs
camelia_obs_raw1 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_1/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_obs_raw2 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_2/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_obs_raw3 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_3/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_obs_raw4 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_4/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_obs_raw <- rbind(camelia_obs_raw1, camelia_obs_raw2, camelia_obs_raw3, camelia_obs_raw4)
camelia_obs_raw <- camelia_obs_raw[rownames(survey_obs), ]

camelia_obs_opth1 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_1/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_obs_opth2 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_2/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_obs_opth3 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_3/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_obs_opth4 <- read.csv("results/CAMELIA_OBS/mlp1_bce_msle/fold_4/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_obs_opth <- rbind(camelia_obs_opth1, camelia_obs_opth2, camelia_obs_opth3, camelia_obs_opth4)
camelia_obs_opth <- camelia_obs_opth[rownames(survey_obs), ]

camelia_obs_prr <- prr(camelia_obs_raw, SR_mem)

### camelia mem
camelia_mem_raw1 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_1/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_mem_raw2 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_2/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_mem_raw3 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_3/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_mem_raw4 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_4/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_mem_raw <- rbind(camelia_mem_raw1, camelia_mem_raw2, camelia_mem_raw3, camelia_mem_raw4)
camelia_mem_raw <- camelia_mem_raw[rownames(survey_obs), ]

camelia_mem_opth1 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_1/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_mem_opth2 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_2/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_mem_opth3 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_3/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_mem_opth4 <- read.csv("results/CAMELIA_mem/mlp1_bce_msle/fold_4/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_mem_opth <- rbind(camelia_mem_opth1, camelia_mem_opth2, camelia_mem_opth3, camelia_mem_opth4)
camelia_mem_opth <- camelia_mem_opth[rownames(survey_obs), ]

camelia_mem_prr <- prr(camelia_mem_raw, SR_mem)

### camelia proxy
camelia_proxy_raw1 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_1/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_proxy_raw2 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_2/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_proxy_raw3 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_3/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_proxy_raw4 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_4/prediction_raw.csv", header = TRUE, row.names = 1)
camelia_proxy_raw <- rbind(camelia_proxy_raw1, camelia_proxy_raw2, camelia_proxy_raw3, camelia_proxy_raw4)
camelia_proxy_raw <- camelia_proxy_raw[rownames(survey_obs), ]

camelia_proxy_opth1 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_1/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_proxy_opth2 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_2/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_proxy_opth3 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_3/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_proxy_opth4 <- read.csv("results/CAMELIA_proxy/mlp1_bce_msle/fold_4/prediction_opth.csv", header = TRUE, row.names = 1)
camelia_proxy_opth <- rbind(camelia_proxy_opth1, camelia_proxy_opth2, camelia_proxy_opth3, camelia_proxy_opth4)
camelia_proxy_opth <- camelia_proxy_opth[rownames(survey_obs), ]

camelia_proxy_prr <- prr(camelia_proxy_raw, SR_mem)

# COMMUNITY INDICES
## OBS
CM_obs <- CM(survey_obs)
colnames(CM_obs)=c("ID","LNC_obs", "SLA_obs","PLH_obs")
CSTD_obs <- CSTD(survey_obs)
colnames(CSTD_obs)=c("ID","LNC_obs", "SLA_obs","PLH_obs")
SR_obs <- data.frame("ID"=rownames(survey_obs), "SR_obs"=rowSums(survey_obs))

## SSDM
CM_sdm_raw <- CM(sdm_raw)
colnames(CM_sdm_raw)=c("ID","LNC_sdm_raw", "SLA_sdm_raw","PLH_sdm_raw")
CSTD_sdm_raw <- CSTD(sdm_raw)
colnames(CSTD_sdm_raw)=c("ID","LNC_sdm_raw", "SLA_sdm_raw","PLH_sdm_raw")
SR_sdm_raw <- data.frame("ID"=rownames(sdm_raw), "SR_sdm_raw"=rowSums(sdm_raw))

CM_sdm_opth <- CM(sdm_opth)
colnames(CM_sdm_opth)=c("ID","LNC_sdm_opth", "SLA_sdm_opth","PLH_sdm_opth")
CSTD_sdm_opth <- CSTD(sdm_opth)
colnames(CSTD_sdm_opth)=c("ID","LNC_sdm_opth", "SLA_sdm_opth","PLH_sdm_opth")
SR_sdm_opth <- data.frame("ID"=rownames(sdm_opth), "SR_sdm_opth"=rowSums(sdm_opth))

CM_sdm_prr <- CM(sdm_prr)
colnames(CM_sdm_prr)=c("ID","LNC_sdm_prr", "SLA_sdm_prr","PLH_sdm_prr")
CSTD_sdm_prr <- CSTD(sdm_prr)
colnames(CSTD_sdm_prr)=c("ID","LNC_sdm_prr", "SLA_sdm_prr","PLH_sdm_prr")
SR_sdm_prr <- data.frame("ID"=rownames(sdm_prr), "SR_sdm_prr"=rowSums(sdm_prr))

## MEM
CM_mem <- read.csv("data/community_indices/CM_mem.csv", header = TRUE, row.names = 1)
colnames(CM_mem) <- c('LNC_mem', 'SLA_mem','PLH_mem')
CM_mem$ID <- rownames(CM_mem)
CSTD_mem <- read.csv("data/community_indices/CSTD_mem.csv", header = TRUE, row.names = 1)
colnames(CSTD_mem) <- c('LNC_mem', 'SLA_mem','PLH_mem')
CSTD_mem$ID <- rownames(CSTD_mem)
SR_mem <- read.csv("data/community_indices/SR_mem.csv", header = TRUE, row.names = 1)
colnames(SR_mem) <- c('SR_mem')
SR_mem$ID <- rownames(SR_mem)

## PROXY
CM_proxy <- read.csv("data/community_indices/CM_proxy.csv", header = TRUE, row.names = 1)
colnames(CM_proxy) <- c('LNC_proxy', 'SLA_proxy','PLH_proxy')
CM_proxy$ID <- rownames(CM_proxy)
CSTD_proxy <- read.csv("data/community_indices/CSTD_proxy.csv", header = TRUE, row.names = 1)
colnames(CSTD_proxy) <- c('LNC_proxy', 'SLA_proxy','PLH_proxy')
CSTD_proxy$ID <- rownames(CSTD_proxy)
SR_proxy <- read.csv("data/community_indices/SR_proxy.csv", header = TRUE, row.names = 1)
colnames(SR_proxy) <- c('SR_proxy')
SR_proxy$ID <- rownames(SR_proxy)

## CAMELIA OBS
CM_camelia_obs_raw <- CM(camelia_obs_raw)
colnames(CM_camelia_obs_raw)=c("ID","LNC_camelia_obs_raw", "SLA_camelia_obs_raw","PLH_camelia_obs_raw")
CSTD_camelia_obs_raw <- CSTD(camelia_obs_raw)
colnames(CSTD_camelia_obs_raw)=c("ID","LNC_camelia_obs_raw", "SLA_camelia_obs_raw","PLH_camelia_obs_raw")
SR_camelia_obs_raw <- data.frame("ID"=rownames(camelia_obs_raw), "SR_camelia_obs_raw"=rowSums(camelia_obs_raw))

CM_camelia_obs_opth <- CM(camelia_obs_opth)
colnames(CM_camelia_obs_opth)=c("ID","LNC_camelia_obs_opth", "SLA_camelia_obs_opth","PLH_camelia_obs_opth")
CSTD_camelia_obs_opth <- CSTD(camelia_obs_opth)
colnames(CSTD_camelia_obs_opth)=c("ID","LNC_camelia_obs_opth", "SLA_camelia_obs_opth","PLH_camelia_obs_opth")
SR_camelia_obs_opth <- data.frame("ID"=rownames(camelia_obs_opth), "SR_camelia_obs_opth"=rowSums(camelia_obs_opth))

CM_camelia_obs_prr <- CM(camelia_obs_prr)
colnames(CM_camelia_obs_prr)=c("ID","LNC_camelia_obs_prr", "SLA_camelia_obs_prr","PLH_camelia_obs_prr")
CSTD_camelia_obs_prr <- CSTD(camelia_obs_prr)
colnames(CSTD_camelia_obs_prr)=c("ID","LNC_camelia_obs_prr", "SLA_camelia_obs_prr","PLH_camelia_obs_prr")
SR_camelia_obs_prr <- data.frame("ID"=rownames(camelia_obs_prr), "SR_camelia_obs_prr"=rowSums(camelia_obs_prr))

## CAMELIA MEM
CM_camelia_mem_raw <- CM(camelia_mem_raw)
colnames(CM_camelia_mem_raw)=c("ID","LNC_camelia_mem_raw", "SLA_camelia_mem_raw","PLH_camelia_mem_raw")
CSTD_camelia_mem_raw <- CSTD(camelia_mem_raw)
colnames(CSTD_camelia_mem_raw)=c("ID","LNC_camelia_mem_raw", "SLA_camelia_mem_raw","PLH_camelia_mem_raw")
SR_camelia_mem_raw <- data.frame("ID"=rownames(camelia_mem_raw), "SR_camelia_mem_raw"=rowSums(camelia_mem_raw))

CM_camelia_mem_opth <- CM(camelia_mem_opth)
colnames(CM_camelia_mem_opth)=c("ID","LNC_camelia_mem_opth", "SLA_camelia_mem_opth","PLH_camelia_mem_opth")
CSTD_camelia_mem_opth <- CSTD(camelia_mem_opth)
colnames(CSTD_camelia_mem_opth)=c("ID","LNC_camelia_mem_opth", "SLA_camelia_mem_opth","PLH_camelia_mem_opth")
SR_camelia_mem_opth <- data.frame("ID"=rownames(camelia_mem_opth), "SR_camelia_mem_opth"=rowSums(camelia_mem_opth))

CM_camelia_mem_prr <- CM(camelia_mem_prr)
colnames(CM_camelia_mem_prr)=c("ID","LNC_camelia_mem_prr", "SLA_camelia_mem_prr","PLH_camelia_mem_prr")
CSTD_camelia_mem_prr <- CSTD(camelia_mem_prr)
colnames(CSTD_camelia_mem_prr)=c("ID","LNC_camelia_mem_prr", "SLA_camelia_mem_prr","PLH_camelia_mem_prr")
SR_camelia_mem_prr <- data.frame("ID"=rownames(camelia_mem_prr), "SR_camelia_mem_prr"=rowSums(camelia_mem_prr))

## CAMELIA PROXY
CM_camelia_proxy_raw <- CM(camelia_proxy_raw)
colnames(CM_camelia_proxy_raw)=c("ID","LNC_camelia_proxy_raw", "SLA_camelia_proxy_raw","PLH_camelia_proxy_raw")
CSTD_camelia_proxy_raw <- CSTD(camelia_proxy_raw)
colnames(CSTD_camelia_proxy_raw)=c("ID","LNC_camelia_proxy_raw", "SLA_camelia_proxy_raw","PLH_camelia_proxy_raw")
SR_camelia_proxy_raw <- data.frame("ID"=rownames(camelia_proxy_raw), "SR_camelia_proxy_raw"=rowSums(camelia_proxy_raw))

CM_camelia_proxy_opth <- CM(camelia_proxy_opth)
colnames(CM_camelia_proxy_opth)=c("ID","LNC_camelia_proxy_opth", "SLA_camelia_proxy_opth","PLH_camelia_proxy_opth")
CSTD_camelia_proxy_opth <- CSTD(camelia_proxy_opth)
colnames(CSTD_camelia_proxy_opth)=c("ID","LNC_camelia_proxy_opth", "SLA_camelia_proxy_opth","PLH_camelia_proxy_opth")
SR_camelia_proxy_opth <- data.frame("ID"=rownames(camelia_proxy_opth), "SR_camelia_proxy_opth"=rowSums(camelia_proxy_opth))

CM_camelia_proxy_prr <- CM(camelia_proxy_prr)
colnames(CM_camelia_proxy_prr)=c("ID","LNC_camelia_proxy_prr", "SLA_camelia_proxy_prr","PLH_camelia_proxy_prr")
CSTD_camelia_proxy_prr <- CSTD(camelia_proxy_prr)
colnames(CSTD_camelia_proxy_prr)=c("ID","LNC_camelia_proxy_prr", "SLA_camelia_proxy_prr","PLH_camelia_proxy_prr")
SR_camelia_proxy_prr <- data.frame("ID"=rownames(camelia_proxy_prr), "SR_camelia_proxy_prr"=rowSums(camelia_proxy_prr))

## TOTAL
CM_tot_raw <- left_join(CM_obs, CM_sdm_raw, by="ID")
CM_tot_raw <- left_join(CM_tot_raw, CM_camelia_mem_raw, by="ID")
CM_tot_raw <- left_join(CM_tot_raw, CM_camelia_obs_raw, by="ID")
CM_tot_raw <- left_join(CM_tot_raw, CM_camelia_proxy_raw, by="ID")
CM_tot_raw <- left_join(CM_tot_raw, CM_mem, by="ID")

CSTD_tot_raw <- left_join(CSTD_obs, CSTD_sdm_raw, by="ID")
CSTD_tot_raw <- left_join(CSTD_tot_raw, CSTD_camelia_mem_raw, by="ID")
CSTD_tot_raw <- left_join(CSTD_tot_raw, CSTD_camelia_obs_raw, by="ID")
CSTD_tot_raw <- left_join(CSTD_tot_raw, CSTD_camelia_proxy_raw, by="ID")
CSTD_tot_raw <- left_join(CSTD_tot_raw, CSTD_mem, by="ID")

SR_tot_raw <- left_join(SR_obs, SR_sdm_raw, by="ID")
SR_tot_raw <- left_join(SR_tot_raw, SR_camelia_mem_raw, by="ID")
SR_tot_raw <- left_join(SR_tot_raw, SR_camelia_obs_raw, by="ID")
SR_tot_raw <- left_join(SR_tot_raw, SR_camelia_proxy_raw, by="ID")
SR_tot_raw <- left_join(SR_tot_raw, SR_mem, by="ID")

CM_tot_opth <- left_join(CM_obs, CM_sdm_opth, by="ID")
CM_tot_opth <- left_join(CM_tot_opth, CM_camelia_mem_opth, by="ID")
CM_tot_opth <- left_join(CM_tot_opth, CM_camelia_obs_opth, by="ID")
CM_tot_opth <- left_join(CM_tot_opth, CM_camelia_proxy_opth, by="ID")
CM_tot_opth <- left_join(CM_tot_opth, CM_mem, by="ID")

CSTD_tot_opth <- left_join(CSTD_obs, CSTD_sdm_opth, by="ID")
CSTD_tot_opth <- left_join(CSTD_tot_opth, CSTD_camelia_mem_opth, by="ID")
CSTD_tot_opth <- left_join(CSTD_tot_opth, CSTD_camelia_obs_opth, by="ID")
CSTD_tot_opth <- left_join(CSTD_tot_opth, CSTD_camelia_proxy_opth, by="ID")
CSTD_tot_opth <- left_join(CSTD_tot_opth, CSTD_mem, by="ID")

SR_tot_opth <- left_join(SR_obs, SR_sdm_opth, by="ID")
SR_tot_opth <- left_join(SR_tot_opth, SR_camelia_mem_opth, by="ID")
SR_tot_opth <- left_join(SR_tot_opth, SR_camelia_obs_opth, by="ID")
SR_tot_opth <- left_join(SR_tot_opth, SR_camelia_proxy_opth, by="ID")
SR_tot_opth <- left_join(SR_tot_opth, SR_mem, by="ID")

CM_tot_prr <- left_join(CM_obs, CM_sdm_prr, by="ID")
CM_tot_prr <- left_join(CM_tot_prr, CM_camelia_mem_prr, by="ID")
CM_tot_prr <- left_join(CM_tot_prr, CM_camelia_obs_prr, by="ID")
CM_tot_prr <- left_join(CM_tot_prr, CM_camelia_proxy_prr, by="ID")
CM_tot_prr <- left_join(CM_tot_prr, CM_mem, by="ID")

CSTD_tot_prr <- left_join(CSTD_obs, CSTD_sdm_prr, by="ID")
CSTD_tot_prr <- left_join(CSTD_tot_prr, CSTD_camelia_mem_prr, by="ID")
CSTD_tot_prr <- left_join(CSTD_tot_prr, CSTD_camelia_obs_prr, by="ID")
CSTD_tot_prr <- left_join(CSTD_tot_prr, CSTD_camelia_proxy_prr, by="ID")
CSTD_tot_prr <- left_join(CSTD_tot_prr, CSTD_mem, by="ID")

SR_tot_prr <- left_join(SR_obs, SR_sdm_prr, by="ID")
SR_tot_prr <- left_join(SR_tot_prr, SR_camelia_mem_prr, by="ID")
SR_tot_prr <- left_join(SR_tot_prr, SR_camelia_obs_prr, by="ID")
SR_tot_prr <- left_join(SR_tot_prr, SR_camelia_proxy_prr, by="ID")
SR_tot_prr <- left_join(SR_tot_prr, SR_mem, by="ID")


#### -- FIGURE 1 -- ####
figure_raw <- data.frame("obs_CM"=rep(CM_tot_raw$PLH_obs,5),
                         "Pred_CM"=c(CM_tot_raw$PLH_sdm_raw,
                                     CM_tot_raw$PLH_mem,
                                     CM_tot_raw$PLH_camelia_obs_raw,
                                     CM_tot_raw$PLH_camelia_mem_raw,
                                     CM_tot_raw$PLH_camelia_proxy_raw),
                         'obs_CSTD'=rep(CSTD_tot_raw$PLH_obs,5),
                         "Pred_CSTD"=c(CSTD_tot_raw$PLH_sdm_raw,
                                       CSTD_tot_raw$PLH_mem,
                                       CSTD_tot_raw$PLH_camelia_obs_raw,
                                       CSTD_tot_raw$PLH_camelia_mem_raw,
                                       CSTD_tot_raw$PLH_camelia_proxy_raw),
                         'obs_SR' =rep(SR_tot_raw$SR_obs,5),
                         'Pred_SR' =c(SR_tot_raw$SR_sdm_raw,
                                      SR_tot_raw$SR_mem,
                                      SR_tot_raw$SR_camelia_obs_raw,
                                      SR_tot_raw$SR_camelia_mem_raw,
                                      SR_tot_raw$SR_camelia_proxy_raw),
                         "Models"=c(rep("SSDM", 4463), rep('MEM', 4463), rep("CAMELIA OBS", 4463), rep("CAMELIA MEM", 4463), rep("CAMELIA PROXY", 4463))
)

figure_raw$Models <- factor(figure_raw$Models, levels = c("SSDM", "MEM", "CAMELIA MEM", "CAMELIA PROXY", "CAMELIA OBS"))

cor_metrics <- data.frame(
  Models = factor(c("SSDM", "MEM", "CAMELIA MEM", "CAMELIA PROXY", "CAMELIA OBS"),
                  levels = c("SSDM", "MEM", "CAMELIA MEM", "CAMELIA PROXY", "CAMELIA OBS")),
  R2_CM = c(cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_sdm_raw)^2,
            cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_mem)^2,
            cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_mem_raw)^2,
            cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_proxy_raw)^2,
            cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_obs_raw)^2),
  RMSE_CM = c(rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_sdm_raw),
              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_mem),
              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_mem_raw),
              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_proxy_raw),
              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_obs_raw)),
  R2_CSTD = c(cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_sdm_raw)^2,
              cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_mem)^2,
              cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_mem_raw)^2,
              cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_proxy_raw)^2,
              cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_obs_raw)^2),
  RMSE_CSTD = c(rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_sdm_raw),
                rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_mem),
                rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_mem_raw),
                rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_proxy_raw),
                rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_obs_raw)),
  R2_SR = c(cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_sdm_raw)^2,
            cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_mem)^2,
            cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_mem_raw)^2,
            cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_proxy_raw)^2,
            cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_obs_raw)^2),
  RMSE_SR = c(rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_sdm_raw),
              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_mem),
              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_mem_raw),
              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_proxy_raw),
              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_obs_raw))
)


CM_raw <- ggplot(figure_raw, aes(x = obs_CM, y = Pred_CM, color = Models, fill = Models)) +
  geom_hdr() +
  geom_point(shape = 21, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("#E1812C", "#3274A1", "darkblue", "#B22222", "darkgreen")) +
  scale_fill_manual(values = c("#E1812C", "#3274A1", "darkblue", "#B22222", "darkgreen")) +
  facet_wrap(vars(Models), nrow = 1) +
  xlab("Observed CM PLH") +
  ylab("Predicted CM PLH") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 10)
  )+
  geom_label(data = cor_metrics,
             aes(x = -Inf, y = Inf,
                 label = sprintf("R² = %.2f\nRMSE = %.2f", round(R2_CM, 2), round(RMSE_CM, 2))),
             color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)

CSTD_raw <- ggplot(figure_raw, aes(x = obs_CSTD, y = Pred_CSTD, color = Models, fill = Models)) +
  geom_hdr() +
  geom_point(shape = 21, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("#E1812C", "#3274A1", "darkblue", "#B22222", "darkgreen")) +
  scale_fill_manual(values = c("#E1812C", "#3274A1", "darkblue", "#B22222", "darkgreen")) +
  facet_wrap(vars(Models), nrow = 1) +
  xlab("Observed CSTD PLH") +
  ylab("Predicted CSTD PLH") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 10)
  ) +
  geom_label(data = cor_metrics,
             aes(x = -Inf, y = Inf,
                 label = sprintf("R² = %.2f\nRMSE = %.2f", round(R2_CSTD, 2), round(RMSE_CSTD, 2))),
             color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)

SR_raw <- ggplot(figure_raw, aes(x = obs_SR, y = Pred_SR, color = Models, fill = Models)) +
  geom_hdr() +
  geom_point(shape = 21, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("#E1812C", "#3274A1", "darkblue", "#B22222", "darkgreen")) +
  scale_fill_manual(values = c("#E1812C", "#3274A1", "darkblue", "#B22222", "darkgreen")) +
  facet_wrap(vars(Models), nrow = 1) +
  xlab("Observed SR") +
  ylab("Predicted SR") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 10)
  ) +
  geom_label(data = cor_metrics,
             aes(x = -Inf, y = Inf,
                 label = sprintf("R² = %.2f\nRMSE = %.2f", round(R2_SR, 2), round(RMSE_SR, 2))),
             color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)

combined_plot <- SR_raw /  CM_raw / CSTD_raw +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("figure_1.pdf", 
       plot = combined_plot, 
       device = "pdf", 
       width = 8.27,
       height = 11.69,
       dpi = 300)


#### -- FIGURE 1 -- ####
## ranking based binarization
eval_rbb <- data.frame("ID"=rep(0, 100*4463),
                       "richness"= rep(0, 100*4463),
                       "CM_LNC"= rep(0, 100*4463),
                       "CM_SLA"= rep(0, 100*4463),
                       "CM_PLH"= rep(0, 100*4463),
                       "CSTD_LNC"= rep(0, 100*4463),
                       "CSTD_SLA"= rep(0, 100*4463),
                       "CSTD_PLH"= rep(0, 100*4463),
                       "richness_ssdm"=rep(0, 100*4463),
                       "richness_camelia_mem"=rep(0, 100*4463), 
                       "richness_camelia_obs"=rep(0, 100*4463),
                       "richness_camelia_proxy"=rep(0, 100*4463),
                       "Precision_ssdm"=rep(0, 100*4463),
                       "Precision_camelia_mem"=rep(0, 100*4463),
                       "Precision_camelia_obs"=rep(0, 100*4463),
                       "Precision_camelia_proxy"=rep(0, 100*4463),
                       "Sorensen_ssdm"=rep(0, 100*4463),
                       "Sorensen_camelia_mem"=rep(0, 100*4463),
                       "Sorensen_camelia_obs"=rep(0, 100*4463),
                       "Sorensen_camelia_proxy"=rep(0, 100*4463),
                       "Recall_ssdm"=rep(0,100*4463),
                       "Recall_camelia_mem"=rep(0, 100*4463),
                       "Recall_camelia_obs"=rep(0, 100*4463),
                       "Recall_camelia_proxy"=rep(0, 100*4463),
                       "Specificite_ssdm"=rep(0,100*4463),
                       "Specificite_camelia_mem"=rep(0, 100*4463),
                       "Specificite_camelia_obs"=rep(0, 100*4463),
                       "Specificite_camelia_proxy"=rep(0, 100*4463),
                       "CM_LNC_ssdm"=rep(0, 100*4463),
                       "CM_SLA_ssdm"=rep(0, 100*4463),
                       "CM_PLH_ssdm"=rep(0, 100*4463),
                       "CM_LNC_camelia_mem"=rep(0, 100*4463),
                       "CM_SLA_camelia_mem"=rep(0, 100*4463),
                       "CM_PLH_camelia_mem"=rep(0, 100*4463),
                       "CM_LNC_camelia_obs"=rep(0, 100*4463),
                       "CM_SLA_camelia_obs"=rep(0, 100*4463),
                       "CM_PLH_camelia_obs"=rep(0, 100*4463),
                       "CM_LNC_camelia_proxy"=rep(0, 100*4463),
                       "CM_SLA_camelia_proxy"=rep(0, 100*4463),
                       "CM_PLH_camelia_proxy"=rep(0, 100*4463),
                       "CSTD_LNC_ssdm"=rep(0, 100*4463),
                       "CSTD_SLA_ssdm"=rep(0, 100*4463),
                       "CSTD_PLH_ssdm"=rep(0, 100*4463),
                       "CSTD_LNC_camelia_mem"=rep(0, 100*4463),
                       "CSTD_SLA_camelia_mem"=rep(0, 100*4463),
                       "CSTD_PLH_camelia_mem"=rep(0, 100*4463),
                       "CSTD_LNC_camelia_obs"=rep(0, 100*4463),
                       "CSTD_SLA_camelia_obs"=rep(0, 100*4463),
                       "CSTD_PLH_camelia_obs"=rep(0, 100*4463),
                       "CSTD_LNC_camelia_proxy"=rep(0, 100*4463),
                       "CSTD_SLA_camelia_proxy"=rep(0, 100*4463),
                       "CSTD_PLH_camelia_proxy"=rep(0, 100*4463))
for (i in 1:4463){
  obs<- as.numeric(survey_obs[i,])
  eval_rbb$ID[((i-1)*100+1):((i-1)*100+100)]=rownames(survey_obs)[i]
  eval_rbb$SR[((i-1)*100+1):((i-1)*100+100)]=length(which(obs==1))
  eval_rbb$CM_LNC[((i-1)*100+1):((i-1)*100+100)]=sum(obs*traits$LNC, na.rm = TRUE)/length(which(obs==1))
  eval_rbb$CM_SLA[((i-1)*100+1):((i-1)*100+100)]=sum(obs*traits$SLA, na.rm = TRUE)/length(which(obs==1))
  eval_rbb$CM_PLH[((i-1)*100+1):((i-1)*100+100)]=sum(obs*traits$PLH, na.rm = TRUE)/length(which(obs==1))
  
  eval_rbb$CSTD_LNC[((i-1)*100+1):((i-1)*100+100)]= sqrt(sum(obs*(traits$LNC-eval_rbb$CM_LNC[(i-1)*100+1])^2, na.rm=TRUE)/length(which(obs==1)))
  eval_rbb$CSTD_SLA[((i-1)*100+1):((i-1)*100+100)]= sqrt(sum(obs*(traits$SLA-eval_rbb$CM_SLA[(i-1)*100+1])^2, na.rm=TRUE)/length(which(obs==1)))
  eval_rbb$CSTD_PLH[((i-1)*100+1):((i-1)*100+100)]= sqrt(sum(obs*(traits$PLH-eval_rbb$CM_PLH[(i-1)*100+1])^2, na.rm=TRUE)/length(which(obs==1)))
  
  predictions_sdm <- as.numeric(ssdm_predictions[i,])
  predictions_camelia_mem <- as.numeric(camelia_mem_raw[i,])
  predictions_camelia_obs <- as.numeric(camelia_obs_raw[i,])
  predictions_camelia_proxy <- as.numeric(camelia_proxy_raw[i,])
  for (j in seq(0.01, 1, by = 0.01)){
    predictions_sdm_j <- rep(0,831)
    predictions_sdm_j[which(predictions_sdm>j)]=1
    predictions_camelia_mem_j <- rep(0,831)
    predictions_camelia_mem_j[which(predictions_camelia_mem>j)]=1
    predictions_camelia_obs_j <- rep(0,831)
    predictions_camelia_obs_j[which(predictions_camelia_obs>j)]=1
    predictions_camelia_proxy_j <- rep(0,831)
    predictions_camelia_proxy_j[which(predictions_camelia_proxy>j)]=1
    
    VP_sdm = length(which(which(predictions_sdm_j==1) %in% which(obs==1)))
    FP_sdm = length(which(which(predictions_sdm_j==1) %in% which(obs==0)))
    FN_sdm = length(which(which(predictions_sdm_j==0) %in% which(obs==1)))
    VN_sdm = length(which(which(predictions_sdm_j==0) %in% which(obs==0)))
    eval_rbb$Precision_ssdm[(i-1)*100+j*100] = VP_sdm/(VP_sdm+FP_sdm)
    eval_rbb$Sorensen_ssdm[(i-1)*100+j*100] = 2*VP_sdm/(2*VP_sdm+FP_sdm+FN_sdm)
    eval_rbb$Recall_ssdm[(i-1)*100+j*100] = VP_sdm/(VP_sdm+FN_sdm)
    eval_rbb$Specificite_ssdm[(i-1)*100+j*100] = FP_sdm/(VN_sdm+FP_sdm)
    eval_rbb$richness_ssdm[(i-1)*100+j*100]=length(which(predictions_sdm_j==1))
    eval_rbb$CM_LNC_ssdm[(i-1)*100+j*100]=sum(predictions_sdm_j*traits$LNC, na.rm = TRUE)/length(which(predictions_sdm_j==1))
    eval_rbb$CM_SLA_ssdm[(i-1)*100+j*100]=sum(predictions_sdm_j*traits$SLA, na.rm = TRUE)/length(which(predictions_sdm_j==1))
    eval_rbb$CM_PLH_ssdm[(i-1)*100+j*100]=sum(predictions_sdm_j*traits$PLH, na.rm = TRUE)/length(which(predictions_sdm_j==1))
    eval_rbb$CSTD_LNC_ssdm[(i-1)*100+j*100]= sqrt(sum(predictions_sdm_j*(traits$LNC-eval_rbb$CM_LNC_ssdm[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_sdm_j==1)))
    eval_rbb$CSTD_SLA_ssdm[(i-1)*100+j*100]= sqrt(sum(predictions_sdm_j*(traits$SLA-eval_rbb$CM_SLA_ssdm[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_sdm_j==1)))
    eval_rbb$CSTD_PLH_ssdm[(i-1)*100+j*100]= sqrt(sum(predictions_sdm_j*(traits$PLH-eval_rbb$CM_PLH_ssdm[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_sdm_j==1)))
    
    VP_camelia_mem = length(which(which(predictions_camelia_mem_j==1) %in% which(obs==1)))
    FP_camelia_mem = length(which(which(predictions_camelia_mem_j==1) %in% which(obs==0)))
    FN_camelia_mem = length(which(which(predictions_camelia_mem_j==0) %in% which(obs==1)))
    VN_camelia_mem = length(which(which(predictions_camelia_mem_j==0) %in% which(obs==0)))
    eval_rbb$Precision_camelia_mem[(i-1)*100+j*100] = VP_camelia_mem/(VP_camelia_mem+FP_camelia_mem)
    eval_rbb$Recall_camelia_mem[(i-1)*100+j*100] = VP_camelia_mem/(VP_camelia_mem+FN_camelia_mem)
    eval_rbb$Specificite_camelia_mem[(i-1)*100+j*100] = FP_camelia_mem/(VN_camelia_mem+FP_camelia_mem)
    eval_rbb$Sorensen_camelia_mem[(i-1)*100+j*100] = 2*VP_camelia_mem/(2*VP_camelia_mem+FP_camelia_mem+FN_camelia_mem)
    eval_rbb$richness_camelia_mem[(i-1)*100+j*100]=length(which(predictions_camelia_mem_j==1))
    eval_rbb$CM_LNC_camelia_mem[(i-1)*100+j*100]=sum(predictions_camelia_mem_j*traits$LNC, na.rm = TRUE)/length(which(predictions_camelia_mem_j==1))
    eval_rbb$CM_SLA_camelia_mem[(i-1)*100+j*100]=sum(predictions_camelia_mem_j*traits$SLA, na.rm = TRUE)/length(which(predictions_camelia_mem_j==1))
    eval_rbb$CM_PLH_camelia_mem[(i-1)*100+j*100]=sum(predictions_camelia_mem_j*traits$PLH, na.rm = TRUE)/length(which(predictions_camelia_mem_j==1))
    eval_rbb$CSTD_LNC_camelia_mem[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_mem_j*(traits$LNC-eval_rbb$CM_LNC_camelia_mem[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_mem_j==1)))
    eval_rbb$CSTD_SLA_camelia_mem[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_mem_j*(traits$SLA-eval_rbb$CM_SLA_camelia_mem[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_mem_j==1)))
    eval_rbb$CSTD_PLH_camelia_mem[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_mem_j*(traits$PLH-eval_rbb$CM_PLH_camelia_mem[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_mem_j==1)))
    
    VP_camelia_obs = length(which(which(predictions_camelia_obs_j==1) %in% which(obs==1)))
    FP_camelia_obs = length(which(which(predictions_camelia_obs_j==1) %in% which(obs==0)))
    FN_camelia_obs = length(which(which(predictions_camelia_obs_j==0) %in% which(obs==1)))
    VN_camelia_obs = length(which(which(predictions_camelia_obs_j==0) %in% which(obs==0)))
    eval_rbb$Precision_camelia_obs[(i-1)*100+j*100] = VP_camelia_obs/(VP_camelia_obs+FP_camelia_obs)
    eval_rbb$Sorensen_camelia_obs[(i-1)*100+j*100] = 2*VP_camelia_obs/(2*VP_camelia_obs+FP_camelia_obs+FN_camelia_obs)
    eval_rbb$Recall_camelia_obs[(i-1)*100+j*100] = VP_camelia_obs/(VP_camelia_obs+FN_camelia_obs)
    eval_rbb$Specificite_camelia_obs[(i-1)*100+j*100] = FP_camelia_obs/(VN_camelia_obs+FP_camelia_obs)
    eval_rbb$richness_camelia_obs[(i-1)*100+j*100]=length(which(predictions_camelia_obs_j==1))
    eval_rbb$CM_LNC_camelia_obs[(i-1)*100+j*100]=sum(predictions_camelia_obs_j*traits$LNC, na.rm = TRUE)/length(which(predictions_camelia_obs_j==1))
    eval_rbb$CM_SLA_camelia_obs[(i-1)*100+j*100]=sum(predictions_camelia_obs_j*traits$SLA, na.rm = TRUE)/length(which(predictions_camelia_obs_j==1))
    eval_rbb$CM_PLH_camelia_obs[(i-1)*100+j*100]=sum(predictions_camelia_obs_j*traits$PLH, na.rm = TRUE)/length(which(predictions_camelia_obs_j==1))
    eval_rbb$CSTD_LNC_camelia_obs[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_obs_j*(traits$LNC-eval_rbb$CM_LNC_camelia_obs[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_obs_j==1)))
    eval_rbb$CSTD_SLA_camelia_obs[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_obs_j*(traits$SLA-eval_rbb$CM_SLA_camelia_obs[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_obs_j==1)))
    eval_rbb$CSTD_PLH_camelia_obs[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_obs_j*(traits$PLH-eval_rbb$CM_PLH_camelia_obs[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_obs_j==1)))
    
    VP_camelia_proxy = length(which(which(predictions_camelia_proxy_j==1) %in% which(obs==1)))
    FP_camelia_proxy = length(which(which(predictions_camelia_proxy_j==1) %in% which(obs==0)))
    FN_camelia_proxy = length(which(which(predictions_camelia_proxy_j==0) %in% which(obs==1)))
    VN_camelia_proxy = length(which(which(predictions_camelia_proxy_j==0) %in% which(obs==0)))
    eval_rbb$Precision_camelia_proxy[(i-1)*100+j*100] = VP_camelia_proxy/(VP_camelia_proxy+FP_camelia_proxy)
    eval_rbb$Sorensen_camelia_proxy[(i-1)*100+j*100] = 2*VP_camelia_proxy/(2*VP_camelia_proxy+FP_camelia_proxy+FN_camelia_proxy)
    eval_rbb$Recall_camelia_proxy[(i-1)*100+j*100] = VP_camelia_proxy/(VP_camelia_proxy+FN_camelia_proxy)
    eval_rbb$Specificite_camelia_proxy[(i-1)*100+j*100] = FP_camelia_proxy/(VN_camelia_proxy+FP_camelia_proxy)
    eval_rbb$richness_camelia_proxy[(i-1)*100+j*100]=length(which(predictions_camelia_proxy_j==1))
    eval_rbb$CM_LNC_camelia_proxy[(i-1)*100+j*100]=sum(predictions_camelia_proxy_j*traits$LNC, na.rm = TRUE)/length(which(predictions_camelia_proxy_j==1))
    eval_rbb$CM_SLA_camelia_proxy[(i-1)*100+j*100]=sum(predictions_camelia_proxy_j*traits$SLA, na.rm = TRUE)/length(which(predictions_camelia_proxy_j==1))
    eval_rbb$CM_PLH_camelia_proxy[(i-1)*100+j*100]=sum(predictions_camelia_proxy_j*traits$PLH, na.rm = TRUE)/length(which(predictions_camelia_proxy_j==1))
    eval_rbb$CSTD_LNC_camelia_proxy[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_proxy_j*(traits$LNC-eval_rbb$CM_LNC_camelia_proxy[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_proxy_j==1)))
    eval_rbb$CSTD_SLA_camelia_proxy[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_proxy_j*(traits$SLA-eval_rbb$CM_SLA_camelia_proxy[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_proxy_j==1)))
    eval_rbb$CSTD_PLH_camelia_proxy[(i-1)*100+j*100]= sqrt(sum(predictions_camelia_proxy_j*(traits$PLH-eval_rbb$CM_PLH_camelia_proxy[(i-1)*100+j*100])^2, na.rm=TRUE)/length(which(predictions_camelia_proxy_j==1)))
  }
  print(i)
}

plot_diff_CSTD_PLH <- ggplot(eval_rbb) + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = 'gray', alpha=0.3) +
  stat_summary(aes(x=log(eval_rbb$richness_ssdm/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_ssdm)/eval_rbb$CSTD_PLH)),fun.y = "mean", geom = "line", color = "#E1812C", alpha=0.2) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_mem/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_camelia_mem)/eval_rbb$CSTD_PLH)),fun.y = "mean", geom = "line", color="darkblue", alpha=0.2) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_camelia_obs)/eval_rbb$CSTD_PLH)),fun.y = "mean", geom = "line", color="darkgreen", alpha=0.2) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_proxy/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_camelia_proxy)/eval_rbb$CSTD_PLH)),fun.y = "mean", geom = "line", color="#B22222", alpha=0.2) +
  
  geom_smooth(aes(x=log(eval_rbb$richness_ssdm/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_ssdm)/eval_rbb$CSTD_PLH)), span = 0.05, color = "#E1812C", fill = "#E1812C") + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_mem/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_camelia_mem)/eval_rbb$CSTD_PLH)),span = 0.05, color="darkblue", fill = "darkblue") + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_camelia_obs)/eval_rbb$CSTD_PLH)),span = 0.05, color="darkgreen", fill = "darkgreen") +  
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_proxy/eval_rbb$richness), y=abs((eval_rbb$CSTD_PLH-eval_rbb$CSTD_PLH_camelia_proxy)/eval_rbb$CSTD_PLH)),span = 0.05, color="#B22222", fill = "#B22222") +  
  scale_y_continuous(trans=log10_trans(), labels = scales::comma, limits = c(0.0000005, 10)) +
  coord_trans(y='reverse') +
  xlab(bquote('log(richness pred/richness obs)')) + 
  ylab(bquote('abs((CSTD obs - CSTD pred)/CSTD obs)')) + 
  theme_bw()

plot_diff_CM_PLH <- ggplot(eval_rbb) +   
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = 'gray', alpha=0.3) +
  stat_summary(aes(x=log(eval_rbb$richness_ssdm/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_ssdm)/eval_rbb$CM_PLH)),fun.y = "mean", geom = "line", color = "#E1812C", alpha=0.2) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_mem/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_camelia_mem)/eval_rbb$CM_PLH)),fun.y = "mean", geom = "line", color="darkblue", alpha=0.2) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_camelia_obs)/eval_rbb$CM_PLH)),fun.y = "mean", geom = "line", color="darkgreen", alpha=0.2) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_proxy/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_camelia_proxy)/eval_rbb$CM_PLH)),fun.y = "mean", geom = "line", color="#B22222", alpha=0.2) + 
  
  geom_smooth(aes(x=log(eval_rbb$richness_ssdm/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_ssdm)/eval_rbb$CM_PLH)), span = 0.05, color = "#E1812C", fill = "#E1812C") + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_mem/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_camelia_mem)/eval_rbb$CM_PLH)),span = 0.05, color="darkblue", fill = "darkblue") + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_camelia_obs)/eval_rbb$CM_PLH)),span = 0.05, color="darkgreen", fill = "darkgreen") +  
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_proxy/eval_rbb$richness), y=abs((eval_rbb$CM_PLH-eval_rbb$CM_PLH_camelia_proxy)/eval_rbb$CM_PLH)),span = 0.05, color="#B22222", fill = "#B22222") +  
  scale_y_continuous(trans = log10_trans(), labels = scales::comma, limits = c(0.0000005, 10))+ #, breaks = c(1e-4, 1e-2, 1)) +
  coord_trans(y='reverse') +
  xlab(bquote('log(richness pred/richness obs)')) + 
  ylab(bquote('abs((CM obs - CM pred)/CM obs)')) + 
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none")

Sorensen_logrichness <- ggplot(eval_rbb) + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = 'gray', alpha=0.3) +
  stat_summary(aes(x=log(eval_rbb$richness_ssdm/eval_rbb$richness), y=eval_rbb$Sorensen_ssdm),fun.y = "mean", geom = "line",color = "#E1812C", alpha=0.3) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_mem/eval_rbb$richness), y=eval_rbb$Sorensen_camelia_mem),fun.y = "mean", geom = "line", color="darkblue", alpha=0.3) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=eval_rbb$Sorensen_camelia_obs),fun.y = "mean", geom = "line", color="darkgreen", alpha=0.3) + 
  stat_summary(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=eval_rbb$Sorensen_camelia_proxy),fun.y = "mean", geom = "line", color="#B22222", alpha=0.3) + 
  
  geom_smooth(aes(x=log(eval_rbb$richness_ssdm/eval_rbb$richness), y=eval_rbb$Sorensen_ssdm), span = 0.05, color = "#E1812C", fill = "#E1812C", se=FALSE) + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_mem/eval_rbb$richness), y=eval_rbb$Sorensen_camelia_mem),span = 0.05, color="darkblue", fill = "darkblue", se=FALSE) + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=eval_rbb$Sorensen_camelia_obs),span = 0.05, color="darkgreen", fill = "darkgreen", se=FALSE) + 
  geom_smooth(aes(x=log(eval_rbb$richness_camelia_obs/eval_rbb$richness), y=eval_rbb$Sorensen_camelia_proxy),span = 0.05, color="#B22222", fill = "#B22222", se=FALSE) + 
  xlab(bquote('log(richness pred/richness obs)')) + 
  ylab(bquote('Sorensen')) + 
  ylim(-0.01, 1)+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none")

## TSS optimised threshold binarization :
figure_sites_opth <- data.frame("ID"=Sites$ID, 
                                "richness"=rowSums(survey_obs),
                                "richness_sdm_opth"=rowSums(sdm_opth),
                                "VP_sdm_opth"=rep(0,4463),
                                "FP_sdm_opth"=rep(0,4463),
                                "VN_sdm_opth"=rep(0,4463),
                                "FN_sdm_opth"=rep(0,4463),
                                "Specificite_sdm_opth"=rep(0,4463),
                                "Sensibilite_sdm_opth"=rep(0,4463),
                                "Precision_sdm_opth"=rep(0,4463),
                                "Sorensen_sdm_opth"=rep(0,4463),
                                "VP_camelia_mem_opth"=rep(0,4463),
                                "FP_camelia_mem_opth"=rep(0,4463),
                                "VN_camelia_mem_opth"=rep(0,4463),
                                "FN_camelia_mem_opth"=rep(0,4463),
                                "richness_camelia_mem_opth"=rowSums(camelia_mem_opth),
                                "Specificite_camelia_mem_opth"=rep(0,4463),
                                "Sensibilite_camelia_mem_opth"=rep(0,4463),
                                "Precision_camelia_mem_opth"=rep(0,4463),
                                "Sorensen_camelia_mem_opth"=rep(0,4463),
                                "VP_camelia_obs_opth"=rep(0,4463),
                                "FP_camelia_obs_opth"=rep(0,4463),
                                "VN_camelia_obs_opth"=rep(0,4463),
                                "FN_camelia_obs_opth"=rep(0,4463),
                                "richness_camelia_obs_opth"=rowSums(camelia_obs_opth),
                                "Specificite_camelia_obs_opth"=rep(0,4463),
                                "Sensibilite_camelia_obs_opth"=rep(0,4463),
                                "Precision_camelia_obs_opth"=rep(0,4463),
                                "Sorensen_camelia_obs_opth"=rep(0,4463),
                                "VP_camelia_proxy_opth"=rep(0,4463),
                                "FP_camelia_proxy_opth"=rep(0,4463),
                                "VN_camelia_proxy_opth"=rep(0,4463),
                                "FN_camelia_proxy_opth"=rep(0,4463),
                                "richness_camelia_proxy_opth"=rowSums(camelia_proxy_opth),
                                "Specificite_camelia_proxy_opth"=rep(0,4463),
                                "Sensibilite_camelia_proxy_opth"=rep(0,4463),
                                "Precision_camelia_proxy_opth"=rep(0,4463),
                                "Sorensen_camelia_proxy_opth"=rep(0,4463))

for (i in 1:4463){
  VP_sdm_opth=length(which(which(sdm_opth[i,]==1) %in% which(survey_obs[i,]==1)))
  FP_sdm_opth=length(which(which(sdm_opth[i,]==1) %in% which(survey_obs[i,]==0)))
  VN_sdm_opth=length(which(which(sdm_opth[i,]==0) %in% which(survey_obs[i,]==0)))
  FN_sdm_opth=length(which(which(sdm_opth[i,]==0) %in% which(survey_obs[i,]==1)))
  figure_sites_opth$VP_sdm_opth[i]=VP_sdm_opth
  figure_sites_opth$FP_sdm_opth[i]=FP_sdm_opth
  figure_sites_opth$VN_sdm_opth[i]=VN_sdm_opth
  figure_sites_opth$FN_sdm_opth[i]=FN_sdm_opth
  figure_sites_opth$Sensibilite_sdm_opth[i]=VP_sdm_opth/(VP_sdm_opth+FN_sdm_opth)
  figure_sites_opth$Specificite_sdm_opth[i]=VN_sdm_opth/(VN_sdm_opth+FP_sdm_opth)
  figure_sites_opth$Precision_sdm_opth[i]=VP_sdm_opth/(VP_sdm_opth+FP_sdm_opth)
  figure_sites_opth$Sorensen_sdm_opth[i]=2*VP_sdm_opth/(2*VP_sdm_opth+FP_sdm_opth+FN_sdm_opth)
  
  VP_camelia_mem_opth=length(which(which(camelia_mem_opth[i,]==1) %in% which(survey_obs[i,]==1)))
  FP_camelia_mem_opth=length(which(which(camelia_mem_opth[i,]==1) %in% which(survey_obs[i,]==0)))
  VN_camelia_mem_opth=length(which(which(camelia_mem_opth[i,]==0) %in% which(survey_obs[i,]==0)))
  FN_camelia_mem_opth=length(which(which(camelia_mem_opth[i,]==0) %in% which(survey_obs[i,]==1)))
  figure_sites_opth$VP_camelia_mem_opth[i]=VP_camelia_mem_opth
  figure_sites_opth$FP_camelia_mem_opth[i]=FP_camelia_mem_opth
  figure_sites_opth$VN_camelia_mem_opth[i]=VN_camelia_mem_opth
  figure_sites_opth$FN_camelia_mem_opth[i]=FN_camelia_mem_opth
  figure_sites_opth$Sensibilite_camelia_mem_opth[i]=VP_camelia_mem_opth/(VP_camelia_mem_opth+FN_camelia_mem_opth)
  figure_sites_opth$Specificite_camelia_mem_opth[i]=VN_camelia_mem_opth/(VN_camelia_mem_opth+FP_camelia_mem_opth)
  figure_sites_opth$Precision_camelia_mem_opth[i]=VP_camelia_mem_opth/(VP_camelia_mem_opth+FP_camelia_mem_opth)
  figure_sites_opth$Sorensen_camelia_mem_opth[i]=2*VP_camelia_mem_opth/(2*VP_camelia_mem_opth+FP_camelia_mem_opth+FN_camelia_mem_opth)
  
  VP_camelia_obs_opth=length(which(which(camelia_obs_opth[i,]==1) %in% which(survey_obs[i,]==1)))
  FP_camelia_obs_opth=length(which(which(camelia_obs_opth[i,]==1) %in% which(survey_obs[i,]==0)))
  VN_camelia_obs_opth=length(which(which(camelia_obs_opth[i,]==0) %in% which(survey_obs[i,]==0)))
  FN_camelia_obs_opth=length(which(which(camelia_obs_opth[i,]==0) %in% which(survey_obs[i,]==1)))
  figure_sites_opth$VP_camelia_obs_opth[i]=VP_camelia_obs_opth
  figure_sites_opth$FP_camelia_obs_opth[i]=FP_camelia_obs_opth
  figure_sites_opth$VN_camelia_obs_opth[i]=VN_camelia_obs_opth
  figure_sites_opth$FN_camelia_obs_opth[i]=FN_camelia_obs_opth
  figure_sites_opth$Sensibilite_camelia_obs_opth[i]=VP_camelia_obs_opth/(VP_camelia_obs_opth+FN_camelia_obs_opth)
  figure_sites_opth$Specificite_camelia_obs_opth[i]=VN_camelia_obs_opth/(VN_camelia_obs_opth+FP_camelia_obs_opth)
  figure_sites_opth$Precision_camelia_obs_opth[i]=VP_camelia_obs_opth/(VP_camelia_obs_opth+FP_camelia_obs_opth)
  figure_sites_opth$Sorensen_camelia_obs_opth[i]=2*VP_camelia_obs_opth/(2*VP_camelia_obs_opth+FP_camelia_obs_opth+FN_camelia_obs_opth)
  
  VP_camelia_proxy_opth=length(which(which(camelia_proxy_opth[i,]==1) %in% which(survey_obs[i,]==1)))
  FP_camelia_proxy_opth=length(which(which(camelia_proxy_opth[i,]==1) %in% which(survey_obs[i,]==0)))
  VN_camelia_proxy_opth=length(which(which(camelia_proxy_opth[i,]==0) %in% which(survey_obs[i,]==0)))
  FN_camelia_proxy_opth=length(which(which(camelia_proxy_opth[i,]==0) %in% which(survey_obs[i,]==1)))
  figure_sites_opth$VP_camelia_proxy_opth[i]=VP_camelia_proxy_opth
  figure_sites_opth$FP_camelia_proxy_opth[i]=FP_camelia_proxy_opth
  figure_sites_opth$VN_camelia_proxy_opth[i]=VN_camelia_proxy_opth
  figure_sites_opth$FN_camelia_proxy_opth[i]=FN_camelia_proxy_opth
  figure_sites_opth$Sensibilite_camelia_proxy_opth[i]=VP_camelia_proxy_opth/(VP_camelia_proxy_opth+FN_camelia_proxy_opth)
  figure_sites_opth$Specificite_camelia_proxy_opth[i]=VN_camelia_proxy_opth/(VN_camelia_proxy_opth+FP_camelia_proxy_opth)
  figure_sites_opth$Precision_camelia_proxy_opth[i]=VP_camelia_proxy_opth/(VP_camelia_proxy_opth+FP_camelia_proxy_opth)
  figure_sites_opth$Sorensen_camelia_proxy_opth[i]=2*VP_camelia_proxy_opth/(2*VP_camelia_proxy_opth+FP_camelia_proxy_opth+FN_camelia_proxy_opth)
  print(i)
}

figure_sites_opth <- left_join(figure_sites_opth, CM_tot_opth, by="ID")
figure_sites_opth <- left_join(figure_sites_opth, CSTD_tot_opth, by="ID")


data_richness_precision_opth=data.frame("ID"=rep(figure_sites_opth$ID,4),
                                        "log(rich_pred/rich_obs)"=c(log(figure_sites_opth$richness_sdm_opth/figure_sites_opth$richness),
                                                                    log(figure_sites_opth$richness_camelia_mem_opth/figure_sites_opth$richness),
                                                                    log(figure_sites_opth$richness_camelia_obs_opth/figure_sites_opth$richness),
                                                                    log(figure_sites_opth$richness_camelia_proxy_opth/figure_sites_opth$richness)),
                                        "Precision"=c(figure_sites_opth$Precision_sdm_opth, figure_sites_opth$Precision_camelia_mem_opth_opth, figure_sites_opth$Precision_camelia_obs_opth, figure_sites_opth$Precision_camelia_proxy_opth),
                                        "Sorensen"=c(figure_sites_opth$Sorensen_sdm_opth, figure_sites_opth$Sorensen_camelia_mem_opth_opth, figure_sites_opth$Sorensen_camelia_obs_opth, figure_sites_opth$Sorensen_camelia_proxy_opth),
                                        "PLH abs((CSTD obs - CSTD pred)/CSTD obs)"=c(abs(figure_sites_opth$CSTD_PLH_obs-figure_sites_opth$CSTD_PLH_sdm_opth)/figure_sites_opth$CSTD_PLH_obs,
                                                                                     abs(figure_sites_opth$CSTD_PLH_obs-figure_sites_opth$CSTD_PLH_camelia_mem_opth_opth)/figure_sites_opth$CSTD_PLH_obs,
                                                                                     abs(figure_sites_opth$CSTD_PLH_obs-figure_sites_opth$CSTD_PLH_camelia_obs_opth)/figure_sites_opth$CSTD_PLH_obs,
                                                                                     abs(figure_sites_opth$CSTD_PLH_obs-figure_sites_opth$CSTD_PLH_camelia_proxy_opth)/figure_sites_opth$CSTD_PLH_obs),
                                        "SLA abs((CSTD obs - CSTD pred)/CSTD obs)"=c(abs(figure_sites_opth$CSTD_SLA_obs-figure_sites_opth$CSTD_SLA_sdm_opth)/figure_sites_opth$CSTD_SLA_obs,
                                                                                     abs(figure_sites_opth$CSTD_SLA_obs-figure_sites_opth$CSTD_SLA_camelia_mem_opth_opth)/figure_sites_opth$CSTD_SLA_obs,
                                                                                     abs(figure_sites_opth$CSTD_SLA_obs-figure_sites_opth$CSTD_SLA_camelia_obs_opth)/figure_sites_opth$CSTD_SLA_obs,
                                                                                     abs(figure_sites_opth$CSTD_SLA_obs-figure_sites_opth$CSTD_SLA_camelia_proxy_opth)/figure_sites_opth$CSTD_SLA_obs),
                                        "LNC abs((CSTD obs - CSTD pred)/CSTD obs)"=c(abs(figure_sites_opth$CSTD_LNC_obs-figure_sites_opth$CSTD_LNC_sdm_opth)/figure_sites_opth$CSTD_LNC_obs,
                                                                                     abs(figure_sites_opth$CSTD_LNC_obs-figure_sites_opth$CSTD_LNC_camelia_mem_opth_opth)/figure_sites_opth$CSTD_LNC_obs,
                                                                                     abs(figure_sites_opth$CSTD_LNC_obs-figure_sites_opth$CSTD_LNC_camelia_obs_opth)/figure_sites_opth$CSTD_LNC_obs,
                                                                                     abs(figure_sites_opth$CSTD_LNC_obs-figure_sites_opth$CSTD_LNC_camelia_proxy_opth)/figure_sites_opth$CSTD_LNC_obs),
                                        "PLH abs((CM obs - CM pred)/CM obs)"=c(abs(figure_sites_opth$CM_PLH_obs-figure_sites_opth$CM_PLH_sdm_opth)/figure_sites_opth$CM_PLH_obs,
                                                                               abs(figure_sites_opth$CM_PLH_obs-figure_sites_opth$CM_PLH_camelia_mem_opth_opth)/figure_sites_opth$CM_PLH_obs,
                                                                               abs(figure_sites_opth$CM_PLH_obs-figure_sites_opth$CM_PLH_camelia_obs_opth)/figure_sites_opth$CM_PLH_obs,
                                                                               abs(figure_sites_opth$CM_PLH_obs-figure_sites_opth$CM_PLH_camelia_proxy_opth)/figure_sites_opth$CM_PLH_obs),
                                        "SLA abs((CM obs - CM pred)/CM obs)"=c(abs(figure_sites_opth$CM_SLA_obs-figure_sites_opth$CM_SLA_sdm_opth)/figure_sites_opth$CM_SLA_obs,
                                                                               abs(figure_sites_opth$CM_SLA_obs-figure_sites_opth$CM_SLA_camelia_mem_opth_opth)/figure_sites_opth$CM_SLA_obs,
                                                                               abs(figure_sites_opth$CM_SLA_obs-figure_sites_opth$CM_SLA_camelia_obs_opth)/figure_sites_opth$CM_SLA_obs,
                                                                               abs(figure_sites_opth$CM_SLA_obs-figure_sites_opth$CM_SLA_camelia_proxy_opth)/figure_sites_opth$CM_SLA_obs),
                                        "LNC abs((CM obs - CM pred)/CM obs)"=c(abs(figure_sites_opth$CM_LNC_obs-figure_sites_opth$CM_LNC_sdm_opth)/figure_sites_opth$CM_LNC_obs,
                                                                               abs(figure_sites_opth$CM_LNC_obs-figure_sites_opth$CM_LNC_camelia_mem_opth_opth)/figure_sites_opth$CM_LNC_obs,
                                                                               abs(figure_sites_opth$CM_LNC_obs-figure_sites_opth$CM_LNC_camelia_obs_opth)/figure_sites_opth$CM_LNC_obs,
                                                                               abs(figure_sites_opth$CM_LNC_obs-figure_sites_opth$CM_LNC_camelia_proxy_opth)/figure_sites_opth$CM_LNC_obs),
                                        "Models"=c(rep("ssdm",4463),rep("camelia mem",4463),rep("camelia obs",4463),rep("camelia proxy",4463))
)
data_richness_precision_opth$Models <- factor(data_richness_precision_opth$Models, levels = c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY"))

### CM
plot2 <- ggplot(data=data_richness_precision_opth, aes(x=Models, y=data_richness_precision_opth$PLH.abs..CM.obs...CM.pred..CM.obs., fill=Models))+
  geom_boxplot()+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) + 
  theme_minimal() +
  labs(x = NULL) +
  theme(legend.position = "none", axis.title.y = element_text(size = 14))
plot1 <- ggplot(data=data_richness_precision_opth, aes(x=Models, y=data_richness_precision_opth$log.rich_pred.rich_obs., fill=Models))+
  geom_boxplot()+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) + 
  theme_minimal() +
  labs(x = NULL) +
  theme(legend.position = "none", axis.title.y = element_text(size = 14))
plot.x <- layer_data(plot1)
plot.y <- layer_data(plot2)
colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY")
df$category <- factor(df$category , levels = c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY"))
df.outliers <- df %>%
  select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
x.outliers <- data.frame("category"=c(rep("SSDM", length(df.outliers$x.outliers[[1]])),
                                      rep("CAMELIA MEM", length(df.outliers$x.outliers[[2]])),
                                      rep("CAMELIA OBS", length(df.outliers$x.outliers[[3]])),
                                      rep("CAMELIA PROXY", length(df.outliers$x.outliers[[4]]))),
                         "x.outliers"=unlist(unlist(df.outliers$x.outliers)),
                         "y.middle"=c(rep(df.outliers$y.middle[1], length(df.outliers$x.outliers[[1]])),
                                      rep(df.outliers$y.middle[2], length(df.outliers$x.outliers[[2]])),
                                      rep(df.outliers$y.middle[3], length(df.outliers$x.outliers[[3]])),
                                      rep(df.outliers$y.middle[4], length(df.outliers$x.outliers[[4]])))
)
y.outliers <- data.frame("category"=c(rep("SSDM", length(df.outliers$y.outliers[[1]])),
                                      rep("CAMELIA MEM", length(df.outliers$y.outliers[[2]])),
                                      rep("CAMELIA OBS", length(df.outliers$y.outliers[[3]])),
                                      rep("CAMELIA PROXY", length(df.outliers$y.outliers[[4]]))),
                         "y.outliers"=unlist(unlist(df.outliers$y.outliers)),
                         "x.middle"=c(rep(df.outliers$x.middle[1], length(df.outliers$y.outliers[[1]])),
                                      rep(df.outliers$x.middle[2], length(df.outliers$y.outliers[[2]])),
                                      rep(df.outliers$x.middle[3], length(df.outliers$y.outliers[[3]])),
                                      rep(df.outliers$x.middle[4], length(df.outliers$y.outliers[[4]]))))
boxplot_CM_PLH <- ggplot()+
  geom_rect(data = df, aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper, fill = category, color = category), alpha = 0.3) +
  geom_rect(data = df, aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper, fill = category, color = category), color = "black", fill = NA) +
  geom_segment(data = df, aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle, color = category)) + #whiskers
  geom_segment(data = df, aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper, color = category)) + #lower end
  geom_segment(data = df, aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper, color = category)) + #upper end
  geom_segment(data = df, aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max, color = category)) + #whiskers
  geom_segment(data = df, aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min, color = category)) + #lower end
  geom_segment(data = df, aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max, color = category)) + #upper end
  geom_point(data = x.outliers, aes(x = x.outliers, y = y.middle, fill = category, color = category), size = 3, shape = 1) + # x-direction
  geom_point(data = y.outliers, aes(x = x.middle, y = y.outliers, fill = category, color = category), size = 3, shape = 1) +  # y-direction
  geom_segment(data = df, aes(x=x.middle, y=y.middle, xend=x.middle, yend=0, color = category), linetype="dashed")+
  geom_segment(data = df, aes(x=x.middle, y=y.middle, xend=-Inf, yend=y.middle, color = category), linetype="dashed")+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) +
  scale_color_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) +
  scale_y_continuous(trans = log10_trans(), labels = scales::comma, limits = c(0.0000005, 10)) +
  coord_trans(y='reverse') +
  xlim(-0.75,4.25)+
  xlab(bquote('log(richness pred/richness obs)')) + 
  ylab(bquote('abs((CM obs - CM pred)/CM obs)')) + 
  theme_bw()+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none")

### CSTD
plot2 <- ggplot(data=data_richness_precision_opth, aes(x=Models, y=data_richness_precision_opth$PLH.abs..CSTD.obs...CSTD.pred..CSTD.obs., fill=Models))+
  geom_boxplot()+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) + 
  theme_minimal() +
  labs(x = NULL) +
  theme(legend.position = "none", axis.title.y = element_text(size = 14))
plot1 <- ggplot(data=data_richness_precision_opth, aes(x=Models, y=data_richness_precision_opth$log.rich_pred.rich_obs., fill=Models))+
  geom_boxplot()+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) + 
  theme_minimal() +
  labs(x = NULL) +
  theme(legend.position = "none", axis.title.y = element_text(size = 14))
plot.x <- layer_data(plot1)
plot.y <- layer_data(plot2)
colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY")
df$category <- factor(df$category , levels = c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY"))
df.outliers <- df %>%
  select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
x.outliers <- data.frame("category"=c(rep("SSDM", length(df.outliers$x.outliers[[1]])),
                                      rep("CAMELIA MEM", length(df.outliers$x.outliers[[2]])),
                                      rep("CAMELIA OBS", length(df.outliers$x.outliers[[3]])),
                                      rep("CAMELIA PROXY", length(df.outliers$x.outliers[[4]]))),
                         "x.outliers"=unlist(unlist(df.outliers$x.outliers)),
                         "y.middle"=c(rep(df.outliers$y.middle[1], length(df.outliers$x.outliers[[1]])),
                                      rep(df.outliers$y.middle[2], length(df.outliers$x.outliers[[2]])),
                                      rep(df.outliers$y.middle[3], length(df.outliers$x.outliers[[3]])),
                                      rep(df.outliers$y.middle[4], length(df.outliers$x.outliers[[4]])))
)
y.outliers <- data.frame("category"=c(rep("SSDM", length(df.outliers$y.outliers[[1]])),
                                      rep("CAMELIA MEM", length(df.outliers$y.outliers[[2]])),
                                      rep("CAMELIA OBS", length(df.outliers$y.outliers[[3]])),
                                      rep("CAMELIA PROXY", length(df.outliers$y.outliers[[4]]))),
                         "y.outliers"=unlist(unlist(df.outliers$y.outliers)),
                         "x.middle"=c(rep(df.outliers$x.middle[1], length(df.outliers$y.outliers[[1]])),
                                      rep(df.outliers$x.middle[2], length(df.outliers$y.outliers[[2]])),
                                      rep(df.outliers$x.middle[3], length(df.outliers$y.outliers[[3]])),
                                      rep(df.outliers$x.middle[4], length(df.outliers$y.outliers[[4]]))))
boxplot_CSTD_PLH <- ggplot()+
  geom_rect(data = df, aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper, fill = category, color = category), alpha = 0.3) +
  geom_rect(data = df, aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper, fill = category, color = category), color = "black", fill = NA) +
  geom_segment(data = df, aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle, color = category)) + #whiskers
  geom_segment(data = df, aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper, color = category)) + #lower end
  geom_segment(data = df, aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper, color = category)) + #upper end
  geom_segment(data = df, aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max, color = category)) + #whiskers
  geom_segment(data = df, aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min, color = category)) + #lower end
  geom_segment(data = df, aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max, color = category)) + #upper end
  geom_point(data = x.outliers, aes(x = x.outliers, y = y.middle, fill = category, color = category), size = 3, shape = 1) + # x-direction
  geom_point(data = y.outliers, aes(x = x.middle, y = y.outliers, fill = category, color = category), size = 3, shape = 1) +  # y-direction
  geom_segment(data = df, aes(x=x.middle, y=y.middle, xend=x.middle, yend=0, color = category), linetype="dashed")+
  geom_segment(data = df, aes(x=x.middle, y=y.middle, xend=-Inf, yend=y.middle, color = category), linetype="dashed")+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) +
  scale_color_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) +
  scale_y_continuous(trans=log10_trans(), labels = scales::comma, limits = c(0.0000005, 10)) +
  coord_trans(y='reverse') +
  xlim(-0.75,4.25)+
  xlab(bquote('log(richness pred/richness obs)')) + 
  ylab(bquote('abs((CSTD obs - CSTD pred)/CSTD obs)')) + 
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none")

## SORENSEN
plot2 <- ggplot(data=data_richness_precision_opth, aes(x=Models, y= Sorensen, fill=Models))+
  geom_boxplot()+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) + 
  theme_minimal() +
  labs(x = NULL) +
  ylim(0,1)+
  theme(legend.position = "none", axis.title.y = element_text(size = 14))
plot1 <- ggplot(data=data_richness_precision_opth, aes(x=Models, y=data_richness_precision_opth$log.rich_pred.rich_obs., fill=Models))+
  geom_boxplot()+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) + 
  theme_minimal() +
  labs(x = NULL) +
  theme(legend.position = "none", axis.title.y = element_text(size = 14))
plot.x <- layer_data(plot1)
plot.y <- layer_data(plot2)
colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY")
df$category <- factor(df$category , levels = c("SSDM","CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY"))
df.outliers <- df %>%
  select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
x.outliers <- data.frame("category"=c(rep("SSDM", length(df.outliers$x.outliers[[1]])),
                                      rep("CAMELIA MEM", length(df.outliers$x.outliers[[2]])),
                                      rep("CAMELIA OBS", length(df.outliers$x.outliers[[3]])),
                                      rep("CAMELIA PROXY", length(df.outliers$x.outliers[[4]]))),
                         "x.outliers"=unlist(unlist(df.outliers$x.outliers)),
                         "y.middle"=c(rep(df.outliers$y.middle[1], length(df.outliers$x.outliers[[1]])),
                                      rep(df.outliers$y.middle[2], length(df.outliers$x.outliers[[2]])),
                                      rep(df.outliers$y.middle[3], length(df.outliers$x.outliers[[3]])),
                                      rep(df.outliers$y.middle[4], length(df.outliers$x.outliers[[4]])))
)
y.outliers <- data.frame("category"=c(rep("SSDM", length(df.outliers$y.outliers[[1]])),
                                      rep("CAMELIA MEM", length(df.outliers$y.outliers[[2]])),
                                      rep("CAMELIA OBS", length(df.outliers$y.outliers[[3]])),
                                      rep("CAMELIA PROXY", length(df.outliers$y.outliers[[4]]))),
                         "y.outliers"=unlist(unlist(df.outliers$y.outliers)),
                         "x.middle"=c(rep(df.outliers$x.middle[1], length(df.outliers$y.outliers[[1]])),
                                      rep(df.outliers$x.middle[2], length(df.outliers$y.outliers[[2]])),
                                      rep(df.outliers$x.middle[3], length(df.outliers$y.outliers[[3]])),
                                      rep(df.outliers$x.middle[4], length(df.outliers$y.outliers[[4]]))))
boxplot_Sorensen <- ggplot()+
  geom_rect(data = df, aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper, fill = category, color = category), alpha = 0.3) +
  geom_rect(data = df, aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper, fill = category, color = category), color = "black", fill = NA) +
  geom_segment(data = df, aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle, color = category)) + #whiskers
  geom_segment(data = df, aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper, color = category)) + #lower end
  geom_segment(data = df, aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper, color = category)) + #upper end
  geom_segment(data = df, aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max, color = category)) + #whiskers
  geom_segment(data = df, aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min, color = category)) + #lower end
  geom_segment(data = df, aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max, color = category)) + #upper end
  geom_point(data = x.outliers, aes(x = x.outliers, y = y.middle, fill = category, color = category), size = 3, shape = 1) + # x-direction
  geom_point(data = y.outliers, aes(x = x.middle, y = y.outliers, fill = category, color = category), size = 3, shape = 1) +  # y-direction
  geom_segment(data = df, aes(x=x.middle, y=y.middle, xend=x.middle, yend=0, color = category), linetype="dashed")+
  geom_segment(data = df, aes(x=x.middle, y=y.middle, xend=-Inf, yend=y.middle, color = category), linetype="dashed")+
  scale_fill_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222")) +
  scale_color_manual(values = c( "#E1812C","darkblue", "darkgreen", "#B22222"))  +
  xlim(-0.75,4.25)+
  ylim(-0.01, 1)+
  xlab(bquote('log(richness pred/richness obs)')) + 
  ylab(bquote('Sorensen')) + 
  theme_bw()+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none")


### FULL PLOT
plot_full <- ggarrange( 
  Sorensen_logrichness, 
  boxplot_Sorensen,
  plot_diff_CM_PLH, 
  boxplot_CM_PLH,
  plot_diff_CSTD_PLH,
  boxplot_CSTD_PLH, 
  labels = c("A", "B", "C", "D","E","F"),
  ncol = 2, nrow = 3, align = "hv",common.legend = TRUE, legend = "bottom", widths = c(3, 2))

ggsave("figure_2.png", plot = plot_full, width = 8.27, height = 11.69, units = "in")
ggsave("figure_2.pdf", plot = plot_full, width = 8.27, height = 11.69, units = "in", dpi = 300)

#### -- Supp Mat -- ####
# CORRELATIONS community indices:
table_S2 <- data.frame("Indices"=c(rep('CM PLH',5),
                                   rep('CM SLA',5),
                                   rep('CM LNC',5),
                                   rep('CSTD PLH',5),
                                   rep('CSTD SLA',5),
                                   rep('CSTD LNC',5),
                                   rep('SR',5)),
                       "Models"=rep(c('MEM', 'SSDM', 'CAMELIA MEM', 'CAMELIA OBS', 'CAMELIA PROXY'), 7), 
                       "Probability R2"=c(cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_mem)^2,
                                          cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_sdm_raw)^2,
                                          cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_mem_raw)^2,
                                          cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_obs_raw)^2,
                                          cor(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_proxy_raw)^2,
                                          cor(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_mem)^2,
                                          cor(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_sdm_raw)^2,
                                          cor(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_camelia_mem_raw)^2,
                                          cor(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_camelia_obs_raw)^2,
                                          cor(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_camelia_proxy_raw)^2,
                                          cor(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_mem)^2,
                                          cor(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_sdm_raw)^2,
                                          cor(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_camelia_mem_raw)^2,
                                          cor(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_camelia_obs_raw)^2,
                                          cor(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_camelia_proxy_raw)^2,
                                          
                                          cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_mem)^2,
                                          cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_sdm_raw)^2,
                                          cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_mem_raw)^2,
                                          cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_obs_raw)^2,
                                          cor(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_proxy_raw)^2,
                                          cor(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_mem)^2,
                                          cor(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_sdm_raw)^2,
                                          cor(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_camelia_mem_raw)^2,
                                          cor(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_camelia_obs_raw)^2,
                                          cor(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_camelia_proxy_raw)^2,
                                          cor(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_mem)^2,
                                          cor(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_sdm_raw)^2,
                                          cor(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_camelia_mem_raw)^2,
                                          cor(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_camelia_obs_raw)^2,
                                          cor(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_camelia_proxy_raw)^2,
                                          
                                          cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_mem)^2,
                                          cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_sdm_raw)^2,
                                          cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_mem_raw)^2,
                                          cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_obs_raw)^2,
                                          cor(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_proxy_raw)^2), 
                       "Probability RMSE" = c(rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_mem),
                                              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_sdm_raw),
                                              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_mem_raw),
                                              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_obs_raw),
                                              rmse(CM_tot_raw$PLH_obs, CM_tot_raw$PLH_camelia_proxy_raw),
                                              rmse(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_mem),
                                              rmse(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_sdm_raw),
                                              rmse(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_camelia_mem_raw),
                                              rmse(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_camelia_obs_raw),
                                              rmse(CM_tot_raw$SLA_obs, CM_tot_raw$SLA_camelia_proxy_raw),
                                              rmse(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_mem),
                                              rmse(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_sdm_raw),
                                              rmse(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_camelia_mem_raw),
                                              rmse(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_camelia_obs_raw),
                                              rmse(CM_tot_raw$LNC_obs, CM_tot_raw$LNC_camelia_proxy_raw),
                                              
                                              rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_mem),
                                              rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_sdm_raw),
                                              rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_mem_raw),
                                              rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_obs_raw),
                                              rmse(CSTD_tot_raw$PLH_obs, CSTD_tot_raw$PLH_camelia_proxy_raw),
                                              rmse(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_mem),
                                              rmse(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_sdm_raw),
                                              rmse(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_camelia_mem_raw),
                                              rmse(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_camelia_obs_raw),
                                              rmse(CSTD_tot_raw$SLA_obs, CSTD_tot_raw$SLA_camelia_proxy_raw),
                                              rmse(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_mem),
                                              rmse(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_sdm_raw),
                                              rmse(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_camelia_mem_raw),
                                              rmse(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_camelia_obs_raw),
                                              rmse(CSTD_tot_raw$LNC_obs, CSTD_tot_raw$LNC_camelia_proxy_raw),
                                              
                                              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_mem),
                                              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_sdm_raw),
                                              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_mem_raw),
                                              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_obs_raw),
                                              rmse(SR_tot_raw$SR_obs, SR_tot_raw$SR_camelia_proxy_raw)), 
                       "Binarised by TSS-threshold R2"=c(cor(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_mem)^2,
                                                         cor(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_sdm_opth)^2,
                                                         cor(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_camelia_mem_opth)^2,
                                                         cor(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_camelia_obs_opth)^2,
                                                         cor(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_camelia_proxy_opth)^2,
                                                         cor(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_mem)^2,
                                                         cor(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_sdm_opth)^2,
                                                         cor(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_camelia_mem_opth)^2,
                                                         cor(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_camelia_obs_opth)^2,
                                                         cor(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_camelia_proxy_opth)^2,
                                                         cor(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_mem)^2,
                                                         cor(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_sdm_opth)^2,
                                                         cor(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_camelia_mem_opth)^2,
                                                         cor(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_camelia_obs_opth)^2,
                                                         cor(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_camelia_proxy_opth)^2,
                                                         
                                                         cor(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_mem)^2,
                                                         cor(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_sdm_opth)^2,
                                                         cor(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_camelia_mem_opth)^2,
                                                         cor(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_camelia_obs_opth)^2,
                                                         cor(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_camelia_proxy_opth)^2,
                                                         cor(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_mem)^2,
                                                         cor(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_sdm_opth)^2,
                                                         cor(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_camelia_mem_opth)^2,
                                                         cor(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_camelia_obs_opth)^2,
                                                         cor(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_camelia_proxy_opth)^2,
                                                         cor(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_mem)^2,
                                                         cor(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_sdm_opth)^2,
                                                         cor(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_camelia_mem_opth)^2,
                                                         cor(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_camelia_obs_opth)^2,
                                                         cor(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_camelia_proxy_opth)^2,
                                                         
                                                         cor(SR_tot_opth$SR_obs, SR_tot_opth$SR_mem)^2,
                                                         cor(SR_tot_opth$SR_obs, SR_tot_opth$SR_sdm_opth)^2,
                                                         cor(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_mem_opth)^2,
                                                         cor(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_obs_opth)^2,
                                                         cor(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_proxy_opth)^2), 
                       "Binarised by TSS-threshold RMSE"=c(rmse(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_mem),
                                                           rmse(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_sdm_opth),
                                                           rmse(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_camelia_mem_opth),
                                                           rmse(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_camelia_obs_opth),
                                                           rmse(CM_tot_opth$PLH_obs, CM_tot_opth$PLH_camelia_proxy_opth),
                                                           rmse(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_mem),
                                                           rmse(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_sdm_opth),
                                                           rmse(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_camelia_mem_opth),
                                                           rmse(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_camelia_obs_opth),
                                                           rmse(CM_tot_opth$SLA_obs, CM_tot_opth$SLA_camelia_proxy_opth),
                                                           rmse(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_mem),
                                                           rmse(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_sdm_opth),
                                                           rmse(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_camelia_mem_opth),
                                                           rmse(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_camelia_obs_opth),
                                                           rmse(CM_tot_opth$LNC_obs, CM_tot_opth$LNC_camelia_proxy_opth),
                                                           
                                                           rmse(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_mem),
                                                           rmse(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_sdm_opth),
                                                           rmse(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_camelia_mem_opth),
                                                           rmse(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_camelia_obs_opth),
                                                           rmse(CSTD_tot_opth$PLH_obs, CSTD_tot_opth$PLH_camelia_proxy_opth),
                                                           rmse(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_mem),
                                                           rmse(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_sdm_opth),
                                                           rmse(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_camelia_mem_opth),
                                                           rmse(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_camelia_obs_opth),
                                                           rmse(CSTD_tot_opth$SLA_obs, CSTD_tot_opth$SLA_camelia_proxy_opth),
                                                           rmse(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_mem),
                                                           rmse(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_sdm_opth),
                                                           rmse(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_camelia_mem_opth),
                                                           rmse(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_camelia_obs_opth),
                                                           rmse(CSTD_tot_opth$LNC_obs, CSTD_tot_opth$LNC_camelia_proxy_opth),
                                                           
                                                           rmse(SR_tot_opth$SR_obs, SR_tot_opth$SR_mem),
                                                           rmse(SR_tot_opth$SR_obs, SR_tot_opth$SR_sdm_opth),
                                                           rmse(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_mem_opth),
                                                           rmse(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_obs_opth),
                                                           rmse(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_proxy_opth)))

# Species predictions evaluation : 
## RAW
AUC_tot <- data.frame("Species"=rep(AUC_tot$Species,3), "Quantiles"=rep(AUC_tot$Quantiles,3), "AUC"=rep(0,831*3), "pr-AUC"=rep(0,831*3), "Approach"=c(rep("SSDM",831), rep("CAMELIA MEM",831), rep("CAMELIA OBS",831), rep("CAMELIA PROXY", 831)))
for (i in 1:831){
  AUC_raw$AUC[i] = auc(survey_obs[,i], as.numeric(unlist(sdm_raw[,i])))
  pr_curve <- pr.curve(survey_obs[,i], as.numeric(unlist(sdm_raw[,i])))
  AUC_raw$pr.AUC[i]= pr_curve$auc.integral
  AUC_raw$AUC[831+i]=auc(survey_obs[,i], camelia_mem_raw[,i])
  pr_curve <- pr.curve(survey_obs[,i], camelia_mem_raw[,i])
  AUC_raw$pr.AUC[831+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*2+i]=auc(survey_obs[,i], camelia_obs_raw[,i])
  pr_curve <- pr.curve(survey_obs[,i], camelia_obs_raw[,i])
  AUC_raw$pr.AUC[831*2+i]= pr_curve$auc.integral
  AUC_raw$AUC[831*3+i]=auc(survey_obs[,i], camelia_proxy_raw[,i])
  pr_curve <- pr.curve(survey_obs[,i], camelia_proxy_raw[,i])
  AUC_raw$pr.AUC[831*3+i]= pr_curve$auc.integral
  print(i)
}

## BINARIZED USING PER SPECIES TSS-OPTIMIZED THRESHOLD
eval_species <- data.frame("species"=colnames(survey_obs), 
                           "prevalence"=colSums(survey_obs),
                           "Specificite_SSDM"=rep(0,831),
                           "Sensibilite_SSDM"=rep(0,831),
                           "Precision_SSDM"=rep(0,831),
                           "TSS_SSDM"=rep(0,831),
                           "Specificite_CAMELIA_MEM"=rep(0,831),
                           "Sensibilite_CAMELIA_MEM"=rep(0,831),
                           "Precision_CAMELIA_MEM"=rep(0,831),
                           "TSS_CAMELIA_MEM"=rep(9,831),
                           "Specificite_CAMELIA_OBS"=rep(0,831),
                           "Sensibilite_CAMELIA_OBS"=rep(0,831),
                           "Precision_CAMELIA_OBS"=rep(0,831),
                           "TSS_CAMELIA_OBS"=rep(0,831),
                           "Specificite_CAMELIA_PROXY"=rep(0,831),
                           "Sensibilite_CAMELIA_PROXY"=rep(0,831),
                           "Precision_CAMELIA_PROXY"=rep(0,831),
                           "TSS_CAMELIA_PROXY"=rep(0,831))

eval_species_opth <- eval_species
for (i in 1:831){
  VP_SSDM=length(which(which(sdm_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_SSDM=length(which(which(sdm_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_SSDM=length(which(which(sdm_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_SSDM=length(which(which(sdm_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_SSDM[i]=VP_SSDM/(VP_SSDM+FN_SSDM)
  eval_species_opth$Specificite_SSDM[i]=VN_SSDM/(VN_SSDM+FP_SSDM)
  eval_species_opth$Precision_SSDM[i]=VP_SSDM/(VP_SSDM+FP_SSDM)
  
  VP_CAMELIA_MEM=length(which(which(camelia_mem_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_CAMELIA_MEM=length(which(which(camelia_mem_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_CAMELIA_MEM=length(which(which(camelia_mem_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_CAMELIA_MEM=length(which(which(camelia_mem_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_MEM[i]=VP_CAMELIA_MEM/(VP_CAMELIA_MEM+FN_CAMELIA_MEM)
  eval_species_opth$Specificite_CAMELIA_MEM[i]=VN_CAMELIA_MEM/(VN_CAMELIA_MEM+FP_CAMELIA_MEM)
  eval_species_opth$Precision_CAMELIA_MEM[i]=VP_CAMELIA_MEM/(VP_CAMELIA_MEM+FP_CAMELIA_MEM)
  
  VP_CAMELIA_OBS=length(which(which(camelia_obs_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_CAMELIA_OBS=length(which(which(camelia_obs_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_CAMELIA_OBS=length(which(which(camelia_obs_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_CAMELIA_OBS=length(which(which(camelia_obs_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_OBS[i]=VP_CAMELIA_OBS/(VP_CAMELIA_OBS+FN_CAMELIA_OBS)
  eval_species_opth$Specificite_CAMELIA_OBS[i]=VN_CAMELIA_OBS/(VN_CAMELIA_OBS+FP_CAMELIA_OBS)
  eval_species_opth$Precision_CAMELIA_OBS[i]=VP_CAMELIA_OBS/(VP_CAMELIA_OBS+FP_CAMELIA_OBS)
  
  VP_CAMELIA_PROXY=length(which(which(camelia_proxy_opth[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_CAMELIA_PROXY=length(which(which(camelia_proxy_opth[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_CAMELIA_PROXY=length(which(which(camelia_proxy_opth[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_CAMELIA_PROXY=length(which(which(camelia_proxy_opth[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_opth$Sensibilite_CAMELIA_PROXY[i]=VP_CAMELIA_PROXY/(VP_CAMELIA_PROXY+FN_CAMELIA_PROXY)
  eval_species_opth$Specificite_CAMELIA_PROXY[i]=VN_CAMELIA_PROXY/(VN_CAMELIA_PROXY+FP_CAMELIA_PROXY)
  eval_species_opth$Precision_CAMELIA_PROXY[i]=VP_CAMELIA_PROXY/(VP_CAMELIA_PROXY+FP_CAMELIA_PROXY)
  print(i)
}
eval_species_opth$TSS_SSDM = eval_species_opth$Specificite_SSDM+eval_species_opth$Sensibilite_SSDM-1
eval_species_opth$TSS_CAMELIA_MEM = eval_species_opth$Specificite_CAMELIA_MEM+eval_species_opth$Sensibilite_CAMELIA_MEM-1
eval_species_opth$TSS_CAMELIA_OBS = eval_species_opth$Specificite_CAMELIA_OBS+eval_species_opth$Sensibilite_CAMELIA_OBS-1
eval_species_opth$TSS_CAMELIA_PROXY = eval_species_opth$Specificite_CAMELIA_PROXY+eval_species_opth$Sensibilite_CAMELIA_PROXY-1

## BINARIZED USING PRR METHOD
eval_species_prr <- eval_species
for (i in 1:831){
  VP_SSDM=length(which(which(sdm_prr[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_SSDM=length(which(which(sdm_prr[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_SSDM=length(which(which(sdm_prr[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_SSDM=length(which(which(sdm_prr[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_prr$Sensibilite_SSDM[i]=VP_SSDM/(VP_SSDM+FN_SSDM)
  eval_species_prr$Specificite_SSDM[i]=VN_SSDM/(VN_SSDM+FP_SSDM)
  eval_species_prr$Precision_SSDM[i]=VP_SSDM/(VP_SSDM+FP_SSDM)
  
  VP_CAMELIA_MEM=length(which(which(camelia_mem_prr[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_CAMELIA_MEM=length(which(which(camelia_mem_prr[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_CAMELIA_MEM=length(which(which(camelia_mem_prr[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_CAMELIA_MEM=length(which(which(camelia_mem_prr[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_prr$Sensibilite_CAMELIA_MEM[i]=VP_CAMELIA_MEM/(VP_CAMELIA_MEM+FN_CAMELIA_MEM)
  eval_species_prr$Specificite_CAMELIA_MEM[i]=VN_CAMELIA_MEM/(VN_CAMELIA_MEM+FP_CAMELIA_MEM)
  eval_species_prr$Precision_CAMELIA_MEM[i]=VP_CAMELIA_MEM/(VP_CAMELIA_MEM+FP_CAMELIA_MEM)
  
  VP_CAMELIA_OBS=length(which(which(camelia_obs_prr[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_CAMELIA_OBS=length(which(which(camelia_obs_prr[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_CAMELIA_OBS=length(which(which(camelia_obs_prr[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_CAMELIA_OBS=length(which(which(camelia_obs_prr[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_prr$Sensibilite_CAMELIA_OBS[i]=VP_CAMELIA_OBS/(VP_CAMELIA_OBS+FN_CAMELIA_OBS)
  eval_species_prr$Specificite_CAMELIA_OBS[i]=VN_CAMELIA_OBS/(VN_CAMELIA_OBS+FP_CAMELIA_OBS)
  eval_species_prr$Precision_CAMELIA_OBS[i]=VP_CAMELIA_OBS/(VP_CAMELIA_OBS+FP_CAMELIA_OBS)
  
  VP_CAMELIA_PROXY=length(which(which(camelia_proxy_prr[,i]==1) %in% which(survey_obs[,i]==1)))
  FP_CAMELIA_PROXY=length(which(which(camelia_proxy_prr[,i]==1) %in% which(survey_obs[,i]==0)))
  VN_CAMELIA_PROXY=length(which(which(camelia_proxy_prr[,i]==0) %in% which(survey_obs[,i]==0)))
  FN_CAMELIA_PROXY=length(which(which(camelia_proxy_prr[,i]==0) %in% which(survey_obs[,i]==1)))
  eval_species_prr$Sensibilite_CAMELIA_PROXY[i]=VP_CAMELIA_PROXY/(VP_CAMELIA_PROXY+FN_CAMELIA_PROXY)
  eval_species_prr$Specificite_CAMELIA_PROXY[i]=VN_CAMELIA_PROXY/(VN_CAMELIA_PROXY+FP_CAMELIA_PROXY)
  eval_species_prr$Precision_CAMELIA_PROXY[i]=VP_CAMELIA_PROXY/(VP_CAMELIA_PROXY+FP_CAMELIA_PROXY)
  print(i)
}
eval_species_prr$TSS_SSDM = eval_species_prr$Specificite_SSDM+eval_species_prr$Sensibilite_SSDM-1
eval_species_prr$TSS_CAMELIA_MEM = eval_species_prr$Specificite_CAMELIA_MEM+eval_species_prr$Sensibilite_CAMELIA_MEM-1
eval_species_prr$TSS_CAMELIA_OBS = eval_species_prr$Specificite_CAMELIA_OBS+eval_species_prr$Sensibilite_CAMELIA_OBS-1
eval_species_prr$TSS_CAMELIA_PROXY = eval_species_prr$Specificite_CAMELIA_PROXY+eval_species_prr$Sensibilite_CAMELIA_PROXY-1

