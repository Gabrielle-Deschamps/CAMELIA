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

## PROXY
### PROxy construction :
#CM_proxy <- data.frame("PLH_proxy"= CM_obs$PLH+runif(4463, min = -30*mean(CM_obs$PLH)/100, max =30*mean(CM_obs$PLH)/100),
#                       "SLA_proxy"=CM_obs$SLA+runif(4463, min = -30*mean(CM_obs$SLA)/100, max =30*mean(CM_obs$SLA)/100),
#                       "LNC_proxy"=CM_obs$LNC+runif(4463, min = -30*mean(CM_obs$LNC)/100, max =30*mean(CM_obs$LNC)/100),
#                       "ID"=CM_obs$ID)
#
#CSTD_proxy <- data.frame("PLH_proxy"=CSTD_obs$PLH+runif(4463, min = -30*mean(CSTD_obs$PLH)/100, max =30*mean(CSTD_obs$PLH)/100),
#                       "SLA_proxy"=CSTD_obs$SLA+runif(4463, min = -30*mean(CSTD_obs$SLA)/100, max =30*mean(CSTD_obs$SLA)/100),
#                       "LNC_proxy"=CSTD_obs$LNC+runif(4463, min = -30*mean(CSTD_obs$LNC)/100, max =30*mean(CSTD_obs$LNC)/100),
#                       "ID"=CSTD_obs$ID)
#
#SR_proxy <- data.frame("SR_proxy"=SR_obs$SR + runif(4463, min=-30*mean(SR_obs$SR)/100,max=30*mean(SR_obs$SR)/100),
#                       "ID"=CSTD_obs$ID)
#
#write.csv(CM_proxy, file="data/community_indices/CM_proxy.csv")
#write.csv(CSTD_proxy, file="data/community_indices/CSTD_proxy.csv")
#write.csv(SR_proxy, file="data/community_indices/SR_proxy.csv")

CM_proxy <- read.csv("data/community_indices/CM_proxy.csv", header = TRUE, row.names = 1)
colnames(CM_proxy) <- paste0(colnames(CM_proxy), "_proxy")
CM_proxy$ID <- rownames(CM_proxy)

CSTD_proxy <- read.csv("data/community_indices/CSTD_proxy.csv", header = TRUE, row.names = 1)
colnames(CSTD_proxy) <- paste0(colnames(CSTD_proxy), "_proxy")
CSTD_proxy$ID <- rownames(CSTD_proxy)

SR_proxy <- read.csv("data/community_indices/SR_proxy.csv", header = TRUE, row.names = 1)
colnames(SR_proxy) <- paste0(colnames(SR_proxy), "_proxy")
SR_proxy$ID <- rownames(SR_proxy)

## SSDM
sdm_raw_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_raw_", i, ".csv"), header = TRUE, row.names = 1))
sdm_raw <- do.call(rbind, sdm_raw_list)
sdm_raw <- sdm_raw[rownames(survey_obs), ]

sdm_opth_list <- lapply(1:4, function(i) read.csv(paste0("results/sdm_bin_opth_", i, ".csv"), header = TRUE, row.names = 1))
sdm_opth <- do.call(rbind, sdm_opth_list)
sdm_opth <- sdm_opth[rownames(survey_obs), ]

sdm_prr <- prr(sdm_raw, SR_mem)

## CAMELIA PREDICTIONS
### camelia obs
camelia_obs_raw <- load_camelia_predictions("OBS", "raw")
camelia_obs_opth <- load_camelia_predictions("OBS", "opth")
camelia_obs_prr <- prr(camelia_obs_raw, SR_mem)

### camelia mem
camelia_mem_raw <- load_camelia_predictions("mem", "raw")
camelia_mem_opth <- load_camelia_predictions("mem", "opth")
camelia_mem_prr <- prr(camelia_mem_raw, SR_mem)

### camelia proxy
camelia_proxy_raw <- load_camelia_predictions("proxy", "raw")
camelia_proxy_opth <- load_camelia_predictions("proxy", "opth")
camelia_proxy_prr <- prr(camelia_proxy_raw, SR_mem)

# COMMUNITY INDICES
# Compute indices for all datasets
indices_obs <- compute_indices(survey_obs, "obs")
indices_sdm_raw <- compute_indices(sdm_raw, "sdm_raw")
indices_sdm_opth <- compute_indices(sdm_opth, "sdm_opth")
indices_sdm_prr <- compute_indices(sdm_prr, "sdm_prr")
indices_camelia_obs_raw <- compute_indices(camelia_obs_raw, "camelia_obs_raw")
indices_camelia_obs_opth <- compute_indices(camelia_obs_opth, "camelia_obs_opth")
indices_camelia_obs_prr <- compute_indices(camelia_obs_prr, "camelia_obs_prr")
indices_camelia_mem_raw <- compute_indices(camelia_mem_raw, "camelia_mem_raw")
indices_camelia_mem_opth <- compute_indices(camelia_mem_opth, "camelia_mem_opth")
indices_camelia_mem_prr <- compute_indices(camelia_mem_prr, "camelia_mem_prr")
indices_camelia_proxy_raw <- compute_indices(camelia_proxy_raw, "camelia_proxy_raw")
indices_camelia_proxy_opth <- compute_indices(camelia_proxy_opth, "camelia_proxy_opth")
indices_camelia_proxy_prr <- compute_indices(camelia_proxy_prr, "camelia_proxy_prr")

# Merge CM indices for all methods
idx <- "CM"
CM_tot_raw <- left_join(indices_obs[[idx]], indices_sdm_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_raw[[idx]], by = "ID") %>%
  left_join(CM_mem, by = "ID") %>%
  left_join(CM_proxy, by = "ID")

CM_tot_opth <- left_join(indices_obs[[idx]], indices_sdm_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_opth[[idx]], by = "ID") %>%
  left_join(CM_mem, by = "ID") %>%
  left_join(CM_proxy, by = "ID")

CM_tot_prr <- left_join(indices_obs[[idx]], indices_sdm_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_prr[[idx]], by = "ID") %>%
  left_join(CM_mem, by = "ID") %>%
  left_join(CM_proxy, by = "ID")

# Merge CSTD indices for all methods
idx <- "CSTD"
CSTD_tot_raw <- left_join(indices_obs[[idx]], indices_sdm_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_raw[[idx]], by = "ID") %>%
  left_join(CSTD_mem, by = "ID") %>%
  left_join(CSTD_proxy, by = "ID")

CSTD_tot_opth <- left_join(indices_obs[[idx]], indices_sdm_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_opth[[idx]], by = "ID") %>%
  left_join(CSTD_mem, by = "ID") %>%
  left_join(CSTD_proxy, by = "ID")

CSTD_tot_prr <- left_join(indices_obs[[idx]], indices_sdm_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_prr[[idx]], by = "ID") %>%
  left_join(CSTD_mem, by = "ID") %>%
  left_join(CSTD_proxy, by = "ID")

# Merge SR indices for all methods
idx <- "SR"
SR_tot_raw <- left_join(indices_obs[[idx]], indices_sdm_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_raw[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_raw[[idx]], by = "ID") %>%
  left_join(SR_mem, by = "ID") %>%
  left_join(SR_proxy, by = "ID")

SR_tot_opth <- left_join(indices_obs[[idx]], indices_sdm_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_opth[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_opth[[idx]], by = "ID") %>%
  left_join(SR_mem, by = "ID") %>%
  left_join(SR_proxy, by = "ID")

SR_tot_prr <- left_join(indices_obs[[idx]], indices_sdm_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_mem_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_obs_prr[[idx]], by = "ID") %>%
  left_join(indices_camelia_proxy_prr[[idx]], by = "ID") %>%
  left_join(SR_mem, by = "ID") %>%
  left_join(SR_proxy, by = "ID")

#### -- FIGURE 3 -- ####
# Create a dataframe containing observed and predicted values for different models
figure_raw <- data.frame("obs_CM"=rep(CM_tot_raw$PLH_obs,4),
                         "Pred_CM"=c(CM_tot_raw$PLH_sdm_raw,
                                     CM_tot_raw$PLH_camelia_obs_raw,
                                     CM_tot_raw$PLH_camelia_mem_raw,
                                     CM_tot_raw$PLH_camelia_proxy_raw),
                         'obs_CSTD'=rep(CSTD_tot_raw$PLH_obs,4),
                         "Pred_CSTD"=c(CSTD_tot_raw$PLH_sdm_raw,
                                       CSTD_tot_raw$PLH_camelia_obs_raw,
                                       CSTD_tot_raw$PLH_camelia_mem_raw,
                                       CSTD_tot_raw$PLH_camelia_proxy_raw),
                         'obs_SR' =rep(SR_tot_raw$SR_obs,4),
                         'Pred_SR' =c(SR_tot_raw$SR_sdm_raw,
                                      SR_tot_raw$SR_camelia_obs_raw,
                                      SR_tot_raw$SR_camelia_mem_raw,
                                      SR_tot_raw$SR_camelia_proxy_raw),
                         "Models"=c(rep("SSDM", 4463), rep("CAMELIA OBS", 4463), rep("CAMELIA MEM", 4463), rep("CAMELIA PROXY", 4463))
)
figure_raw$Models <- factor(figure_raw$Models, levels = c("SSDM", "CAMELIA MEM", "CAMELIA PROXY", "CAMELIA OBS"))

# Compute R² and RMSE for each model
cor_metrics <- data.frame(
  Models = factor(c("SSDM", "CAMELIA MEM", "CAMELIA PROXY", "CAMELIA OBS"), levels = c("SSDM", "CAMELIA MEM", "CAMELIA PROXY", "CAMELIA OBS")),
  R2_CM = sapply(list(CM_tot_raw$PLH_sdm_raw, CM_tot_raw$PLH_camelia_mem_raw, CM_tot_raw$PLH_camelia_proxy_raw, CM_tot_raw$PLH_camelia_obs_raw), function(x) cor(CM_tot_raw$PLH_obs, x)^2),
  RMSE_CM = sapply(list(CM_tot_raw$PLH_sdm_raw, CM_tot_raw$PLH_camelia_mem_raw, CM_tot_raw$PLH_camelia_proxy_raw, CM_tot_raw$PLH_camelia_obs_raw), function(x) rmse(CM_tot_raw$PLH_obs, x)),
  R2_CSTD = sapply(list(CSTD_tot_raw$PLH_sdm_raw, CSTD_tot_raw$PLH_camelia_mem_raw, CSTD_tot_raw$PLH_camelia_proxy_raw, CSTD_tot_raw$PLH_camelia_obs_raw), function(x) cor(CSTD_tot_raw$PLH_obs, x)^2),
  RMSE_CSTD = sapply(list(CSTD_tot_raw$PLH_sdm_raw, CSTD_tot_raw$PLH_camelia_mem_raw, CSTD_tot_raw$PLH_camelia_proxy_raw, CSTD_tot_raw$PLH_camelia_obs_raw), function(x) rmse(CSTD_tot_raw$PLH_obs, x)),
  R2_SR = sapply(list(SR_tot_raw$SR_sdm_raw, SR_tot_raw$SR_camelia_mem_raw, SR_tot_raw$SR_camelia_proxy_raw, SR_tot_raw$SR_camelia_obs_raw), function(x) cor(SR_tot_raw$SR_obs, x)^2),
  RMSE_SR = sapply(list(SR_tot_raw$SR_sdm_raw, SR_tot_raw$SR_camelia_mem_raw, SR_tot_raw$SR_camelia_proxy_raw, SR_tot_raw$SR_camelia_obs_raw), function(x) rmse(SR_tot_raw$SR_obs, x))
)

plot_fig_1 <- function(obs, pred, lab_text, cor_column, rmse_column) {
  ggplot(figure_raw, aes(y = obs, x = pred, color = Models, fill = Models)) +
    geom_hdr() +
    geom_point(shape = 21, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("#E1812C", "darkblue", "#B22222", "darkgreen")) +
    scale_fill_manual(values = c("#E1812C", "darkblue", "#B22222", "darkgreen")) +
    facet_wrap(vars(Models), nrow = 1) +
    ylab(paste0("Observed ",lab_text)) +
    xlab(paste0("Predicted ",lab_text)) +
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
                   label = sprintf("R² = %.2f\nRMSE = %.2f", round(get(cor_column), 2), round(get(rmse_column), 2))),
               color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)
}

# Generate plots
CM_raw <- plot_fig_1(figure_raw$obs_CM, figure_raw$Pred_CM, "CM PLH", "R2_CM", "RMSE_CM")
CSTD_raw <- plot_fig_1(figure_raw$obs_CSTD, figure_raw$Pred_CSTD, "CSTD PLH", "R2_CSTD", "RMSE_CSTD")
SR_raw <- plot_fig_1(figure_raw$obs_SR, figure_raw$Pred_SR, "SR", "R2_SR", "RMSE_SR")

# Combine plots
combined_plot <- SR_raw / CM_raw / CSTD_raw +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("figure_1.pdf", 
       plot = combined_plot, 
       device = "pdf", 
       width = 8.27,
       height = 11.69,
       dpi = 300)

## -- FIGURE S1 -- ##

# Create a dataframe for MEM and PROXY predictions
figure_comm <- data.frame(
  "obs_CM" = rep(CM_tot_raw$PLH_obs, 2),
  "Pred_CM" = c(CM_tot_raw$PLH_mem, CM_tot_raw$PLH_proxy),
  "obs_CSTD" = rep(CSTD_tot_raw$PLH_obs, 2),
  "Pred_CSTD" = c(CSTD_tot_raw$PLH_mem, CSTD_tot_raw$PLH_proxy),
  "obs_SR" = rep(SR_tot_raw$SR_obs, 2),
  "Pred_SR" = c(SR_tot_raw$SR_mem, SR_tot_raw$SR_proxy),
  "Models" = rep(c("MEM", "PROXY"), each = nrow(CM_tot_raw))
)
figure_comm$Models <- factor(figure_comm$Models, levels = c("MEM", "PROXY"))

# Compute R² and RMSE for MEM and PROXY
cor_metrics_comm <- data.frame(
  Models = factor(c("MEM", "PROXY"), levels = c("MEM", "PROXY")),
  R2_CM = sapply(list(CM_tot_raw$PLH_mem, CM_tot_raw$PLH_proxy), function(x) cor(CM_tot_raw$PLH_obs, x)^2),
  RMSE_CM = sapply(list(CM_tot_raw$PLH_mem, CM_tot_raw$PLH_proxy), function(x) rmse(CM_tot_raw$PLH_obs, x)),
  R2_CSTD = sapply(list(CSTD_tot_raw$PLH_mem, CSTD_tot_raw$PLH_proxy), function(x) cor(CSTD_tot_raw$PLH_obs, x)^2),
  RMSE_CSTD = sapply(list(CSTD_tot_raw$PLH_mem, CSTD_tot_raw$PLH_proxy), function(x) rmse(CSTD_tot_raw$PLH_obs, x)),
  R2_SR = sapply(list(SR_tot_raw$SR_mem, SR_tot_raw$SR_proxy), function(x) cor(SR_tot_raw$SR_obs, x)^2),
  RMSE_SR = sapply(list(SR_tot_raw$SR_mem, SR_tot_raw$SR_proxy), function(x) rmse(SR_tot_raw$SR_obs, x))
)

plot_fig_S1 <- function(obs, pred, lab_text, cor_column, rmse_column) {
  ggplot(figure_comm, aes(y = obs, x = pred, color = Models, fill = Models)) +
    geom_hdr() +
    geom_point(shape = 21, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  
    scale_color_manual(values = c("#3274A1", "#FFC0CB")) +
    scale_fill_manual(values = c("#3274A1", "#FFC0CB")) +
    facet_wrap(vars(Models), nrow = 1) +
    ylab(paste0("Observed ",lab_text)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 10)
    ) +
    geom_label(data = cor_metrics_comm,
               aes(x = -Inf, y = Inf,
                   label = sprintf("R² = %.2f\nRMSE = %.2f", round(get(cor_column), 2), round(get(rmse_column), 2))),
               color = "black", fill = "white", hjust = -0.1, vjust = 1.1, size = 4, inherit.aes = FALSE)
}
# Generate plots for Figure S1
CM_comm <- plot_fig_S1(figure_comm$obs_CM, figure_comm$Pred_CM, "CM PLH", "R2_CM", "RMSE_CM")
CSTD_comm <- plot_fig_S1(figure_comm$obs_CSTD, figure_comm$Pred_CSTD, "CSTD PLH", "R2_CSTD", "RMSE_CSTD")
SR_comm <- plot_fig_S1(figure_comm$obs_SR, figure_comm$Pred_SR, "SR", "R2_SR", "RMSE_SR") 

# Combine plots
comm_plot <- SR_comm / CM_comm / CSTD_comm +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

# Save Figure S1
ggsave("figure_S1.pdf",
       plot = comm_plot,
       device = "pdf",
       width = 8.27,
       height = 11.69,
       dpi = 300)


# -- FIGURE 4: Ranking-based binarization --
# Creating an empty data frame to store results
eval_rbb <- data.frame(
  "ID" = rep(0, 100 * 4463),
  "richness" = rep(0, 100 * 4463),
  "CM_LNC" = rep(0, 100 * 4463),
  "CM_SLA" = rep(0, 100 * 4463),
  "CM_PLH" = rep(0, 100 * 4463),
  "CSTD_LNC" = rep(0, 100 * 4463),
  "CSTD_SLA" = rep(0, 100 * 4463),
  "CSTD_PLH" = rep(0, 100 * 4463),
  "richness_ssdm" = rep(0, 100 * 4463),
  "richness_camelia_mem" = rep(0, 100 * 4463),
  "richness_camelia_obs" = rep(0, 100 * 4463),
  "richness_camelia_proxy" = rep(0, 100 * 4463),
  "Precision_ssdm" = rep(0, 100 * 4463),
  "Precision_camelia_mem" = rep(0, 100 * 4463),
  "Precision_camelia_obs" = rep(0, 100 * 4463),
  "Precision_camelia_proxy" = rep(0, 100 * 4463),
  "Sorensen_ssdm" = rep(0, 100 * 4463),
  "Sorensen_camelia_mem" = rep(0, 100 * 4463),
  "Sorensen_camelia_obs" = rep(0, 100 * 4463),
  "Sorensen_camelia_proxy" = rep(0, 100 * 4463),
  "Recall_ssdm" = rep(0, 100 * 4463),
  "Recall_camelia_mem" = rep(0, 100 * 4463),
  "Recall_camelia_obs" = rep(0, 100 * 4463),
  "Recall_camelia_proxy" = rep(0, 100 * 4463),
  "Specificite_ssdm" = rep(0, 100 * 4463),
  "Specificite_camelia_mem" = rep(0, 100 * 4463),
  "Specificite_camelia_obs" = rep(0, 100 * 4463),
  "Specificite_camelia_proxy" = rep(0, 100 * 4463),
  "CM_LNC_ssdm" = rep(0, 100 * 4463),
  "CM_SLA_ssdm" = rep(0, 100 * 4463),
  "CM_PLH_ssdm" = rep(0, 100 * 4463),
  "CM_LNC_camelia_mem" = rep(0, 100 * 4463),
  "CM_SLA_camelia_mem" = rep(0, 100 * 4463),
  "CM_PLH_camelia_mem" = rep(0, 100 * 4463),
  "CM_LNC_camelia_obs" = rep(0, 100 * 4463),
  "CM_SLA_camelia_obs" = rep(0, 100 * 4463),
  "CM_PLH_camelia_obs" = rep(0, 100 * 4463),
  "CM_LNC_camelia_proxy" = rep(0, 100 * 4463),
  "CM_SLA_camelia_proxy" = rep(0, 100 * 4463),
  "CM_PLH_camelia_proxy" = rep(0, 100 * 4463),
  "CSTD_LNC_ssdm" = rep(0, 100 * 4463),
  "CSTD_SLA_ssdm" = rep(0, 100 * 4463),
  "CSTD_PLH_ssdm" = rep(0, 100 * 4463),
  "CSTD_LNC_camelia_mem" = rep(0, 100 * 4463),
  "CSTD_SLA_camelia_mem" = rep(0, 100 * 4463),
  "CSTD_PLH_camelia_mem" = rep(0, 100 * 4463),
  "CSTD_LNC_camelia_obs" = rep(0, 100 * 4463),
  "CSTD_SLA_camelia_obs" = rep(0, 100 * 4463),
  "CSTD_PLH_camelia_obs" = rep(0, 100 * 4463),
  "CSTD_LNC_camelia_proxy" = rep(0, 100 * 4463),
  "CSTD_SLA_camelia_proxy" = rep(0, 100 * 4463),
  "CSTD_PLH_camelia_proxy" = rep(0, 100 * 4463)
)

for (i in 1:4463) {
  obs <- as.numeric(survey_obs[i, ])
  # Calculate richness and community metrics for each species
  eval_rbb$ID[((i-1)*100+1):((i-1)*100+100)] = rownames(survey_obs)[i]
  eval_rbb$richness[((i-1)*100+1):((i-1)*100+100)] = sum(obs == 1)
  
  # Calculate community metrics (CM) for each trait: LNC, SLA, PLH
  eval_rbb$CM_LNC[((i-1)*100+1):((i-1)*100+100)] = sum(obs * traits$LNC, na.rm = TRUE) / sum(obs == 1)
  eval_rbb$CM_SLA[((i-1)*100+1):((i-1)*100+100)] = sum(obs * traits$SLA, na.rm = TRUE) / sum(obs == 1)
  eval_rbb$CM_PLH[((i-1)*100+1):((i-1)*100+100)] = sum(obs * traits$PLH, na.rm = TRUE) / sum(obs == 1)
  
  # Calculate community standard deviations (CSTD) for each trait
  eval_rbb$CSTD_LNC[((i-1)*100+1):((i-1)*100+100)] = sqrt(sum(obs * (traits$LNC - eval_rbb$CM_LNC[(i-1)*100+1])^2, na.rm = TRUE) / sum(obs == 1))
  eval_rbb$CSTD_SLA[((i-1)*100+1):((i-1)*100+100)] = sqrt(sum(obs * (traits$SLA - eval_rbb$CM_SLA[(i-1)*100+1])^2, na.rm = TRUE) / sum(obs == 1))
  eval_rbb$CSTD_PLH[((i-1)*100+1):((i-1)*100+100)] = sqrt(sum(obs * (traits$PLH - eval_rbb$CM_PLH[(i-1)*100+1])^2, na.rm = TRUE) / sum(obs == 1))
  
  # Predictions
  predictions_sdm <- as.numeric(sdm_raw[i, ])
  predictions_camelia_mem <- as.numeric(camelia_mem_raw[i, ])
  predictions_camelia_obs <- as.numeric(camelia_obs_raw[i, ])
  predictions_camelia_proxy <- as.numeric(camelia_proxy_raw[i, ])
  
  # Loop over thresholds (from 0.01 to 1)
  for (j in seq(0.01, 1, by = 0.01)) {
    # Apply threshold and binarize predictions
    predictions_sdm_j <- as.numeric(predictions_sdm > j)
    predictions_camelia_mem_j <- as.numeric(predictions_camelia_mem > j)
    predictions_camelia_obs_j <- as.numeric(predictions_camelia_obs > j)
    predictions_camelia_proxy_j <- as.numeric(predictions_camelia_proxy > j)
    
    # Calculate performance metrics: Precision, Sorensen, Recall, Specificity for each model
    eval_rbb$Precision_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j == 1 & obs == 1) / sum(predictions_sdm_j == 1)
    eval_rbb$Sorensen_ssdm[(i-1)*100 + j*100] = 2 * sum(predictions_sdm_j == 1 & obs == 1) / (sum(predictions_sdm_j == 1) + sum(obs == 1))
    eval_rbb$Recall_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j == 1 & obs == 1) / sum(obs == 1)
    eval_rbb$Specificite_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j == 0 & obs == 0) / sum(predictions_sdm_j == 0)

    eval_rbb$Precision_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j == 1 & obs == 1) / sum(predictions_camelia_mem_j == 1)
    eval_rbb$Sorensen_camelia_mem[(i-1)*100 + j*100] = 2 * sum(predictions_camelia_mem_j == 1 & obs == 1) / (sum(predictions_camelia_mem_j == 1) + sum(obs == 1))
    eval_rbb$Recall_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j == 1 & obs == 1) / sum(obs == 1)
    eval_rbb$Specificite_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j == 0 & obs == 0) / sum(predictions_camelia_mem_j == 0)
    
    eval_rbb$Precision_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j == 1 & obs == 1) / sum(predictions_camelia_obs_j == 1)
    eval_rbb$Sorensen_camelia_obs[(i-1)*100 + j*100] = 2 * sum(predictions_camelia_obs_j == 1 & obs == 1) / (sum(predictions_camelia_obs_j == 1) + sum(obs == 1))
    eval_rbb$Recall_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j == 1 & obs == 1) / sum(obs == 1)
    eval_rbb$Specificite_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j == 0 & obs == 0) / sum(predictions_camelia_obs_j == 0)
    
    eval_rbb$Precision_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j == 1 & obs == 1) / sum(predictions_camelia_proxy_j == 1)
    eval_rbb$Sorensen_camelia_proxy[(i-1)*100 + j*100] = 2 * sum(predictions_camelia_proxy_j == 1 & obs == 1) / (sum(predictions_camelia_proxy_j == 1) + sum(obs == 1))
    eval_rbb$Recall_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j == 1 & obs == 1) / sum(obs == 1)
    eval_rbb$Specificite_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j == 0 & obs == 0) / sum(predictions_camelia_proxy_j == 0)
    
    # Calculate richness and CM for SDM predictions
    eval_rbb$richness_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j == 1)
    eval_rbb$CM_LNC_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j * traits$LNC, na.rm = TRUE) / sum(predictions_sdm_j == 1)
    eval_rbb$CM_SLA_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j * traits$SLA, na.rm = TRUE) / sum(predictions_sdm_j == 1)
    eval_rbb$CM_PLH_ssdm[(i-1)*100 + j*100] = sum(predictions_sdm_j * traits$PLH, na.rm = TRUE) / sum(predictions_sdm_j == 1)
    eval_rbb$CSTD_LNC_ssdm[(i-1)*100 + j*100] = sqrt(sum(predictions_sdm_j * (traits$LNC - eval_rbb$CM_LNC_ssdm[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_sdm_j == 1)))
    eval_rbb$CSTD_SLA_ssdm[(i-1)*100 + j*100] = sqrt(sum(predictions_sdm_j * (traits$SLA - eval_rbb$CM_SLA_ssdm[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_sdm_j == 1)))
    eval_rbb$CSTD_PLH_ssdm[(i-1)*100 + j*100] = sqrt(sum(predictions_sdm_j * (traits$PLH - eval_rbb$CM_PLH_ssdm[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_sdm_j == 1)))
    
    eval_rbb$richness_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j == 1)
    eval_rbb$CM_LNC_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j * traits$LNC, na.rm = TRUE) / sum(predictions_camelia_mem_j == 1)
    eval_rbb$CM_SLA_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j * traits$SLA, na.rm = TRUE) / sum(predictions_camelia_mem_j == 1)
    eval_rbb$CM_PLH_camelia_mem[(i-1)*100 + j*100] = sum(predictions_camelia_mem_j * traits$PLH, na.rm = TRUE) / sum(predictions_camelia_mem_j == 1)
    eval_rbb$CSTD_LNC_camelia_mem[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_mem_j * (traits$LNC - eval_rbb$CM_LNC_camelia_mem[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_mem_j == 1)))
    eval_rbb$CSTD_SLA_camelia_mem[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_mem_j * (traits$SLA - eval_rbb$CM_SLA_camelia_mem[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_mem_j == 1)))
    eval_rbb$CSTD_PLH_camelia_mem[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_mem_j * (traits$PLH - eval_rbb$CM_PLH_camelia_mem[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_mem_j == 1)))
    
    eval_rbb$richness_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j == 1)
    eval_rbb$CM_LNC_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j * traits$LNC, na.rm = TRUE) / sum(predictions_camelia_obs_j == 1)
    eval_rbb$CM_SLA_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j * traits$SLA, na.rm = TRUE) / sum(predictions_camelia_obs_j == 1)
    eval_rbb$CM_PLH_camelia_obs[(i-1)*100 + j*100] = sum(predictions_camelia_obs_j * traits$PLH, na.rm = TRUE) / sum(predictions_camelia_obs_j == 1)
    eval_rbb$CSTD_LNC_camelia_obs[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_obs_j * (traits$LNC - eval_rbb$CM_LNC_camelia_obs[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_obs_j == 1)))
    eval_rbb$CSTD_SLA_camelia_obs[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_obs_j * (traits$SLA - eval_rbb$CM_SLA_camelia_obs[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_obs_j == 1)))
    eval_rbb$CSTD_PLH_camelia_obs[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_obs_j * (traits$PLH - eval_rbb$CM_PLH_camelia_obs[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_obs_j == 1)))
    
    eval_rbb$richness_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j == 1)
    eval_rbb$CM_LNC_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j * traits$LNC, na.rm = TRUE) / sum(predictions_camelia_proxy_j == 1)
    eval_rbb$CM_SLA_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j * traits$SLA, na.rm = TRUE) / sum(predictions_camelia_proxy_j == 1)
    eval_rbb$CM_PLH_camelia_proxy[(i-1)*100 + j*100] = sum(predictions_camelia_proxy_j * traits$PLH, na.rm = TRUE) / sum(predictions_camelia_proxy_j == 1)
    eval_rbb$CSTD_LNC_camelia_proxy[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_proxy_j * (traits$LNC - eval_rbb$CM_LNC_camelia_proxy[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_proxy_j == 1)))
    eval_rbb$CSTD_SLA_camelia_proxy[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_proxy_j * (traits$SLA - eval_rbb$CM_SLA_camelia_proxy[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_proxy_j == 1)))
    eval_rbb$CSTD_PLH_camelia_proxy[(i-1)*100 + j*100] = sqrt(sum(predictions_camelia_proxy_j * (traits$PLH - eval_rbb$CM_PLH_camelia_proxy[(i-1)*100+j*100])^2, na.rm = TRUE) / length(which(predictions_camelia_proxy_j == 1)))
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
  scale_y_continuous(trans = log10_trans(), labels = scales::comma, limits = c(0.0000005, 10))+
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
# Create a dataframe to store evaluation metrics for each site
figure_sites_opth <- data.frame(
  ID = rownames(survey_obs),
  richness = rowSums(survey_obs),
  richness_sdm_opth = rowSums(sdm_opth),
  richness_camelia_mem_opth = rowSums(camelia_mem_opth),
  richness_camelia_obs_opth = rowSums(camelia_obs_opth),
  richness_camelia_proxy_opth = rowSums(camelia_proxy_opth)
)

# Define models to evaluate
models <- list(
  "sdm_opth" = sdm_opth,
  "camelia_mem_opth" = camelia_mem_opth,
  "camelia_obs_opth" = camelia_obs_opth,
  "camelia_proxy_opth" = camelia_proxy_opth
)

# Loop through models and compute performance metrics
for (model_name in names(models)) {
  pred <- models[[model_name]]
  
  VP <- rowSums((pred == 1) & (survey_obs == 1)) # True Positives
  FP <- rowSums((pred == 1) & (survey_obs == 0)) # False Positives
  VN <- rowSums((pred == 0) & (survey_obs == 0)) # True Negatives
  FN <- rowSums((pred == 0) & (survey_obs == 1)) # False Negatives
  
  # Add results to the dataframe
  figure_sites_opth[[paste0("VP_", model_name)]] <- VP
  figure_sites_opth[[paste0("FP_", model_name)]] <- FP
  figure_sites_opth[[paste0("VN_", model_name)]] <- VN
  figure_sites_opth[[paste0("FN_", model_name)]] <- FN
  figure_sites_opth[[paste0("Sensitivity_", model_name)]] <- VP / (VP + FN) # Recall
  figure_sites_opth[[paste0("Specificity_", model_name)]] <- VN / (VN + FP)
  figure_sites_opth[[paste0("Precision_", model_name)]] <- VP / (VP + FP)
  figure_sites_opth[[paste0("Sorensen_", model_name)]] <- 2 * VP / (2 * VP + FP + FN)
}

# Merge community indices data
figure_sites_opth <- left_join(figure_sites_opth, CM_tot_opth, by = "ID")
colnames(figure_sites_opth)[(ncol(figure_sites_opth) - length(colnames(CM_tot_opth)) + 2):ncol(figure_sites_opth)] <- paste0("CM_", colnames(figure_sites_opth)[(ncol(figure_sites_opth) - length(colnames(CM_tot_opth)) + 2):ncol(figure_sites_opth)])

figure_sites_opth <- left_join(figure_sites_opth, CSTD_tot_opth, by = "ID")
colnames(figure_sites_opth)[(ncol(figure_sites_opth) - length(colnames(CSTD_tot_opth)) + 2):ncol(figure_sites_opth)] <- paste0("CSTD_", colnames(figure_sites_opth)[(ncol(figure_sites_opth) - length(colnames(CSTD_tot_opth)) + 2):ncol(figure_sites_opth)])

data_richness_precision_opth <- data.frame(
  ID = rep(figure_sites_opth$ID, 4),
  `log(rich_pred/rich_obs)` = c(
    log(figure_sites_opth$richness_sdm_opth / figure_sites_opth$richness),
    log(figure_sites_opth$richness_camelia_mem_opth / figure_sites_opth$richness),
    log(figure_sites_opth$richness_camelia_obs_opth / figure_sites_opth$richness),
    log(figure_sites_opth$richness_camelia_proxy_opth / figure_sites_opth$richness)
  ),
  Precision = c(
    figure_sites_opth$Precision_sdm_opth,
    figure_sites_opth$Precision_camelia_mem_opth,
    figure_sites_opth$Precision_camelia_obs_opth,
    figure_sites_opth$Precision_camelia_proxy_opth
  ),
  Sorensen = c(
    figure_sites_opth$Sorensen_sdm_opth,
    figure_sites_opth$Sorensen_camelia_mem_opth,
    figure_sites_opth$Sorensen_camelia_obs_opth,
    figure_sites_opth$Sorensen_camelia_proxy_opth
  ),
  `PLH abs((CSTD obs - CSTD pred)/CSTD obs)` = c(
    abs(figure_sites_opth$CSTD_PLH_obs - figure_sites_opth$CSTD_PLH_sdm_opth) / figure_sites_opth$CSTD_PLH_obs,
    abs(figure_sites_opth$CSTD_PLH_obs - figure_sites_opth$CSTD_PLH_camelia_mem_opth) / figure_sites_opth$CSTD_PLH_obs,
    abs(figure_sites_opth$CSTD_PLH_obs - figure_sites_opth$CSTD_PLH_camelia_obs_opth) / figure_sites_opth$CSTD_PLH_obs,
    abs(figure_sites_opth$CSTD_PLH_obs - figure_sites_opth$CSTD_PLH_camelia_proxy_opth) / figure_sites_opth$CSTD_PLH_obs
  ),
  `SLA abs((CSTD obs - CSTD pred)/CSTD obs)` = c(
    abs(figure_sites_opth$CSTD_SLA_obs - figure_sites_opth$CSTD_SLA_sdm_opth) / figure_sites_opth$CSTD_SLA_obs,
    abs(figure_sites_opth$CSTD_SLA_obs - figure_sites_opth$CSTD_SLA_camelia_mem_opth) / figure_sites_opth$CSTD_SLA_obs,
    abs(figure_sites_opth$CSTD_SLA_obs - figure_sites_opth$CSTD_SLA_camelia_obs_opth) / figure_sites_opth$CSTD_SLA_obs,
    abs(figure_sites_opth$CSTD_SLA_obs - figure_sites_opth$CSTD_SLA_camelia_proxy_opth) / figure_sites_opth$CSTD_SLA_obs
  ),
  `LNC abs((CSTD obs - CSTD pred)/CSTD obs)` = c(
    abs(figure_sites_opth$CSTD_LNC_obs - figure_sites_opth$CSTD_LNC_sdm_opth) / figure_sites_opth$CSTD_LNC_obs,
    abs(figure_sites_opth$CSTD_LNC_obs - figure_sites_opth$CSTD_LNC_camelia_mem_opth) / figure_sites_opth$CSTD_LNC_obs,
    abs(figure_sites_opth$CSTD_LNC_obs - figure_sites_opth$CSTD_LNC_camelia_obs_opth) / figure_sites_opth$CSTD_LNC_obs,
    abs(figure_sites_opth$CSTD_LNC_obs - figure_sites_opth$CSTD_LNC_camelia_proxy_opth) / figure_sites_opth$CSTD_LNC_obs
  ),
  `PLH abs((CM obs - CM pred)/CM obs)` = c(
    abs(figure_sites_opth$CM_PLH_obs - figure_sites_opth$CM_PLH_sdm_opth) / figure_sites_opth$CM_PLH_obs,
    abs(figure_sites_opth$CM_PLH_obs - figure_sites_opth$CM_PLH_camelia_mem_opth) / figure_sites_opth$CM_PLH_obs,
    abs(figure_sites_opth$CM_PLH_obs - figure_sites_opth$CM_PLH_camelia_obs_opth) / figure_sites_opth$CM_PLH_obs,
    abs(figure_sites_opth$CM_PLH_obs - figure_sites_opth$CM_PLH_camelia_proxy_opth) / figure_sites_opth$CM_PLH_obs
  ),
  `SLA abs((CM obs - CM pred)/CM obs)` = c(
    abs(figure_sites_opth$CM_SLA_obs - figure_sites_opth$CM_SLA_sdm_opth) / figure_sites_opth$CM_SLA_obs,
    abs(figure_sites_opth$CM_SLA_obs - figure_sites_opth$CM_SLA_camelia_mem_opth) / figure_sites_opth$CM_SLA_obs,
    abs(figure_sites_opth$CM_SLA_obs - figure_sites_opth$CM_SLA_camelia_obs_opth) / figure_sites_opth$CM_SLA_obs,
    abs(figure_sites_opth$CM_SLA_obs - figure_sites_opth$CM_SLA_camelia_proxy_opth) / figure_sites_opth$CM_SLA_obs
  ),
  `LNC abs((CM obs - CM pred)/CM obs)` = c(
    abs(figure_sites_opth$CM_LNC_obs - figure_sites_opth$CM_LNC_sdm_opth) / figure_sites_opth$CM_LNC_obs,
    abs(figure_sites_opth$CM_LNC_obs - figure_sites_opth$CM_LNC_camelia_mem_opth) / figure_sites_opth$CM_LNC_obs,
    abs(figure_sites_opth$CM_LNC_obs - figure_sites_opth$CM_LNC_camelia_obs_opth) / figure_sites_opth$CM_LNC_obs,
    abs(figure_sites_opth$CM_LNC_obs - figure_sites_opth$CM_LNC_camelia_proxy_opth) / figure_sites_opth$CM_LNC_obs
  ),
  Models = rep(c("SSDM", "CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY"), each = nrow(figure_sites_opth))
)

data_richness_precision_opth$Models <- factor(data_richness_precision_opth$Models, levels = c("SSDM", "CAMELIA MEM", "CAMELIA OBS", "CAMELIA PROXY"))

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
                                                           rmse(SR_tot_opth$SR_obs, SR_tot_opth$SR_camelia_proxy_opth)),
                       "Binarised by PRR R2"=c(cor(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_mem)^2,
                                                         cor(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_sdm_prr)^2,
                                                         cor(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_camelia_mem_prr)^2,
                                                         cor(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_camelia_obs_prr)^2,
                                                         cor(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_camelia_proxy_prr)^2,
                                                         cor(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_mem)^2,
                                                         cor(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_sdm_prr)^2,
                                                         cor(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_camelia_mem_prr)^2,
                                                         cor(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_camelia_obs_prr)^2,
                                                         cor(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_camelia_proxy_prr)^2,
                                                         cor(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_mem)^2,
                                                         cor(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_sdm_prr)^2,
                                                         cor(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_camelia_mem_prr)^2,
                                                         cor(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_camelia_obs_prr)^2,
                                                         cor(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_camelia_proxy_prr)^2,
                                                         
                                                         cor(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_mem)^2,
                                                         cor(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_sdm_prr)^2,
                                                         cor(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_camelia_mem_prr)^2,
                                                         cor(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_camelia_obs_prr)^2,
                                                         cor(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_camelia_proxy_prr)^2,
                                                         cor(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_mem)^2,
                                                         cor(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_sdm_prr)^2,
                                                         cor(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_camelia_mem_prr)^2,
                                                         cor(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_camelia_obs_prr)^2,
                                                         cor(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_camelia_proxy_prr)^2,
                                                         cor(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_mem)^2,
                                                         cor(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_sdm_prr)^2,
                                                         cor(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_camelia_mem_prr)^2,
                                                         cor(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_camelia_obs_prr)^2,
                                                         cor(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_camelia_proxy_prr)^2,
                                                         
                                                         cor(SR_tot_prr$SR_obs, SR_tot_prr$SR_mem)^2,
                                                         cor(SR_tot_prr$SR_obs, SR_tot_prr$SR_sdm_prr)^2,
                                                         cor(SR_tot_prr$SR_obs, SR_tot_prr$SR_camelia_mem_prr)^2,
                                                         cor(SR_tot_prr$SR_obs, SR_tot_prr$SR_camelia_obs_prr)^2,
                                                         cor(SR_tot_prr$SR_obs, SR_tot_prr$SR_camelia_proxy_prr)^2),  
                       "Binarised by PRR RMSE"=c(rmse(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_mem),
                                                           rmse(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_sdm_prr),
                                                           rmse(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_camelia_mem_prr),
                                                           rmse(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_camelia_obs_prr),
                                                           rmse(CM_tot_prr$PLH_obs, CM_tot_prr$PLH_camelia_proxy_prr),
                                                           rmse(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_mem),
                                                           rmse(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_sdm_prr),
                                                           rmse(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_camelia_mem_prr),
                                                           rmse(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_camelia_obs_prr),
                                                           rmse(CM_tot_prr$SLA_obs, CM_tot_prr$SLA_camelia_proxy_prr),
                                                           rmse(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_mem),
                                                           rmse(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_sdm_prr),
                                                           rmse(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_camelia_mem_prr),
                                                           rmse(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_camelia_obs_prr),
                                                           rmse(CM_tot_prr$LNC_obs, CM_tot_prr$LNC_camelia_proxy_prr),
                                                           
                                                           rmse(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_mem),
                                                           rmse(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_sdm_prr),
                                                           rmse(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_camelia_mem_prr),
                                                           rmse(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_camelia_obs_prr),
                                                           rmse(CSTD_tot_prr$PLH_obs, CSTD_tot_prr$PLH_camelia_proxy_prr),
                                                           rmse(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_mem),
                                                           rmse(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_sdm_prr),
                                                           rmse(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_camelia_mem_prr),
                                                           rmse(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_camelia_obs_prr),
                                                           rmse(CSTD_tot_prr$SLA_obs, CSTD_tot_prr$SLA_camelia_proxy_prr),
                                                           rmse(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_mem),
                                                           rmse(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_sdm_prr),
                                                           rmse(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_camelia_mem_prr),
                                                           rmse(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_camelia_obs_prr),
                                                           rmse(CSTD_tot_prr$LNC_obs, CSTD_tot_prr$LNC_camelia_proxy_prr),
                                                           
                                                           rmse(SR_tot_prr$SR_obs, SR_tot_prr$SR_mem),
                                                           rmse(SR_tot_prr$SR_obs, SR_tot_prr$SR_sdm_prr),
                                                           rmse(SR_tot_prr$SR_obs, SR_tot_prr$SR_camelia_mem_prr),
                                                           rmse(SR_tot_prr$SR_obs, SR_tot_prr$SR_camelia_obs_prr),
                                                           rmse(SR_tot_prr$SR_obs, SR_tot_prr$SR_camelia_proxy_prr)))
write.csv(table_S2, file="table_S2.csv")
# Species predictions evaluation : 
## RAW
AUC_raw <- data.frame("Species"=rep(colnames(survey_obs),4), "AUC"=rep(0,831*4), "pr-AUC"=rep(0,831*4), "Approach"=c(rep("SSDM",831), rep("CAMELIA MEM",831), rep("CAMELIA OBS",831), rep("CAMELIA PROXY", 831)))
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

table_S4 <- data.frame("Stacking processus"=c(rep("Probability", 2), rep("Binarised by TSS threshold", 4), rep("Binarised by PRR",4)),
                       "Indices" = c("mean AUC", "mean pr-AUC", rep(c("mean TSS", "mean sensibility", "mean specificity", "mean precision"),2)),
                       "S-SDM"=c(mean(AUC_raw$AUC[which(AUC_raw$Approach=="SSDM")]), 
                                 mean(AUC_raw$pr.AUC[which(AUC_raw$Approach=="SSDM")]),
                                 mean(eval_species_opth$TSS_SSDM, na.rm = TRUE),
                                 mean(eval_species_opth$Sensibilite_SSDM, na.rm = TRUE),
                                 mean(eval_species_opth$Specificite_SSDM, na.rm = TRUE),
                                 mean(eval_species_opth$Precision_SSDM, na.rm = TRUE),
                                 mean(eval_species_prr$TSS_SSDM, na.rm = TRUE),
                                 mean(eval_species_prr$Sensibilite_SSDM, na.rm = TRUE),
                                 mean(eval_species_prr$Specificite_SSDM, na.rm = TRUE),
                                 mean(eval_species_prr$Precision_SSDM, na.rm = TRUE)),
                       "CAMELIA MEM"=c(mean(AUC_raw$AUC[which(AUC_raw$Approach=="CAMELIA MEM")]), 
                                 mean(AUC_raw$pr.AUC[which(AUC_raw$Approach=="CAMELIA MEM")]),
                                 mean(eval_species_opth$TSS_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_opth$Sensibilite_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_opth$Specificite_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_opth$Precision_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_prr$TSS_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_prr$Sensibilite_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_prr$Specificite_CAMELIA_MEM, na.rm = TRUE),
                                 mean(eval_species_prr$Precision_CAMELIA_MEM, na.rm = TRUE)),
                       "CAMELIA OBS"=c(mean(AUC_raw$AUC[which(AUC_raw$Approach=="CAMELIA OBS")]), 
                                       mean(AUC_raw$pr.AUC[which(AUC_raw$Approach=="CAMELIA OBS")]),
                                       mean(eval_species_opth$TSS_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_opth$Sensibilite_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_opth$Specificite_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_opth$Precision_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_prr$TSS_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_prr$Sensibilite_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_prr$Specificite_CAMELIA_OBS, na.rm = TRUE),
                                       mean(eval_species_prr$Precision_CAMELIA_OBS, na.rm = TRUE)),
                       "CAMELIA PROXY"=c(mean(AUC_raw$AUC[which(AUC_raw$Approach=="CAMELIA PROXY")]), 
                                       mean(AUC_raw$pr.AUC[which(AUC_raw$Approach=="CAMELIA PROXY")]),
                                       mean(eval_species_opth$TSS_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_opth$Sensibilite_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_opth$Specificite_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_opth$Precision_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_prr$TSS_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_prr$Sensibilite_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_prr$Specificite_CAMELIA_PROXY, na.rm = TRUE),
                                       mean(eval_species_prr$Precision_CAMELIA_PROXY, na.rm = TRUE))
                       
)

write.csv(table_S4, file="table_S4.csv")