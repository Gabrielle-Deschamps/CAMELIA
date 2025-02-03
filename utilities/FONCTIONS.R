#This script hold some useful functions to calculated community indices and binarised predicted community using PRR :

#### Packages
library(ggplot2)
library(tidyverse)
library(readxl)
library(Metrics)
library(ggh4x)
library(randomForestSRC)
library(pROC)


#### -- Data preparation -- ####
# This function returns an array of data that will be used to calculate CM and CSTD.
data_table <- function(data){
  data_table <- data.frame("ID"=rownames(data), data)
  colnames(data_table)<-c("ID", colnames(data))
  data_table <- data_table %>% pivot_longer(-ID)
  colnames(data_table)<- c("ID", "Taxa", "Abundance")
  data_table$Abundance <- as.numeric(data_table$Abundance)
  traits = data.frame("Taxa"=rownames(traits), traits)
  data_table <- left_join(data_table, traits, by="Taxa")
  data_table <- cbind(data_table, "A_LNC"=data_table$Abundance*data_table$LNC,
                   "A_SLA"=data_table$Abundance*data_table$SLA,
                   "A_PLH"=data_table$Abundance*data_table$PLH,
                   "pa_LNC"=data_table$Abundance*(data_table$LNC/data_table$LNC),
                   "pa_SLA"=data_table$Abundance*(data_table$SLA/data_table$SLA),
                   "pa_PLH"=data_table$Abundance*(data_table$PLH/data_table$PLH))
  data_table <- data_table %>% group_by(ID)
  return(data_table)
}
#### -- Calcul of the CM -- ####
# This function returns an array of data containing the CM for LNC, SLA and Plant Height.

CM <- function(data){
  site <- data.frame("ID"=rownames(data))
  data <- data_table(data)
  LNC = data %>%
    summarise("LNC"=c(sum(A_LNC, na.rm = TRUE)/sum(pa_LNC, na.rm = TRUE)))
  SLA = data %>% 
    summarise("SLA"=c(sum(A_SLA, na.rm = TRUE)/sum(pa_SLA, na.rm = TRUE)))
  PLH = data %>%
    summarise("PLH"=c(sum(A_PLH, na.rm = TRUE)/sum(pa_PLH, na.rm = TRUE)))
  CM = data.frame("ID"=LNC$ID, "LNC"=LNC$LNC, 'SLA'=SLA$SLA, "PLH"=PLH$PLH)
  CM = left_join(site, CM, by = "ID")
  rownames(CM)=CM$ID
  return(CM)
}


#### -- Calcul of the CSTD -- ####
# This function returns an array of data containing the CSTD for LNC, SLA and Plant Height.

CSTD <-  function(data){
  site <- data.frame("ID"=rownames(data))
  data <- data_table(data)
  LNC = data %>%
    summarise("LNC"=c(sqrt(sum(pa_LNC*(LNC-(sum(A_LNC,na.rm = TRUE)/sum(pa_LNC,na.rm = TRUE)))^2,na.rm = TRUE)/sum(pa_LNC,na.rm = TRUE))))
  SLA = data %>% 
    summarise("SLA"=c(sqrt(sum(pa_SLA*(SLA-(sum(A_SLA,na.rm = TRUE)/sum(pa_SLA,na.rm = TRUE)))^2,na.rm = TRUE)/sum(pa_SLA,na.rm = TRUE))))
  PLH = data %>%
    summarise("PLH"=c(sqrt(sum(pa_PLH*(PLH-(sum(A_PLH,na.rm = TRUE)/sum(pa_PLH,na.rm = TRUE)))^2,na.rm = TRUE)/sum(pa_PLH,na.rm = TRUE))))
  CSTD = data.frame("ID"=LNC$ID, "LNC"=LNC$LNC, 'SLA'=SLA$SLA, "PLH"=PLH$PLH)
  CSTD = left_join(site, CSTD, by = "ID")
  rownames(CSTD)=CSTD$ID
  return(CSTD)
}

#### -- PRR -- ####
prr <- function(Y_hat, richesse) {
  sites <- data.frame("ID"=rownames(Y_hat))
  richesse <- left_join(sites, richesse, by="ID")
  data_prr <- Y_hat
  for (s in 1:nrow(Y_hat)) {
    max <- length(which(Y_hat[s,] > 0))
    if (round(richesse[s,2]) < max) {
      sorted_row <- order(as.numeric(Y_hat[s,]), decreasing = TRUE)[1:round(richesse[s,2])]
    } else {
      sorted_row <- which(Y_hat[s,] > 0)
    }
    data_prr[s,] <- 0
    data_prr[s,sorted_row] <- 1
    print(s)
  }
  return(data_prr)
}

