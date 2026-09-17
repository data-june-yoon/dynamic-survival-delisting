
seedn<<-123

#-------------------------------
# library ready 
library(dplyr) 
library(ggplot2)
library(tidyr)
library(rsample)
library(survival)  # Surv function  
library(risksetROC)
library(purrr)
library(tibble)


#------- [installing LTRCforests package]
# (1) we need these packages before installing LTRCforest (need Rtolls installation)
#install.packages(c('ipred', 'prodlim', 'partykit'))

#(2) Rtools install => download tar.gz file from r cran archive and install the package 
library(partykit)
library(LTRCforests) # install from tar.gz. install ipred package first. Still a proper version of Rtools is needed for installation. 

#------ [installing pcoxtime]
# we need these packages before installation
# install.packages(c('doParallel', 'foreach', 'riskRegression', 'PermAlgo', 'pec', 'RcppArmadillo', 'remotes'))
# remotes::install_github("cran/pcoxtime") 

library(pcoxtime) # install from tar.gz

#----------- modified files ------  

source("R/getSurv_modify.R")
source("R/shat_modify.R")
source('R/predictProb_modify.R')
source('R/ltrccif_functions.R')
source('R/mindepth_varimp_functions.R')

#---------------------------------

ids<-c("permno", "hcomnam", "begdat", "public_date", "enddat", "event_1y", "event_2y", 
        "permco", "hsiccd", "hshrcd", "event_group_new")

ratio<-c("capital_ratio", "at_turn", "inv_turn", "rect_turn",
           "quick_ratio", "aftret_eq", "gpm", "opmad",
           "roa", "de_ratio", "debt_assets", "pe_inc") 

econ<-c("gdp_logR", "cpi_logR", "sp500_logR", "ff_diff") 

varset<-c(ratio, econ)

Formula<-formula(paste("Surv(Start, Stop, Event)~", paste(varset, collapse=" + ")))

startyr<-2000; endyr<-2022
filename<-paste0('data_', startyr, '_', endyr)

