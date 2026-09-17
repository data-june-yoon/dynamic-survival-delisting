
source('R/setup.R')
source('R/functions_evaluation.R')

#-------------------------

# loading data sets and model fits 
startyr<-2000; endyr<-2022
filename<-paste0('data_', startyr, '_', endyr)
load(file = paste0('Data/win_imp_sc_surv_', filename, '.Rdata')) # train_surv, test_surv

# Load Penalized Models
load(file=paste0('Data/lasso5cv_', filename, '.Rdata'))
load(file=paste0('Data/ridge5cv_', filename, '.Rdata'))
load(file=paste0('Data/elnet5cv_', filename, '.Rdata'))

# Load CIF Model
ntree_n <- 500
load(file = paste0('Data/cif_', filename, "_ntree", ntree_n, '.Rdata')) 
cifmodel_fit <- cif_fit$ret

# Setup Landmark Times
units <- 365.25
#eval_years <- c(3, 5, 7, 9, 11)
eval_years<-3:12
eval_times <- eval_years * units 

# Initialize Results Matrix
tableBS_TV <- matrix(nrow = 4, ncol = length(eval_times))
rownames(tableBS_TV) <- c(paste0("CIF(ntree = ", ntree_n, ")"), "Lasso", "Ridge", "Elnet")
colnames(tableBS_TV) <- eval_years

# ==============================================================================
# calculate brier scores
# ==============================================================================

# --- CIF Model ---
print("Calculating BS for CIF...")
tableBS_TV[1, ] <- calc_Brier_At_Times_Clean(cifmodel_fit, eval_times, test_surv)


# --- Penalized Linear Models ---
models <- list(Lasso = lasso_fit, Ridge = ridge_fit, Elnet = elnet_fit)

for(m_name in names(models)) {
    print(paste("Calculating BS for", m_name, "..."))
    model_fit<-models[[m_name]]
    tableBS_TV[m_name, ] <- calc_Brier_At_Times_Penalized(model_fit, eval_times, test_surv)
}

# ==============================================================================
# final output
# ==============================================================================

print(round(tableBS_TV, 4))
save(tableBS_TV, file = paste0('Data/new_BS_TV_5cv', filename, '.Rdata'))
