
#----------- training regularized cox models (pcox) on training dataset
# (1) hyperparameter tuning (2) retraining on the entire training dataset 

rm(list=ls())

source('R/setup.R')
set.seed(seedn)

#----- load pre-processed training / test datasets ---------
load(file = paste0('Data/win_imp_sc_surv_', filename, '.Rdata')) # train_surv, test_surv

parallel::detectCores() # 20 

#----------------------------------
# 5 fold CV
n_cv<-5; c_clu<-5

### [penalized cox : LASSO]
print("==== LASSO tuning ======")

cv_pcox <- pcoxtimecv(Formula, data = train_surv, alphas = 1, lambdas = NULL, nlamdas = 20, 
                      nfolds = n_cv, nclusters = c_clu) 

alp_pcox<-1; lam_pcox <- cv_pcox$lambda.min
print(paste('LASSO : alpha = ', alp_pcox, ', lambda = ', lam_pcox))

lasso_fit <- pcoxtime(Formula, data = train_surv, alpha = alp_pcox, lambda = lam_pcox)

save(alp_pcox, lam_pcox, lasso_fit, file=paste0('Data/lasso5cv_', filename, '.Rdata'))

### [penalized cox : Ridge]
print("==== Ridge tuning ======")
cv_pcox <- pcoxtimecv(Formula, data = train_surv, alphas = 0, lambdas = NULL, nlamdas = 20, 
                      nfolds = n_cv, nclusters = c_clu) 

alp_pcox<-0; lam_pcox <- cv_pcox$lambda.min

print(paste('ridge : alpha = ', alp_pcox, ', lambda = ', lam_pcox))

ridge_fit <- pcoxtime(Formula, data = train_surv, alpha = alp_pcox, lambda = lam_pcox)

save(alp_pcox, lam_pcox, ridge_fit, file=paste0('Data/ridge5cv_', filename, '.Rdata'))


### [penalized cox : elastic net] (alpha, lambda all flexible)
print("==== elastic net tuning ======")
cv_pcox <- pcoxtimecv(Formula, data = train_surv, alphas = seq(from = 0.1, to = 0.8, by = 0.1), 
                      lambdas = NULL, nlamdas = 20,
                      nfolds = n_cv, nclusters = c_clu) 

alp_pcox <- cv_pcox$alpha.optimal
lam_pcox <- cv_pcox$lambda.min

print(paste('elestic net : alpha = ', alp_pcox, ', lambda = ', lam_pcox))

elnet_fit <- pcoxtime(Formula, data = train_surv, alpha = alp_pcox, lambda = lam_pcox)

save(alp_pcox, lam_pcox, elnet_fit, file=paste0('Data/elnet5cv_', filename, '.Rdata'))

