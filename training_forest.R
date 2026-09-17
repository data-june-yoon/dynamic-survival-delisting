
#-----------training cif model on training dataset --------
# ltrccif function does mtry tuning and returns the best model 

rm(list=ls())

source('R/setup.R')
set.seed(seedn)

#----- load training / test datasets ---------
#load(file = paste0('Data/win_imp_sc_surv_', filename, '.Rdata')) # train_surv, test_surv


ntree_v<-c(100, 300, 500)
for(i in 1:length(ntree_v)){
    ntree_n<-ntree_v[i]
    print(paste0('===== ntree = ', ntree_n, '==========='))
    cif_fit = new_ltrccif(formula = Formula, data = train_surv, id = ID, ntree = ntree_n, mtry = NULL) 
    print(paste0('tuned mtry = ', cif_fit$ret$info$control$mtry))
    save(cif_fit, file = paste0('Data/cif_', filename, "_ntree", ntree_n, '.Rdata'))
}







