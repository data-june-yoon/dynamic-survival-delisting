

library(prodlim)
library(truncdist) # rtrunc function 

#--------------------------------------------------------------

source('R/setup.R')
source('R/functions_eval_vimp_ranking.R')
source("R/functions_TV_PH_linear.R")

#--------------------------------------------------------------

set.seed(123)
vars<-c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8') # 8 covariates
Formula<-formula(paste("Surv(Start, Stop, Event)~", paste(vars, collapse=" + ")))

ntree_n<-500 
nsim<-100
nsample<-100  # the number of subjects
n_mtry<-3 # the default conditional inference forest : sqrt(6) = 2.44949

#--------------------------------------
ones<-c(1,1,1,1)

for(set in c(1, 2, 3)){
    # Beta order : (beta1-beta4, beta0)
    if(set == 1){
        Beta <- c(10*ones, 1)
    } else if(set == 2){
        Beta<-c(5*ones, 1)
    } else if(set == 3) {
        Beta<-c(1*ones, 1)
    } # end
    
    simData <- vector("list", nsim)
    
    for(i in 1:nsim){
        print(paste("set=", set, ", wi : i = ", i))
        wi_ret <- simple_PH_linear2(N = nsample, Beta)
        simData[[i]]<-wi_ret$fullData
    } # end i
    
    save(simData, file=paste0('Output/WI_rv_dat', set, '.Rdata'))
    
} # end set
#----------------------------
# running LTRC-CIF and calculating minimal depth on the simulated data 

for(set in 1:3){
    load(file=paste0('Output/WI_rv_dat', set, '.Rdata')) # simData
    min_d <- vector("list", nsim)
    for(i in 1:nsim){
        print(paste("set=", set, ", wi : i = ", i))
        wi_result<-new_ltrccif(formula = Formula, data = simData[[i]], id = ID, mtry = n_mtry, ntree = ntree_n)
        wi_md<-LTRCCIF_min_Depth(wi_result$ltrc_cif_obj, vars, ntree = ntree_n)
        min_d[[i]]<-wi_md
    } # end i
    
    save(min_d, file=paste0('Output/WI_rv_md', set, '.Rdata'))
    
} # end set
#----------------------------

for(set in 1:3){
    print(paste("set=", set))
    load(file=paste0('Output/WI_rv_md', set, '.Rdata')) # min_d
    r<-evaluate_mindepth_simulations(min_d, true_labels = c(1, 1, 1, 1, 0, 0, 0, 0))
    print(r)
} # end set




