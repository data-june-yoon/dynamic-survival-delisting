


calc_Brier_At_Times_Clean <- function(LTRCCIF_model, eval_times, newdata) {
    # Internally add Time 0 to prevent IPCW length 0 error
    eval_times_with_zero <- c(0, eval_times)
    
    Survobj_full <- Surv(newdata$Start, newdata$Stop, newdata$Event)
    Predobj <- predictProb(object = LTRCCIF_model, 
                           newdata = newdata, 
                           newdata.id = ID, 
                           time.eval = eval_times_with_zero)
    
    bs_res <- sbrier_ltrc(obj = Survobj_full, id = newdata$ID, pred = Predobj, type = "BS")
    
    return(bs_res$BScore[-1]) # Drop Time 0 before returning
}

#---------------------------------------------------------

calc_Brier_At_Times_Penalized <- function(model_fit, eval_times, newdata) {
    # Internally add Time 0 to prevent IPCW length 0 error
    eval_times_with_zero <- c(0, eval_times)
    
    pred_obj <- pcoxsurvfit(model_fit, newdata)
    prob_matrix <- shat_pcox(newdata, pred = pred_obj, tpnt = eval_times)
    
    N_unique <- ncol(prob_matrix)
    ones_row <- matrix(1, nrow = 1, ncol = N_unique)
    prob_matrix_extended <- rbind(ones_row, prob_matrix)
    
    pred_list <- list(
        survival.probs = prob_matrix_extended,
        survival.times = eval_times_with_zero,
        survival.tau = rep(max(eval_times_with_zero), N_unique)
    )
    
    Survobj_full <- Surv(newdata$Start, newdata$Stop, newdata$Event)
    bs_res <- sbrier_ltrc(obj = Survobj_full, id = newdata$ID, pred = pred_list, type = "BS")
    
    return(bs_res$BScore[-1]) # Drop Time 0 before returning
}

#----------------------------------------------------

calc_BS_forests<-function(model_fit, newdata){
    
    Tpnt <- seq(0, max(newdata$Stop), length.out = 50)
    Survobj<-Surv(newdata$Start, newdata$Stop, newdata$Event)
    pred <- predictProb(model_fit, newdata = newdata, newdata.id = ID, time.eval = Tpnt)
    
    BS<-sbrier_ltrc(obj = Survobj, id = newdata$ID, pred = pred, type = "BS")$BScore
    IBS<-sbrier_ltrc(obj = Survobj, id = newdata$ID, pred = pred, type = "IBS")
    print(paste0('IBS = ', IBS))
    return(list(pred = pred, IBS = IBS, BS = BS, Tpnt=Tpnt))
} # end function calc_BS_forests 


#---------------------------------------------------

calc_BS_pcox<-function(model, model_fit, alpha, lambda, newdata){
    
    Tpnt <- seq(0, max(newdata$Stop), length.out = 50)
    Survobj<-Surv(newdata$Start, newdata$Stop, newdata$Event)
    
    pcox_pred<- pcoxsurvfit(model_fit, newdata) # This function creates survival curves
    
    # shat does not work with 'pcox_pred' so modified
    pcox_shat = shat_pcox(newdata, pred = pcox_pred, tpnt = Tpnt)
    pred = list(survival.probs = pcox_shat,
                survival.times = Tpnt,
                survival.tau = rep(max(Tpnt), length(unique(newdata$ID))))
    
    BS<-sbrier_ltrc(obj = Survobj, id = newdata$ID, pred = pred, type = "BS")$BScore
    IBS<-sbrier_ltrc(obj = Survobj, id = newdata$ID, pred = pred, type = "IBS")
    print(paste0(model, '_IBS = ', IBS))
    return(list(pred = pred, IBS = IBS, BS = BS, Tpnt = Tpnt))
} # end function calc_BS_pcox


