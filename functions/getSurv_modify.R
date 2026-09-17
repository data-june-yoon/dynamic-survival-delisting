
getsurv <- function(obj, times)
{
    # get the survival probability for times from KM curve `obj'
    
    #if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
    # <FIXME: methods may have problems with that>
    #obj = pred[II]; times = tpntLod[r.ID == jall[j]]
  
    class(obj) <- NULL
    # </FIXME>
    lt <- length(times)
    nsurv <- times
    
    # if the times are the same, return the km-curve
    
    if(length(times) == length(obj$time)) {
        if (all(times == obj$time)) return(obj$surv)
    }

    # otherwise get the km-value for every element of times separatly
    
    inside <- times %in% obj$time
    for (i in (1:lt)) {
        if (inside[i]){
            nsurv[i] <- obj$surv[obj$time == times[i]]
        }   else  {
            less <- obj$time[obj$time < times[i]]
            if (length(less) == 0) {
                nsurv[i] <- 1
            } else {
                nsurv[i] <- obj$surv[obj$time == max(less)]
            } # end else
        } # end else
    }
    nsurv
}
#-------------------------------
getsurv_pcox <- function(obj_surv, obj_time, times)
{
  # get the survival probability for times from KM curve `obj'
  
  #if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
  # <FIXME: methods may have problems with that>
  #obj = pred[II]; times = tpntLod[r.ID == jall[j]]
  
  #class(obj) <- NULL
  # </FIXME>
  lt <- length(times)
  nsurv <- times
  
  # if the times are the same, return the km-curve
  
  if(length(times) == length(obj_time)) {
    if (all(times == obj_time)) return(obj_surv)
  }
  
  # otherwise get the km-value for every element of times separatly
  inside <- times %in% obj_time
  for (i in (1:lt)) {
    if (inside[i]){
      nsurv[i] <- obj_surv[obj_time == times[i]]
    }   else  {
      less <- obj_time[obj_time < times[i]]
      if (length(less) == 0) {
        nsurv[i] <- 1
      } else {
        nsurv[i] <- obj_surv[obj_time == max(less)]
      } # end else
    } # end else
  }
  nsurv
}
