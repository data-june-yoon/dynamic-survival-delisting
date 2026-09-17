predictProb.ltrccif <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                time.eval, time.tau = NULL){
    
    pred <- partykit::predict.cforest(object = object, newdata = newdata, OOB = OOB, type = "prob",
                                      FUN = .pred_Surv_nolog)
    xvar.names <- attr(object$terms,"term.labels")
    yvar.names <- as.character(object$formulaLTRC[[2]])[2:4]
    idname <- "id"
    
    # missing values can be present in the prediction
    if (is.null(newdata) || OOB){
        # first column: Surv(tleft,tright,event), second column: (id)
        newdata <- as.data.frame(as.matrix(object$data[, c(1, ncol(object$data)), drop = FALSE]))
        names(newdata) = c(yvar.names, idname)
    } else {
        if (missing(newdata.id)){
            newdata$id <- 1:nrow(newdata)
        } else {
            names(newdata)[names(newdata) == deparse(substitute(newdata.id))] <- idname
        }
        newdata <- as.data.frame(newdata[, c(yvar.names, idname)])
    }
    
    rm(object)
    
    
    N <- length(unique(newdata[, "id"])) # number of subjects
    
    if (is.null(time.tau)){
        time.tau <- rep(max(time.eval), N)
    } else {
        if (N != length(time.tau)) stop("time.tau should be a vector of length equaling to number of SUBJECT observation! \n
                                     In the time-varying case, check whether newdata.id has been correctly specified!")
    }
    
    Shat <- sapply(1:N, function(Ni) .shatfunc(Ni, data = newdata, pred = pred, tpnt = time.eval, tau = time.tau))
    obj <- Surv(newdata[, yvar.names[1]],
                newdata[, yvar.names[2]],
                newdata[, yvar.names[3]])
    RES <- list(survival.probs = Shat,
                survival.times = time.eval,
                survival.tau = time.tau,
                survival.obj = obj,
                survival.id = newdata$id)
    rm(newdata)
    rm(Shat)
    rm(time.eval)
    rm(time.tau)
    return(RES)
}