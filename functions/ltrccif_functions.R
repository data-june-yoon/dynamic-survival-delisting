# ltrccif function needed modification to store individual trees for min-depth based variable importance function. 

.logrank_trafo2 <- function(x2){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    unique.times <- unique(x2[,2][which(x2[, 3] == 1)])
    D <- rep(NA, length(unique.times))
    R <- rep(NA, length(unique.times))

    for(j in 1:length(unique.times)){
      D[j] = sum(unique.times[j] == x2[, 2])
    }

    for(k in 1:length(unique.times) ){
      value <- unique.times[k]
      R[k] <- sum(apply(x2[, 1:2], 1, function(interval){interval[1] < value & value <= interval[2]}))
    }

    Ratio <- D / R

    Ratio <- Ratio[order(unique.times)]
    Nelson.Aalen <- cumsum(Ratio)
    Event.time <- unique.times[order(unique.times)]
    Left <- sapply(x2[, 1], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
    Right <- sapply(x2[, 2], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})

    result<- x2[, 3] - (Right - Left)
  }
  return(as.double(result))
} # end function 

#--------------------------------------
#' Fit a LTRC conditional inference forest

new_ltrccif <- function(formula, data, id,
                    mtry = NULL, ntree = 100L,
                    bootstrap = c("by.sub","by.root","by.user","none"),
                    samptype = c("swor","swr"),
                    sampfrac = 0.632,
                    samp = NULL,
                    na.action = "na.omit",
                    stepFactor = 2,
                    trace = TRUE,
                    applyfun = NULL, cores = NULL,
                    control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                      minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                      minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                                      minprob = 0.01,
                                                      mincriterion = 0, saveinfo = FALSE)){

  #requireNamespace("inum")

  Call <- match.call()
  Call[[1]] <- as.name('ltrccif')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'data', 'id'),
                names(Call), nomatch = 0L)
  if (indx[1] == 0) stop("a formula argument is required")
  Call$formula <- eval(formula)

  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)

  mf <- eval.parent(temp)
  y <- model.extract(mf, 'response')
  if (!is.Surv(y)) stop("Response must be a survival object")
  if (!attr(y, "type") == "counting") stop("The Surv object must be of type 'counting'.")
  rm(y)

  # pull y-variable names
  yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]

  if (length(yvar.names) == 4){
    yvar.names = yvar.names[2:4]
  }

  Status <- data[, yvar.names[3]]
  # Times <- data[, yvar.names[2]]
  if (sum(Status) == 0) stop("All observations are right-censored with event = 0!")

  n <- nrow(data)

  ## if not specified, the first one will be used as default
  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)

  ## The following code to define id does not work since it could not handle missing values
  # id <- model.extract(mf, 'id')

  # this is a must, otherwise id cannot be passed to the next level in tune.ltrccif
  if (indx[3] == 0){
    ## If id is not present, then we add one more variable
    # mf$id <- 1:nrow(mf) ## Relabel
    data$id <- 1:n
  } else {
    ## If id is present, then we rename the column to be id
    # names(mf)[names(mf) == "(id)"] <- "id"
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }

  # extract the x-variable names
  xvar.names <- attr(terms(formula), 'term.labels')
  rm(temp)
  data <- data[, c("id", yvar.names, xvar.names)]

  ## bootstrap case
  if (length(data$id) == length(unique(data$id))){ # time-invariant LTRC data
    # it includes the case 1) when id = NULL, which is that id is not specified
    #                      2) when id is specified, but indeed LTRC time-invariant
    if (bootstrap == "by.sub") bootstrap <- "by.root"
  } else { # time-varying subject data
    id.sub <- unique(data$id)
    ## number of subjects
    n.sub <- length(id.sub)
  }

  if (samptype == "swor"){
    perturb = list(replace = FALSE, fraction = sampfrac)
  } else if (samptype == "swr"){
    perturb = list(replace = TRUE)
  } else {
    stop("samptype must set to be either 'swor' or 'swr'\n")
  }

  if (bootstrap == "by.sub"){
    size <- n.sub
    if (!perturb$replace) size <- floor(n.sub * perturb$fraction)
    samp <- replicate(ntree,
                      sample(id.sub, size = size,
                             replace = perturb$replace),
                      simplify = FALSE) # a list of length ntree
    samp <- lapply(samp, function(y) unlist(sapply(y, function(x) which(data$id %in% x), simplify = FALSE)))
    samp <- sapply(samp, function(x) as.integer(tabulate(x, nbins = n))) # n x ntree
  } else if (bootstrap == "none"){
    samp <- matrix(1, nrow = n, ncol = ntree)
  } else if (bootstrap == "by.user") {
    if (is.null(samp)) {
      stop("samp must not be NULL when bootstrapping by user\n")
    }
    if (is.matrix(samp)){
      if (!is.matrix(samp)) stop("samp must be a matrx\n")
      if (any(!is.finite(samp))) stop("samp must be finite\n")
      if (any(samp < 0)) stop("samp must be non-negative\n")
      if (all(dim(samp) != c(n, ntree))) stop("dimension of samp must be n x ntree\n")
      samp <- as.matrix(samp)  # transform into matrix
    }
  } else if (bootstrap == "by.root"){
    samp <- rep(1, n)
  } else {
    stop("Wrong bootstrap is given!\n ")
  }

  if (is.null(mtry)){
    mtry <- tune.ltrccif(formula = formula, data = data, id = id,
                         control = control, ntreeTry = ntree,
                         bootstrap = "by.user",
                         samptype = samptype,
                         sampfrac = sampfrac,
                         samp = samp,
                         na.action = na.action,
                         stepFactor = stepFactor,
                         applyfun = applyfun,
                         cores = cores,
                         trace = trace)
    print(sprintf("mtry is tuned to be %1.0f", mtry))
  }

  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (all(is.na(weights)) == 1) weights <- rep(1, NROW(y))
    s <- .logrank_trafo2(y[weights > 0, , drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  ret <- partykit::cforest(formula, data,
                           weights = samp,
                           perturb = perturb,
                           ytrafo = h2,
                           control = control,
                           na.action = na.action,
                           mtry = mtry,
                           ntree = ntree,
                           applyfun = applyfun,
                           cores = cores)
  ltrc_cif_obj<-ret
  
  ret$formulaLTRC <- formula
  ret$info$call <- Call
  ret$info$bootstrap <- bootstrap
  ret$info$samptype <- samptype
  ret$info$sampfrac <- sampfrac
  if (na.action == "na.omit"){
    ret$data$id <- data$id[complete.cases(data) == 1]
  } else {
    ret$data$id <- data$id
  }
  class(ret) <- "ltrccif"
  # revise here because we need individual trees
  return(list(ret=ret, ltrc_cif_obj=ltrc_cif_obj))
  
} # end function 

#-----------------------

tune.ltrccif  <- function(formula, data, id,
                          mtryStart = NULL, stepFactor = 2,
                          time.eval = NULL, time.tau = NULL,
                          ntreeTry = 100L,
                          bootstrap = c("by.sub", "by.root", "none", "by.user"),
                          samptype = c("swor","swr"),
                          sampfrac = 0.632,
                          samp = NULL,
                          na.action = "na.omit",
                          trace = TRUE,
                          doBest = FALSE,
                          plot = FALSE,
                          applyfun = NULL, cores = NULL,
                          control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                            mincriterion = 0, saveinfo = FALSE,
                                                            minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                            minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                                            minprob = 0.01)) {
    
    Call <- match.call()
    
    indx <- match(c('formula', 'id'), names(Call), nomatch = 0)
    if (indx[1] == 0) stop("a formula argument is required")
    
    # pull y-variable names
    yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
    yvar.names <- yvar.names[-length(yvar.names)]
    
    if (length(yvar.names) == 4){
        yvar.names = yvar.names[2:4]
    }
    n <- nrow(data)
    
    ## if not specified, the first one will be used as default
    bootstrap <- match.arg(bootstrap)
    samptype <- match.arg(samptype)
    
    # right-censored time from all observations
    Rtimes <- data[, yvar.names[2]]
    
    # extract the x-variable names
    xvar.names <- attr(terms(formula), 'term.labels')
    nvar <- length(xvar.names)
    
    if (is.null(mtryStart)){
        mtryStart <- ceiling(sqrt(nvar))
    }
    
    # this is a must, otherwise id cannot be passed to the next level in tune.ltrccif
    if (indx[2] == 0){
        ## If id is not present, then we add one more variable
        # mf$`(id)` <- 1:nrow(mf) ## No relabel, due to missing value problem do not need this for tuning output
        data$id <- 1:n # this is a must, otherwise id cannot be passed to the next level
    } else {
        ## If id is present, then we rename the column to be id
        names(data)[names(data) == deparse(substitute(id))] <- "id" # this is a must, otherwise id cannot be passed to the next level
    }
    
    data <- data[, c("id", yvar.names, xvar.names)]
    
    if (na.action == "na.omit") {
        takeid = which(complete.cases(data) == 1)
    } else if (na.action == "na.pass") {
        takeid = 1:n
    } else {
        stop("na.action can only be either 'na.omit' or 'na.pass'.")
    }
    
    id.sub <- unique(data$id[takeid])
    n.seu <- length(takeid)
    ## number of subjects
    n.sub <- length(id.sub)
    
    Rtimes <- Rtimes[takeid]
    ## This is to determine time.tau and time.eval, so is different from ltrccif.R
    if (n.seu == n.sub){ # time-invariant LTRC data
        # it includes the case 1) when id = NULL, which is that id.seu is not specified
        #                      2) when id is specified, but indeed LTRC time-invariant
        if (is.null(time.eval)){
            # estimated survival probabilities will be calculated at (a subset of) time.eval
            time.eval <- c(0, sort(unique(Rtimes)))
        }
    } else { # time-varying subject data
        if (is.null(time.eval)){
            # estimated survival probabilities will be calculated at (a subset of) time.eval
            time.eval <- c(0, sort(unique(Rtimes)), seq(max(Rtimes), 1.5 * max(Rtimes), length.out = 50)[-1])
        }
        if (is.null(time.tau)){
            # For i-th data, estimated survival probabilities only calculated up time.tau[i]
            time.tau <- sapply(1:n.sub, function(ii){
                1.5 * max(Rtimes[data$id[takeid] == id.sub[ii]])
            })
        }
    }
    
    # integrated Brier score of out-of-bag samples for a mtry value at test
    errorOOB_mtry <- function(eformula, edata, id,
                              emtryTest,
                              etpnt, etau,
                              entreeTry, econtrol,
                              ebootstrap,
                              esamptype,
                              esampfrac,
                              esamp,
                              ena.action, eapplyfun, ecores){
        cfOOB <- ltrccif(formula = eformula, data = edata, id = id,
                         mtry = emtryTest,
                         ntree = entreeTry,
                         control = econtrol,
                         bootstrap = ebootstrap,
                         samptype = esamptype,
                         sampfrac = esampfrac,
                         samp = esamp,
                         na.action = ena.action,
                         applyfun = eapplyfun,
                         cores = ecores)
        predOOB <- predictProb(object = cfOOB, time.eval = etpnt, time.tau = etau, OOB = TRUE)
        errorOOB <- sbrier_ltrc(obj = predOOB$survival.obj, id = predOOB$survival.id,
                                pred = predOOB, type = "IBS")
        rm(cfOOB)
        rm(predOOB)
        gc()
        return(errorOOB)
    }
    
    # # errorOld
    errorOld <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                              emtryTest = mtryStart,
                              etpnt = time.eval,
                              etau = time.tau,
                              entreeTry = ntreeTry,
                              econtrol = control,
                              ebootstrap = bootstrap,
                              esamptype = samptype,
                              esampfrac = sampfrac,
                              esamp = samp,
                              ena.action = na.action,
                              eapplyfun = applyfun,
                              ecores = cores)
    if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
    if (trace) {
        cat("mtry = ", mtryStart, " OOB Brier score = ",
            errorOld, "\n")
    }
    
    oobError <- list()
    oobError[[1]] <- errorOld
    names(oobError)[1] <- mtryStart
    
    for (direction in c("left", "right")) {
        if (trace) cat("Searching", direction, "...\n")
        mtryCur <- mtryStart
        while (mtryCur != nvar) {
            mtryOld <- mtryCur
            mtryCur <- if (direction == "left") {
                max(1, ceiling(mtryCur / stepFactor))
            } else {
                min(nvar, floor(mtryCur * stepFactor))
            }
            if (mtryCur == mtryOld) break
            
            errorCur <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                                      emtryTest = mtryCur,
                                      etpnt = time.eval,
                                      etau = time.tau,
                                      entreeTry = ntreeTry,
                                      econtrol = control,
                                      ebootstrap = bootstrap,
                                      esamptype = samptype,
                                      esampfrac = sampfrac,
                                      esamp = samp,
                                      ena.action = na.action,
                                      eapplyfun = applyfun,
                                      ecores = cores)
            
            if (trace) {
                cat("mtry = ", mtryCur, "\tOOB error = ", errorCur, "\n")
            }
            oobError[[as.character(mtryCur)]] <- errorCur
            errorOld <- errorCur
        }
    }
    mtry <- sort(as.numeric(names(oobError)))
    res_all <- unlist(oobError[as.character(mtry)])
    res_all <- cbind(mtry = mtry, OOBError = res_all)
    res <- res_all[which.min(res_all[, 2]), 1]
    
    if (plot) {
        res <- res_all
        plot(res_all, xlab = expression(m[try]), ylab = "OOB Error", type = "o", log = "x", xaxt = "n")
        axis(1, at=res_all[, "mtry"])
    }
    
    if (doBest)
        res <- ltrccif(formula = formula, data = data, id = id,
                       mtry = res, ntree = ntreeTry,
                       control = control,
                       bootstrap = bootstrap,
                       samptype = samptype,
                       sampfrac = sampfrac,
                       samp = samp,
                       na.action = na.action,
                       applyfun = applyfun,
                       cores = cores)
    
    return(res)
    
} # end function 

