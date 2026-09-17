# This file is a modification of TimeVaryingData_LTRCforests/analysis/utils/Loss_funct_tvary.R

shat_cox <- function(data, pred = NULL, tpnt = NULL, obj.roc = NULL){
    if (is.null(tpnt)){
        tpnt = c(0, sort(unique(data$Stop)))
    }
    N = length(unique(data$ID))
    Shatt = sapply(1:N, function(Ni) shat_funct_cox(Ni = Ni, data = data, pred = pred, tpnt = tpnt, obj.roc = obj.roc))
    return(Shatt)
} # end function 

#-------------------------------
## This function is for Cox and TSF
# if ("survfit.cox" %in% class(pred) || "survfitcox" %in% class(pred) || "pcoxsurvfit" %in% class(pred))
shat_funct_cox <- function(Ni, data, pred = NULL, tpnt, obj.roc = NULL){
## This function is to compute the estimated survival probability of the Ni-th subject

    id_uniq <- unique(data$ID)
    
    ## the i-th data
    TestData <- data[data$ID == id_uniq[Ni], ]
    
    TestT <- c(TestData[1, "Start"], TestData[, "Stop"])
    TestTIntN <- nrow(TestData)
    
    tpntL <- c(TestT, tpnt)
    torder <- order(tpntL)
    tpntLod <- tpntL[torder]
    tlen <- length(tpntLod)
    
    ## Compute the estimated survival probability of the Ni-th subject
    Shat_temp <- matrix(0, nrow = 1, ncol = tlen)
    
    r.ID <- findInterval(tpntLod, TestT)
    r.ID[r.ID > TestTIntN] <- TestTIntN
    
    jall <- unique(r.ID[r.ID > 0])
    nj <- length(jall)
    
    ## Deal with left-truncation
    Shat_temp[1, r.ID == 0] <- 1
    
    if(nj == 1){
        ## Get the index of the Pred to compute Shat
        II = which(data$ID == id_uniq[Ni])[jall[nj]]
        Shat_i = getsurv(pred[II], tpntLod[r.ID == jall[nj]])
        Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
    } else {
        ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
        ShatR_temp[1, 1] = 1
        # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
        qL = rep(0, nj)
        for (j in 1:nj){
            ## Get the index of the Pred to compute Shat
            II <- which(data$ID == id_uniq[Ni])[1] + jall[j] - 1
            Shat_j = getsurv(pred[II], tpntLod[r.ID == jall[j]])
            
            qL[j] <- Shat_j[1]
            # S_{j}(R_{j}), j=1,...nj-1
            jR = getsurv(pred[II], TestT[j + 1])
            ShatR_temp[1, j + 1] = jR / qL[j]
            Shat_temp[1, r.ID == jall[j]] <- Shat_j / qL[j]
        }
        
        ql0 <- which(qL == 0)
        if (length(ql0) > 0){
            if (any(qL > 0)){
                maxqlnot0 <- max(which(qL > 0))
                
                ql0lmax <- ql0[ql0 < maxqlnot0]
                ql0mmax <- ql0[ql0 >= maxqlnot0]
                ShatR_temp[1, ql0lmax + 1] <- 1
                Shat_temp[1, r.ID %in% jall[ql0lmax]] <- 1
                ShatR_temp[1, ql0mmax + 1] <- 0
                Shat_temp[1, r.ID %in% jall[ql0mmax]] <- 0
            } else {
                ShatR_b[1, 2:(nj + 1)] <- 0
                Shat_temp[1, r.ID %in% jall] <- 0
            }
        }
        m <- cumprod(ShatR_temp[1, 1:nj])
        for (j in 1:nj){
            Shat_temp[1, r.ID == jall[j]] <- Shat_temp[1, r.ID == jall[j]] * m[j]
        }
    }
    Shat <- Shat_temp[1, -match(TestT, tpntLod)]
    Shat

} # end function

#---------------------------------

shat_pcox <- function(data, pred = NULL, tpnt = NULL, obj.roc = NULL){
    if (is.null(tpnt)){
        tpnt = c(0, sort(unique(data$Stop)))
    }
    N = length(unique(data$ID))
    Shatt = sapply(1:N, function(Ni) shat_funct_pcox(Ni = Ni, data = data, pred = pred, tpnt = tpnt, obj.roc = obj.roc))
    return(Shatt)
} # end function
#-------------------------------

shat_funct_pcox <- function(Ni, data, pred = NULL, tpnt, obj.roc = NULL){
  ## This function is to compute the estimated survival probability of the Ni-th subject
  id_uniq <- unique(data$ID)
  
  ## the i-th data
  TestData <- data[data$ID == id_uniq[Ni], ]
  
  TestT <- c(TestData[1, "Start"], TestData[, "Stop"])
  TestTIntN <- nrow(TestData)
  
  tpntL <- c(TestT, tpnt)
  torder <- order(tpntL)
  tpntLod <- tpntL[torder]
  tlen <- length(tpntLod)
  
  ## Compute the estimated survival probability of the Ni-th subject
  Shat_temp <- matrix(0, nrow = 1, ncol = tlen)
  
  r.ID <- findInterval(tpntLod, TestT)
  r.ID[r.ID > TestTIntN] <- TestTIntN
  
  jall <- unique(r.ID[r.ID > 0])
  nj <- length(jall)
  
  ## Deal with left-truncation
  Shat_temp[1, r.ID == 0] <- 1
  
  if(nj == 1){
    ## Get the index of the Pred to compute Shat
    II = which(data$ID == id_uniq[Ni])[jall[nj]]
    Shat_i = getsurv_pcox(pred$surv[, II], pred$time, tpntLod[r.ID == jall[nj]])
    Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
  } else {
    ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
    ShatR_temp[1, 1] = 1
    # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
    qL = rep(0, nj)
    for (j in 1:nj){
      ## Get the index of the Pred to compute Shat
      II <- which(data$ID == id_uniq[Ni])[1] + jall[j] - 1
      Shat_j = getsurv_pcox(pred$surv[, II], pred$time, tpntLod[r.ID == jall[j]])
      
      qL[j] <- Shat_j[1]
      # S_{j}(R_{j}), j=1,...nj-1
      jR = getsurv_pcox(pred$surv[, II], pred$time, TestT[j + 1])
      ShatR_temp[1, j + 1] = jR / qL[j]
      Shat_temp[1, r.ID == jall[j]] <- Shat_j / qL[j]
    }
    
    ql0 <- which(qL == 0)
    if (length(ql0) > 0){
      if (any(qL > 0)){
        maxqlnot0 <- max(which(qL > 0))
        
        ql0lmax <- ql0[ql0 < maxqlnot0]
        ql0mmax <- ql0[ql0 >= maxqlnot0]
        ShatR_temp[1, ql0lmax + 1] <- 1
        Shat_temp[1, r.ID %in% jall[ql0lmax]] <- 1
        ShatR_temp[1, ql0mmax + 1] <- 0
        Shat_temp[1, r.ID %in% jall[ql0mmax]] <- 0
      } else {
        ShatR_b[1, 2:(nj + 1)] <- 0
        Shat_temp[1, r.ID %in% jall] <- 0
      }
    }
    m <- cumprod(ShatR_temp[1, 1:nj])
    for (j in 1:nj){
      Shat_temp[1, r.ID == jall[j]] <- Shat_temp[1, r.ID == jall[j]] * m[j]
    }
  }
  Shat <- Shat_temp[1, -match(TestT, tpntLod)]
  Shat
  
} # end function

