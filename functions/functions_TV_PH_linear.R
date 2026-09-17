
T_linr2 <- function(TALL, X, U, Beta){  
    
    u = U
    x1 = X$X1
    x2 = X$X2
    x3 = X$X3
    x4 = X$X4
    
    tlen = length(TALL)
    R0 = exp(Beta[1]*x1 + Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5])
    
    Lambda = 0.001 * exp(-sum(Beta[1:4]) * 0.5) 
    
    V = 2
    Alpha = 0
    
    R = Lambda*R0[-tlen] * (TALL[-1]^V - TALL[-tlen]^V)
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID] + TALL[R.ID]^V)^(1/V)
    
    result = list(Time = TT, Row = R.ID, Lambda = Lambda, Beta = Beta, Alpha = Alpha, V = V, Xi = R0)
    return(result)
    
} # T_linr

#----------------------------------------------------


simple_PH_linear2 <- function(N = 100, Beta) {
    
    npseu = 10
    # Pre-allocate data frame
    Data <- as.data.frame(matrix(NA, npseu * N, 15))
    names(Data) <- c("I", "ID", "X1", "X2", "X3", "X4", "X5", "X6", 
                     "X7", "X8", "Start", "Stop", "C", "Event", "Xi")
    Data$ID <- rep(1:N, each = npseu)
    
    ##  Generate Covariates -----------------------------------------
    # Time-invariant
    Data$X1 <- rep(runif(N, 0, 1), each = npseu)
    Data$X5 <- rep(runif(N, 0, 1), each = npseu)
    
    # Time-varying
    Data$X2 <- runif(npseu * N, 0, 1)
    Data$X6 <- runif(npseu * N, 0, 1)
    
    ## Generate Event Times ----------------------------------------
    Count = 1
    while(Count <= N) {
        TS <- rep(0, npseu - 1)
        
        while(any(diff(sort(TS)) < 0.4)) {
            TS <- sort(rtrunc(npseu-1, spec="beta", a=0.001, b=1, shape1=0.05, shape2=5)) * 120
        }
        
        u = runif(1)
        k = runif(2); k_1 = runif(2)
        
        # --- CRITICAL FIX 2: BASELINE VARIANCE ---
        # Generate categorical time-varying X3, X7
        # Sampling with replacement and sorting creates a perfect step function,
        # but ensures they don't all start strictly at 0.
        Data[Data$ID == Count, ]$X3 <- sort(sample(0:2, npseu, replace = TRUE))
        Data[Data$ID == Count, ]$X7 <- sort(sample(0:2, npseu, replace = TRUE))
        
        # Generate continuous time-varying X4, X8
        
        Data[Data$ID == Count, ]$X4 <- k[1] * c(0, TS/120) + k[2]
        Data[Data$ID == Count, ]$X8 <- k_1[1] * c(0, TS/120) + k_1[2]
        
        Data[Data$ID == Count, ]$Start <- c(0, TS)
        Data[Data$ID == Count, ]$Stop <- c(TS, NA)
        
        # Calculate True Survival Time based on PH model
        RT <- T_linr2(TALL = c(0, TS), X = Data[Data$ID == Count, ], U = u, Beta)
        
        t = RT$Time
        rID = RT$Row
        Data[Data$ID == Count, ]$Xi = RT$Xi
        
        # Apply Event to the correct interval
        if(rID == 1) {
            Data[Data$ID == Count, ][1, ]$Stop = t
            Data[Data$ID == Count, ][1, ]$Event = 1
        } else {
            Data[Data$ID == Count, ][1:(rID - 1), ]$Event = 0
            Data[Data$ID == Count, ][rID, ]$Event = 1
            Data[Data$ID == Count, ][rID, ]$Stop = t
        }
        Count = Count + 1
    } # end while 
    
    # Remove rows beyond the event time
    DATA <- Data[!is.na(Data$Event), ]
    rm(Data); gc()
    
    # Add Censoring 
    DATA$C <- 0
    censor_rate <- 0.2
    is_forced_censor <- rbinom(N, 1, censor_rate) 
    Censor.time = rexp(N, rate = 1/30) + 0.05 # Add buffer to avoid 0 time
    
    for(j in 1:N) {
        subj_rows <- which(DATA$ID == j)
        if(length(subj_rows) == 0) next
        
        # The time the event WAS supposed to happen
        actual_event_time <- DATA$Stop[max(subj_rows)]
        
        # Censor if forced OR if the administrative censoring time is earlier than event
        if(is_forced_censor[j] == 1 || actual_event_time > Censor.time[j]) {
            
            # The new observed time is the earlier of the two
            end_time <- min(actual_event_time, Censor.time[j])
            
            # Determine which interval the censoring falls into
            Vec <- c(0, DATA[subj_rows, ]$Stop)
            ID_in_subj <- findInterval(end_time, Vec)
            
            # Update the target row and remove any rows that were supposed to come after
            safe_ID <- max(1, min(ID_in_subj, length(subj_rows)))
            target_idx <- subj_rows[safe_ID]
            
            DATA[target_idx, ]$C <- 1
            DATA[target_idx, ]$Event <- 0
            DATA[target_idx, ]$Stop <- end_time
            
            if(target_idx < max(subj_rows)) {
                DATA[(target_idx + 1):max(subj_rows), ]$Event <- NA
            } # end target_idx
        } # end if 
    } # end j 
    
    # Final cleaning
    Data <- DATA[!is.na(DATA$Event), ]
    
    ## Finalize Outputs 
    RET <- list()
    Data$I <- 1:nrow(Data)
    
    # Clean up tiny precision errors for Cox models
    Data$Start <- round(Data$Start, 3)
    Data$Stop <- round(Data$Stop, 3)
    Data$Stop[Data$Start == Data$Stop] <- Data$Stop[Data$Start == Data$Stop] + 0.001
    
    RET$fullData <- Data
    RET$Info = list(Coeff = list(Lambda = RT$Lambda, Alpha = RT$Alpha, Beta = RT$Beta, V = RT$V), Set = "PH")
    
    # Create Baseline/Static Dataset (One row per ID)
    DATA_BASE = data.frame(matrix(0, nrow = N, ncol = ncol(Data)))
    names(DATA_BASE) = names(Data)
    for (ii in 1:N) {
        subj_subset <- Data[Data$ID == ii, ]
        ni = nrow(subj_subset)
        # Take first row for covariates, but update Stop and Event to the final state
        DATA_BASE[ii, ] = subj_subset[1, ]
        DATA_BASE[ii, ]$Stop = subj_subset$Stop[ni]
        DATA_BASE[ii, ]$Event = subj_subset$Event[ni]
    }
    
    DATA_BASE$I = 1:nrow(DATA_BASE)
    RET$baselineData = DATA_BASE
    
    return(RET)
} # end function 