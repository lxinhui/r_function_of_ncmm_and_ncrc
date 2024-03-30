#Function3----------------------------------------------------#
NCcorrected <- function(mismeasuredvariable,precisievarable,dataused,datavalidate,datatype,aim,boottime){
  #'@param mismeasuredvariable vector of variables measured with error, e.g., precisievarable = c("X","C").
  #'@param precisievarable vector of variables without error, e.g., precisievarable = c("Y","M").
  #'@param dataused the dataframe of main dataset, e.g, data = data2, data2 <- data.frame(X1,C1,Y,M).
  #'@param datavalidate if datatype = "replicate", then datavalidation = NA; otherwise, datavalidation is the dataframe of internal validation data, e.g., datavalidate = data3, data3 <- data.frame(vX,vC,X2,C2).
  #'@param datatype "replicate" or "validation".
  #'@param aim correction, detection, or reduction.
  #'@param boottime: times for bootstrap.
  
  if(datatype=="replicate"){
    estimates <- NCreplicate(mismeasuredvar = mismeasuredvariable,precisievar = precisievarable,data = dataused,aim = aim)
    estimates_NCRC <- estimates[1,2]
    estimates_NCMM <- estimates[1,1]
    
    ###estimate SEs through bootstrap###
    estimateboot_NCRC <- rep(NA,boottime)
    estimateboot_NCMM <- rep(NA,boottime)
    for(TT in 1:boottime){
      cat(TT)
      samp <- sample(1:nrow(dataused),nrow(dataused),replace=T)
      data_boot <- dataused[samp,]
      rownames(data_boot) <- 1:nrow(data_boot)
      estimates_boot <- NCreplicate(mismeasuredvar = mismeasuredvariable,precisievar = precisievarable,data = data_boot,aim = aim)
      estimateboot_NCRC[TT] <- estimates_boot[1,2]
      estimateboot_NCMM[TT] <- estimates_boot[1,1]
    }
  }
  
  if(datatype=="validation"){
    estimates <- NCvalidation(mismeasuredvar = mismeasuredvariable,precisievar = precisievarable,datamain = dataused,datavalidation = datavalidate,aim = aim)
    estimates_NCRC <- estimates[1,2]
    estimates_NCMM <- estimates[1,1]
    
    ###estimate variances through bootstrap###
    estimateboot_NCRC <- rep(NA,boottime)
    estimateboot_NCMM <- rep(NA,boottime)
    for(TT in 1:boottime){
      cat(TT)
      samp <- sample(1:nrow(dataused),nrow(dataused),replace=T)
      data_boot <- dataused[samp,]
      rownames(data_boot) <- 1:nrow(data_boot)
      sampv <- sample(1:nrow(datavalidate),nrow(datavalidate),replace=T)
      data_bootv <- datavalidate[sampv,]
      rownames(data_bootv) <- 1:nrow(data_bootv)
      estimates_boot <- NCvalidation(mismeasuredvar = mismeasuredvariable,precisievar = precisievarable,datamain = data_boot,datavalidation = data_bootv,aim = aim)
      estimateboot_NCRC[TT] <- estimates_boot[1,2]
      estimateboot_NCMM[TT] <- estimates_boot[1,1]
    }
  }
  
  SE_NCRC <- sqrt((sum((estimateboot_NCRC-estimates_NCRC)^2))/(boottime-1))
  SE_NCMM <- sqrt((sum((estimateboot_NCMM-estimates_NCMM)^2))/(boottime-1))
  
  CILOW_RC <- estimates_NCRC - 1.96*(SE_NCRC)
  CIUP_RC <- estimates_NCRC + 1.96*(SE_NCRC)
  
  CILOW_MM <- estimates_NCMM - 1.96*(SE_NCMM)
  CIUP_MM <- estimates_NCMM + 1.96*(SE_NCMM)
  
  results <- data.frame(matrix(c(estimates_NCRC,estimates_NCMM,SE_NCRC,SE_NCMM,CILOW_RC,CILOW_MM,CIUP_RC,CIUP_MM),2,4))
  colnames(results) <- c("Estimate","SE","95%CI_lower","95%CI_upper")
  rownames(results) <- c("NC-RC","NC-MM")
  return(results)
}

# Examples:
#   1# Both exposure and negative control exposure are measured with error, perform unmeasured confounding correction using replication subset when additive equi-confounding assumption hold:
# dataA <- data.frame(X1,X2,C1,C2,Y,M)
# NCcorrected(mismeasuredvariable=c("X","C"),precisievarable=c("Y","M"),dataused=dataA,datavalidate=NA,datatype="replicate",aim="correction",boottime=100)
# 
# 2# Both exposure and negative control exposure are measured with error, perform unmeasured confounding detection using validation dataset:
# dataA <- data.frame(X1,C1,Y,M)
# dataB <- data.frame(vX,vC,X2,C2)
# NCcorrected(mismeasuredvariable=c("X","C"),precisievarable=c("Y","M"),dataused=dataA,datavalidate=dataB,datatype="validation",aim="detection",boottime=100)
