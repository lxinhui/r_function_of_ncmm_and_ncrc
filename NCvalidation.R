#Function2----------------------------------------------------#
NCvalidation <- function(mismeasuredvar,precisievar,datamain,datavalidation,aim){
  #'@param mismeasuredvar vector of the names of error-prone variables. e.g., mismeasuredvar=c("X","C"), the order should be the same as that in data.
  #'@param precisievar vector of the names of preciously measured variables (must include Y). e.g., precisievar=c("Y","M")
  #'@param datamain dataframe of main dataset including all observed variables. e.g., data = data2, data2 <- data.frame(X1,C1,Y,M)
  #'@param datavalidation dataframe of internal validation dataset including the true value and observed value of error-prone variables. e.g., data = data3, data3 <- data.frame(vX,vC,X2,C2)
  #'@param aim "correction", "detection", or "reduction"
  
  NMIS <- length(mismeasuredvar)
  
  ###calculate parameters###
  deta_ME <- matrix(0,NMIS,NMIS)
  
  for(ii in 1:NMIS){
    for(jj in 1:NMIS){
      if(ii==jj){
        deta_ME[ii,jj] <- var(datavalidation[,grep(mismeasuredvar[ii],colnames(datavalidation))[2]]-datavalidation[,grep(mismeasuredvar[ii],colnames(datavalidation))[1]])
      }
      if(ii!=jj){
        deta_ME[ii,jj] <- cov(datavalidation[,grep(mismeasuredvar[ii],colnames(datavalidation))[2]],datavalidation[,grep(mismeasuredvar[jj],colnames(datavalidation))[2]])-
          cov(datavalidation[,grep(mismeasuredvar[ii],colnames(datavalidation))[1]],datavalidation[,grep(mismeasuredvar[jj],colnames(datavalidation))[1]])
      }
    }
  }
  
  if(length(precisievar)!=1){
    ###
    A1 <- as.matrix(datamain[,grep(1,colnames(datamain))])
    uA <- as.matrix(apply(A1, 2, mean))
    ###
    B <- as.matrix(datamain[,-c(grep(1,colnames(datamain)),which(colnames(datamain)=="Y"))])
    uB <- as.matrix(apply(B,2,mean))
    ###
    detaA <- matrix(0,NMIS,NMIS)
    for(ii in 1:NMIS){
      for(jj in 1:NMIS){
        if(ii==jj){
          detaA[ii,jj] <- var(datamain[,grep(mismeasuredvar[ii],colnames(datamain))])
        }
        if(ii!=jj){
          detaA[ii,jj] <- cov(datamain[,grep(mismeasuredvar[ii],colnames(datamain))],datamain[,grep(mismeasuredvar[jj],colnames(datamain))])
        }
      }
    }  
    detaA <- detaA-deta_ME
    ###
    detaAB <- matrix(0,NMIS,length(precisievar)-1)
    for(ii in 1:NMIS){
      for(jj in 1:(length(precisievar)-1)){
        detaAB[ii,jj] <- cov(datamain[,grep(mismeasuredvar[ii],colnames(datamain))],datamain[,grep(precisievar[jj+1],colnames(datamain))])
      }
    }
    ###
    detaB <- matrix(0,length(precisievar)-1,length(precisievar)-1)
    for(ii in 1:length(precisievar)-1){
      for(jj in 1:(length(precisievar)-1)){
        detaB[ii,jj] <- cov(datamain[,grep(precisievar[ii+1],colnames(datamain))],datamain[,grep(precisievar[jj+1],colnames(datamain))])
      }
    }
    
    #####NC-MM#####
    detaAY <- matrix(0,NMIS,1)
    uY <- mean(datamain$Y) 
    for(ii in 1:NMIS){
      detaAY[ii,1] <- cov(datamain[,grep(mismeasuredvar[ii],colnames(datamain))],datamain[,which(colnames(datamain)=="Y")])
    }
    
    detaBY <- matrix(0,length(precisievar)-1,1)
    for(ii in 1:(length(precisievar)-1)){
      detaBY[ii,1] <- cov(datamain[,grep(precisievar[ii+1],colnames(datamain))],datamain[,which(colnames(datamain)=="Y")])
    }
    
    detaAY <- matrix(rbind(as.matrix(detaAY),as.matrix(detaBY)))
    
    beta <- ginv(as.matrix(rbind(cbind(as.matrix(detaA),as.matrix(detaAB)),
                                 cbind(t(as.matrix(detaAB)),as.matrix(detaB)))))%*%detaAY
    
    betac_AD <- beta[grep("C",c(colnames(A1),colnames(B))),1]
    betax_AD <- beta[grep("X",c(colnames(A1),colnames(B))),1]
    
    if(aim=="correction"){
      beta_AD <- betax_AD-betac_AD 
    }
    if(aim=="detection"){
      beta_AD <- betac_AD 
    }
    if(aim=="reduction"){
      beta_AD <- betax_AD 
    }
    
    #######NC-RC#########
    EC <- data.frame(matrix(NA,0,NMIS))
    for(ii in 1:nrow(datamain)){
      vect1 <- as.matrix(cbind(as.matrix(detaA),as.matrix(detaAB)))
      vect2 <- ginv(as.matrix(rbind(cbind(as.matrix(detaA+deta_ME),as.matrix(detaAB)),cbind(t(as.matrix(detaAB)),as.matrix(detaB)))))
      vect3 <- as.matrix(rbind(as.matrix(A1[ii,]-uA),as.matrix(B[ii]-uB)))
      
      linshi <- uA+(vect1%*%vect2)%*%vect3
      EC <- rbind(EC,t(linshi))
    }
    
    colnames(EC) <- mismeasuredvar
    O <- data.frame(datamain$Y)
    colnames(O) <- "Y"
    linshi <- data.frame(EC,B,O)
    func <- "Y~"
    for(ii in 1:(ncol(linshi)-2)){
      func <- paste0(func,colnames(linshi)[ii],"+")
    }
    func <- paste0(func,colnames(linshi)[(ncol(linshi)-1)])
    
    # regree <- paste0("Y~",colnames(EC)[1],"+",colnames(EC)[2],"+M")
    regression2 <- lm(func,data = linshi)
    betax_RC <- regression2$coefficients[grep("X",colnames(linshi))+1]
    betac_RC <- regression2$coefficients[grep("C",colnames(linshi))+1]
    
    if(aim=="correction"){
      beta_RC <- betax_RC-betac_RC 
    }
    if(aim=="detection"){
      beta_RC <- betac_RC 
    }
    if(aim=="reduction"){
      beta_RC <- betax_RC 
    }
  }
  
  if(length(precisievar)==1){
    ###
    A1 <- as.matrix(datamain[,grep(1,colnames(datamain))])
    uA <- as.matrix(apply(A1, 2, mean))
    ###
    detaA <- matrix(0,NMIS,NMIS)
    for(ii in 1:NMIS){
      for(jj in 1:NMIS){
        if(ii==jj){
          detaA[ii,jj] <- var(datamain[,grep(mismeasuredvar[ii],colnames(datamain))])
        }
        if(ii!=jj){
          detaA[ii,jj] <- cov(datamain[,grep(mismeasuredvar[ii],colnames(datamain))],datamain[,grep(mismeasuredvar[jj],colnames(datamain))])
        }
      }
    }  
    detaA <- detaA-deta_ME
    
    #####NC-MM#####
    detaAY <- matrix(0,NMIS,1)
    uY <- mean(datamain$Y) 
    for(ii in 1:NMIS){
      detaAY[ii,1] <- cov(datamain[,grep(mismeasuredvar[ii],colnames(datamain))],datamain[,which(colnames(datamain)=="Y")])
    }
    
    beta <- ginv(as.matrix(detaA))%*%detaAY
    
    betac_AD <- beta[grep("C",c(colnames(A1))),1]
    betax_AD <- beta[grep("X",c(colnames(A1))),1]
    
    if(aim=="correction"){
      beta_AD <- betax_AD-betac_AD 
    }
    if(aim=="detection"){
      beta_AD <- betac_AD 
    }
    if(aim=="reduction"){
      beta_AD <- betax_AD 
    }
    
    #######NC-RC#########
    EC <- data.frame(matrix(NA,0,NMIS))
    for(ii in 1:nrow(datamain)){
      vect1 <- as.matrix(cbind(as.matrix(detaA)))
      vect2 <- ginv(as.matrix(rbind(cbind(as.matrix(detaA+deta_ME)))))
      vect3 <- as.matrix(rbind(as.matrix(A1[ii,]-uA)))
      
      linshi <- uA+(vect1%*%vect2)%*%vect3
      EC <- rbind(EC,t(linshi))
    }
    
    colnames(EC) <- mismeasuredvar
    O <- data.frame(datamain$Y)
    colnames(O) <- "Y"
    linshi <- data.frame(EC,O)
    func <- "Y~"
    for(ii in 1:(ncol(linshi)-2)){
      func <- paste0(func,colnames(linshi)[ii],"+")
    }
    func <- paste0(func,colnames(linshi)[(ncol(linshi)-1)])
    
    # regree <- paste0("Y~",colnames(EC)[1],"+",colnames(EC)[2],"+M")
    regression2 <- lm(func,data = linshi)
    betax_RC <- regression2$coefficients[grep("X",colnames(linshi))+1]
    betac_RC <- regression2$coefficients[grep("C",colnames(linshi))+1]
    
    if(aim=="correction"){
      beta_RC <- betax_RC-betac_RC 
    }
    if(aim=="detection"){
      beta_RC <- betac_RC 
    }
    if(aim=="reduction"){
      beta_RC <- betax_RC 
    }
  }
  
  
  result <- data.frame(beta_AD,beta_RC)
  colnames(result) <- c("NC-MM estimate","NC-RC estimate")
  return(result)
}
