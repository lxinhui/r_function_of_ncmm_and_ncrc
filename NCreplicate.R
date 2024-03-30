#Function1----------------------------------------------------#
library(MASS)
NCreplicate <- function(mismeasuredvar,precisievar,data,aim){
  #'@param mismeasuredvar vector of the names of error-prone variables. e.g., mismeasuredvar=c("X","C"), the order should be the same as that in data.
  #'@param precisievar vector of the names of preciously measured variables (must include Y). e.g., precisievar=c("Y","M")
  #'@param data dataframe of main dataset including all observed variables. e.g., data <- data1, data1 <- data.frame(X1,X2,C1,C2,Y,M)
  #'@param aim "correction", "detection", or "reduction"
  
  NMIS <- length(mismeasuredvar)
  
  for(ii in 1:NMIS){
    if(ii==1){
      linshi <- data[,grep(mismeasuredvar[1],colnames(data))]
      data$nrep <- NA
      data$linshi <- NA
      for(jj in 1:nrow(linshi)){
        data[jj,which(colnames(data)=="nrep")] <- nrow(na.omit(t(linshi[jj,])))
        data[jj,which(colnames(data)=="linshi")] <- (1/data[jj,which(colnames(data)=="nrep")])*sum(na.omit(t(linshi[jj,])))
      }
      colnames(data)[which(colnames(data)=="linshi")] <- paste0(mismeasuredvar[ii],"_AVE")
    }
    if(ii!=1){
      linshi <- data[,grep(mismeasuredvar[ii],colnames(data))]
      data$linshi <- NA
      for(jj in 1:nrow(linshi)){
        data[jj,which(colnames(data)=="linshi")] <- (1/data[jj,which(colnames(data)=="nrep")])*sum(na.omit(t(linshi[jj,])))
      }
      colnames(data)[which(colnames(data)=="linshi")] <- paste0(mismeasuredvar[ii],"_AVE")
    }
  }
  
  A_AVE <- as.matrix(data[,grep("AVE",colnames(data))])
  
  if(length(precisievar)!=1){
    ###calculate parameters###
    deta_ME <- matrix(0,NMIS,NMIS)
    
    for(ii in rownames(data[data$nrep>1,])){
      for(jj in 1:data[as.numeric(ii),which(colnames(data)=="nrep")]){
        Ai <- matrix(as.numeric(data[as.numeric(ii),which(colnames(data)%in%c(paste0(mismeasuredvar,jj)))]))
        A_AVEi <- as.matrix(A_AVE[as.numeric(ii),])
        linshi <- (Ai-A_AVEi)%*%(t(Ai-A_AVEi))
        deta_ME <- deta_ME+linshi
      }
    }
    deta_ME <- deta_ME*(1/sum(data$nrep-1))
    ###
    MM <- data$nrep
    ###
    uA <- matrix(0,NMIS,1)
    for(ii in 1:nrow(A_AVE)){
      linshi <- MM[ii]*as.matrix(A_AVE[ii,])
      uA <- uA+linshi
    }
    uA <- uA/sum(MM)
    ###
    B <- data.frame(matrix(NA,nrow(data),0))
    for(ii in precisievar[-which(precisievar=="Y")]){
      linshi <- as.matrix(data[,grep(ii,colnames(data))])
      colnames(linshi) <- ii
      B <- cbind(B,linshi)
    }
    uB <- apply(B,2,mean)
    ###
    v <- (sum(MM))-(sum(MM^2))/(sum(MM))
    ###
    detaA <- matrix(0,NMIS,NMIS)
    for(ii in 1:nrow(data)){
      A_AVEi <- as.matrix(A_AVE[ii,])
      linshi <- MM[ii]*(as.matrix(A_AVEi-uA)%*%(t(A_AVEi-uA)))
      detaA <- detaA+linshi
    }
    detaA <- (detaA-(nrow(data)-1)*(deta_ME))/v
    ###
    detaAB <- matrix(0,NMIS,length(precisievar)-1)
    for(ii in 1:nrow(data)){
      A_AVEi <- as.matrix(A_AVE[ii,])
      Bi <- as.matrix(B[ii,])
      linshi <- MM[ii]*((A_AVEi-uA)%*%(t(Bi-uB)))
      detaAB <- detaAB+linshi
    }
    detaAB <- detaAB/v
    ###
    detaB <- matrix(0,length(precisievar)-1,length(precisievar)-1)
    for(ii in 1:nrow(data)){
      Bi <- as.matrix(B[ii,])
      linshi <- ((Bi-uB)%*%(t(Bi-uB)))
      detaB <- detaB+linshi
    }
    detaB <- detaB/(nrow(data)-1)
    
    #####NC-MM#####
    detaAY <- matrix(0,NMIS,1)
    uY <- mean(data$Y) 
    for(ii in 1:nrow(data)){
      A_AVEi <- as.matrix(A_AVE[ii,])
      linshi <- MM[ii]*((A_AVEi-uA)*(data$Y[ii]-uY))
      detaAY <- detaAY+linshi
    }
    detaAY <- detaAY/v
    
    detaBY <- matrix(0,length(precisievar)-1,1)
    for(ii in 1:nrow(data)){
      Bi <- as.matrix(B[ii,])
      linshi <- ((Bi-uB)*(data$Y[ii]-uY))
      detaBY <- detaBY+linshi
    }
    detaBY <- detaBY/(nrow(data)-1)
    
    detaAY <- matrix(rbind(as.matrix(detaAY),as.matrix(detaBY)))
    
    beta <- ginv(as.matrix(rbind(cbind(as.matrix(detaA),as.matrix(detaAB)),
                                 cbind(t(as.matrix(detaAB)),as.matrix(detaB)))))%*%detaAY
    
    betac_AD <- beta[grep("C",c(mismeasuredvar,colnames(B))),1]
    betax_AD <- beta[grep("X",c(mismeasuredvar,colnames(B))),1]
    
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
    for(ii in 1:nrow(data)){
      vect1 <- as.matrix(cbind(as.matrix(detaA),as.matrix(detaAB)))
      vect2 <- ginv(as.matrix(rbind(cbind(as.matrix(detaA+(deta_ME/MM[ii])),as.matrix(detaAB)),cbind(t(as.matrix(detaAB)),as.matrix(detaB)))))
      vect3 <- as.matrix(rbind(as.matrix(A_AVE[ii,]-uA),as.matrix(B[ii,]-uB)))
      
      linshi <- uA+(vect1%*%vect2)%*%vect3
      EC <- rbind(EC,t(linshi))
    }
    
    colnames(EC) <- mismeasuredvar
    O <- data.frame(data$Y)
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
    ###calculate parameters###
    deta_ME <- matrix(0,NMIS,NMIS)
    
    for(ii in rownames(data[data$nrep>1,])){
      for(jj in 1:data[as.numeric(ii),which(colnames(data)=="nrep")]){
        Ai <- matrix(as.numeric(data[as.numeric(ii),which(colnames(data)%in%c(paste0(mismeasuredvar,jj)))]))
        A_AVEi <- as.matrix(A_AVE[as.numeric(ii),])
        linshi <- (Ai-A_AVEi)%*%(t(Ai-A_AVEi))
        deta_ME <- deta_ME+linshi
      }
    }
    
    deta_ME <- deta_ME*(1/sum(data$nrep-1))
    ###
    MM <- data$nrep
    ###
    uA <- matrix(0,NMIS,1)
    for(ii in 1:nrow(A_AVE)){
      linshi <- MM[ii]*as.matrix(A_AVE[ii,])
      uA <- uA+linshi
    }
    uA <- uA/sum(MM)
    ###
    v <- (sum(MM))-(sum(MM^2))/(sum(MM))
    ###
    detaA <- matrix(0,NMIS,NMIS)
    for(ii in 1:nrow(data)){
      A_AVEi <- as.matrix(A_AVE[ii,])
      linshi <- MM[ii]*(as.matrix(A_AVEi-uA)%*%(t(A_AVEi-uA)))
      detaA <- detaA+linshi
    }
    detaA <- (detaA-(nrow(data)-1)*(deta_ME))/v
    #####NC-MM#####
    detaAY <- matrix(0,NMIS,1)
    uY <- mean(data$Y) 
    for(ii in 1:nrow(data)){
      A_AVEi <- as.matrix(A_AVE[ii,])
      linshi <- MM[ii]*((A_AVEi-uA)*(data$Y[ii]-uY))
      detaAY <- detaAY+linshi
    }
    detaAY <- detaAY/v
    
    beta <- ginv(as.matrix(detaA))%*%detaAY
    
    betac_AD <- beta[grep("C",c(mismeasuredvar)),1]
    betax_AD <- beta[grep("X",c(mismeasuredvar)),1]
    
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
    for(ii in 1:nrow(data)){
      vect1 <- as.matrix(detaA)
      vect2 <- ginv(as.matrix(detaA+(deta_ME/MM[ii])))
      vect3 <- as.matrix(as.matrix(A_AVE[ii,]-uA))
      
      linshi <- uA+(vect1%*%vect2)%*%vect3
      EC <- rbind(EC,t(linshi))
    }
    
    colnames(EC) <- mismeasuredvar
    O <- data.frame(data$Y)
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
