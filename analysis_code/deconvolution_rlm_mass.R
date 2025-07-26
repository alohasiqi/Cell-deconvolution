rlm_func <- function(mixData,sigMat){
  library(som)
  library(MASS)
  library(pheatmap)
  library(e1071)
  library(methods)
  mixDataSort <- merge(mixData,sigMat,by="GeneSymbol")
  mixDN <- (ncol(mixData)-1)
  sigMatN <- (ncol(sigMat)-1)
  
  mixData <- mixDataSort[,c(2:(mixDN+1)),drop=FALSE]
  rownames(mixData) <- mixDataSort[,1]
  
  sigMat <- mixDataSort[,c((mixDN+2):ncol(mixDataSort)),drop=FALSE]
  rownames(sigMat) <- mixDataSort[,1]
  
  #####normalized the mixData and sigMat to zero mean and unit variance to decrease running time
  ###For mixdata, each sample is normalized with mean 0 and sd 1
  ###For sigMat, all samples together are normalized with mean 0 and sd 1
  mixData <- as.data.frame(normalize(data.matrix(mixData),byrow=FALSE))
  sigMat <- (sigMat-mean(as.matrix(sigMat)))/sd(as.vector(as.matrix(sigMat)))
   
  dataSeq <- paste("data[,1] ~ ",sep="")
  for(i in c(2:(sigMatN+1))){
    dataSeq <- paste(dataSeq,"data[,",i,"]+",sep="")
  }
  dataSeq <- substr(dataSeq,0,nchar(dataSeq)-1) 
  
  resultMatrix <- array(0,dim=c(mixDN,sigMatN+3))
  rownames(resultMatrix) <- colnames(mixData)
  
  colnames(resultMatrix) <- c(colnames(sigMat),"RMSE","Cor","Pvalue")
  cat("Total number of samples: ",mixDN,"\n")
  for (j in 1:mixDN) {
  	cat("Sample: ",colnames(mixData)[j],"\n")
    data = cbind(mixData[,j], sigMat)
    sampleName = colnames(mixData)[j]
    tunedModel = rlm(eval(parse(text= dataSeq)), method = "MM")
    predicted_tunedModel <- predict(tunedModel, newdata = data) 
    
    errorT <- data[,1] - predicted_tunedModel 
        
    tunedModelRMSE <- sqrt(mean(errorT^2))
    
    tst <- cor.test(data[,1], predicted_tunedModel)
    corr1 <- round(tst$estimate, 3)
    pval1 <- round(tst$p.value,3)   
    
    coefTuned <- tunedModel$coefficients[2: length(tunedModel$coefficients)]
    coefTuned[coefTuned <0 ] <- 0
    #normalize the coefs
    coefTunedNormed <- coefTuned /sum(coefTuned)
    x <- c(coefTunedNormed, tunedModelRMSE, corr1, pval1)
    resultMatrix[j,] <- x
  }
  
  return(resultMatrix)
}
