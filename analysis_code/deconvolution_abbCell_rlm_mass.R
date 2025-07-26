library(MASS)
library(pheatmap)
library(e1071)
library(som)
library(methods)

abbCell_rlm <- function(method="rlm_func", mixDataDir=NULL, mixData=NULL, cellMatrixDir=NULL, dataDescription="", projectName=NULL, outputDirectory="./", annotationLegend=TRUE) {
  ##first column of mixData should be gene symbol, do not perform log transformation
  ##first column of reference matrix should be gene symbol
  
  mixDataDir <- .testNull(mixDataDir)
  mixData <- .testNull(mixData)
  cellMatrixDir <- .testNull(cellMatrixDir)
  projectName <- .testNull(projectName)
  
  
  if(is.null(projectName)){
			projectName <- gsub("\\.","_",as.character(as.numeric(Sys.time())))
  }
   
  projectDir <- file.path(outputDirectory,paste("Project_",projectName,sep=""))
  dir.create(projectDir)
  
  if(!is.null(mixDataDir)){   ##if the mixData file is not null, use file first
  	mixData <- read.csv(mixDataDir,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
  }else{
  	if(is.null(mixData)){
  		re <- "Error: Please input a file directory of the mix data or upload an R object.\n"
  		cat(re)
  		return(re)
  	}
  }
  
  if(nrow(mixData)<2 || is.null(nrow(mixData))){
  	re <- "Error: There should be at least one sample in the mix data.\n"
  	cat(re)
  	return(re)
  }
  
  if(!is.null(cellMatrixDir)){
  	sigMat <- read.csv(cellMatrixDir,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
  }else{
  	re <- "Error: Please input a file directory of the cell reference matrix.\n"
  	cat(re)
  	return(re)
  }

  colnames(mixData)[1] <- "GeneSymbol"
  colnames(sigMat)[1] <- "GeneSymbol"
  
  ov <- intersect(mixData[,1], sigMat[,1])
  if(length(ov)<10){
  	re <- "Error: The mix data only contain less than 10 genes included in the cell reference matrix. Please check the mix data.\n"
  	cat(re)
  	return(re)
  }

	x <- data.frame(sampleID=rownames(resultMatrix),resultMatrix,stringsAsFactors=F)
  
  write.table(x,file=paste(outputDirectory,"/Project_",projectName,"/CellFraction_",projectName,"_abbCell.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
	
	cat("Finish!\n")
}

.testNull <- function(parameter){
	if(!is.null(parameter)){
		if(length(parameter)==1 && length(which(parameter=="NULL"))==1){
			parameter <- NULL
		}
	}
	return(parameter)
}
