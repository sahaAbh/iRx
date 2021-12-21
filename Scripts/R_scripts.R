library("sva")
library("preprocessCore")

## Copyright: Paul Geeleher
## 

## a function to do variable selection on expression matrices.
## I.e. remove genes with low variation
## It returns a vector of row ids to keep. Note, rownames() must be specified.
doVariableSelection <- function(exprMat, removeLowVaryingGenes)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}


## Copyright: Paul Geeleher
## All rights Reserved
homogenizeData <- function(testExprMat, trainExprMat, batchCorrect="eb", selection=-1, printOutput=TRUE)
{
  # Take two expression matrices from separate batches, or microarray platforms, and given that the have the same gene ids, return homogenized versions of the matrices.
  # Args:
  #    testExprData: Gene expression matrix for samples on which we wish to predict a phenotype. Gene names as rows, samples names as columns.
  #    trainingExprData: Gene expression matrix for samples for which we the phenotype is already known.
  #    batchCorrect" The type of batch correction to be used. Options are "eb", "none", .....
  #    selection: Parameter can be used to specify how duplicates are handled, by default value NULL means ask the user.
  # Returns:
  #    A list containing two entries $train and $test, which are the homogenized input matrices.
  
  
  # Check batchCorrect paramter
  if(!(batchCorrect %in% c("eb", "qn", "none")))
    stop("\"batchCorrect\" must be one of \"eb\", \"qn\" or \"none\"")
  
  # check if both row and column names have been specified
  if(is.null(rownames(trainExprMat)) || is.null(rownames(testExprMat)))
  {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both training and test expression matrices. Both matices must have the same type of gene identifiers.")
  }
  
  # check that some of the row names overlap between both datasets (print an error if none overlap.
  if(sum(rownames(trainExprMat) %in% rownames(testExprMat)) == 0)
  {
    stop("ERROR: The rownames() of the supplied expression matrices do not match. Note that these are case-sensitive.")
  }
  else
  {
    if(printOutput) cat(paste("\n", sum(rownames(trainExprMat) %in% rownames(testExprMat)), " gene identifiers overlap between the supplied expression matrices... \n", paste=""));
  }
  
  # if there are duplicate gene names, give the option of removing them or summarizing them by their mean.
  if((sum(duplicated(rownames(trainExprMat))) > 0) || sum(sum(duplicated(rownames(testExprMat))) > 0))
  {
    if(selection == -1) #print the following if we're asking the user how to proceed.
    {  
      cat("\nExpression matrix contain duplicated gene identifiers (i.e. duplicate rownames()), how would you like to proceed:")
      cat("\n1. Summarize duplicated gene ids by their mean value (acceptable in most cases)")
      cat("\n2. Disguard all duplicated genes (recommended if unsure)")
      cat("\n3. Abort (if you want to deal with duplicate genes ids manually)\n")
    }
    
    while(is.na(selection) | selection <= 0 | selection > 3 )
    {
      selection <- readline("Selection: ")
      selection <- ifelse(grepl("[^1-3.]", selection), -1 , as.numeric(selection))
    }
    
    cat("\n")
    
    if(selection == 1) # summarize duplicates by their mean
    {
      if((sum(duplicated(rownames(trainExprMat))) > 0))
      {
        trainExprMat <- summarizeGenesByMean(trainExprMat)
      }
      if((sum(duplicated(rownames(testExprMat))) > 0))
      {
        testExprMat <- summarizeGenesByMean(testExprMat)
      }
    }
    else if(selection == 2) # disguard all duplicated genes
    {
      if((sum(duplicated(rownames(trainExprMat))) > 0))
      {
        keepGenes <- names(which(table(rownames(trainExprMat)) == 1))
        trainExprMat <- trainExprMat[keepGenes, ]
      }
      
      if((sum(duplicated(rownames(testExprMat))) > 0))
      {
        keepGenes <- names(which(table(rownames(testExprMat)) == 1))
        testExprMat <- testExprMat[keepGenes, ]
      }      
    }
    else
    {
      stop("Exectution Aborted!")
    }
    
  }
  
  # subset and order gene ids on the expression matrices
  commonGenesIds <- rownames(trainExprMat)[rownames(trainExprMat) %in% rownames(testExprMat)]
  trainExprMat <- trainExprMat[commonGenesIds, ]
  testExprMat <- testExprMat[commonGenesIds, ]
  
  # subset and order the two expresison matrices
  if(batchCorrect == "eb")
  {
    # library("sva")
    
    # subset to common genes andbatch correct using ComBat
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame("(Intercept)"=rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    combatout <- ComBat(dataMat, whichbatch, mod=mod)
    return(list(train=combatout[, whichbatch=="train"], test=combatout[, whichbatch=="test"], selection=selection))
  }
  else if(batchCorrect == "qn")
  {
    # library("preprocessCore")
    dataMat <- cbind(trainExprMat, testExprMat)
    dataMatNorm <- normalize.quantiles(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    return(list(train=dataMatNorm[, whichbatch=="train"], test=dataMatNorm[, whichbatch=="test"], selection=selection))
  } 
  else
    return(list(train=trainExprMat, test=testExprMat, selection=selection))
}


## Paul Geeleher
## All rights Reserved
## February 6, 2014

## This function will take a matrix with duplicate rownames()
## and summarize those rownames() by their mean.
summarizeGenesByMean <- function(exprMat)
{
  geneIds <- rownames(exprMat)
  t <- table(geneIds) # how many times is each gene name duplicated
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]
  
  # create a *new* gene expression matrix with everything in the correct order....
  # start by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), ]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]
  
  # add all the duplicated genes to the bottom of "exprMatUniqueHuman", summarizing as you go
  for(numDups in allNumDups) 
  {
    geneList <- names(which(t == numDups))
    
    for(i in 1:length(geneList))
    {
      exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == geneList[i]), ]))
      gnamesUnique <- c(gnamesUnique, geneList[i])
      #print(i)
    }
  }
  
  rownames(exprMatUnique) <- gnamesUnique
  return(exprMatUnique)
}


## Copyright: Abhisek saha
## Computes pvalues of a gene with output variable
QuickPvalue<-function(trainingExprData, trainingPtype, printOutput=F)
{ 
  pval=rep(0,nrow(trainingExprData))
  for(ii in 1:nrow(trainingExprData)){
  trainFrame <- data.frame(Resp=trainingPtype, trainingExprData[ii,])
  lmModel <- lm(Resp ~ ., data=trainFrame)
  pval[ii]<-summary(lmModel)$coefficients[2,4]}
  if(printOutput) cat("Done\n\nCalculating predicted phenotype...");
  return(pval)
}


## Copyright: Abhisek saha

ClusterDist<-function(ClusterData, Col1=1, Col2=5, ClusterID= 11 ,Tag)
{
ClusterTemp=ClusterData[ClusterData[,Col1]==ClusterID,];
xmin= min(ClusterTemp[,5])-1;xmax= max(ClusterTemp[,5])+1;
hist(ClusterTemp[ClusterTemp$V3==3,Col2], col=rgb(0,0,1,0.5),xlim=c(xmin,xmax), main=Tag ,xlab='Sensitivity')
hist(ClusterTemp[ClusterTemp$V3!=3,Col2], col=rgb(0,1,0,0.5), add=T)
legend('topright',legend=c('Cell lines','Patients'),fill=c('blue','green'))
}

## Copyright: Abhisek saha
CircClust<- function(Data, Cutoff= 45, sense=5)
{
S1= 1*(Data[,sense] > Cutoff); 
S1[Data[,sense] < -1*Cutoff ]=-1;
Data1=cbind(Data, Sen= S1);
TotCnt<-aggregate(V2~Sen+V1,data=Data1,"length") # Gives total size
Data2<-reshape(as.data.frame(TotCnt), direction="wide", v.names = 'V2', timevar = "Sen",
        idvar = "V1")
rownames(Data2)= paste0('G',Data2[,1]);
Data3= Data2[Mixed1,];
print(Data3)
Data3[is.na(Data3)]=0;
names.arg= paste0('C',1:11);
#names.arg[Cl1]=Mixed2;
titl1 = paste0("Split by Cluster, cutoff=",Cutoff);
titl2 = "Split by Cluster";
barplot2(t(as.matrix(Data3[,-1])), beside = F ,names.arg=names.arg,xlab="Cluster Index",main= "Split by Cluster", ylab="Frequency",col=c("blue","grey","red","beige"))
legend("topright",legend=c("Responder","Neutral","Non-responder","Cell lines"),fill=c("blue","grey","red","beige"))
}

## Copyright: Abhisek saha
ClstCirc<- function(idx5Sum.1, idx4.1, RespTag=c("Responder", "Neutral", "Non-responder"), ALLNames=ALLNames)
{
# requires circlize library:  
#----idx4.1 MUST have V1 (CLUSTER INDEX), V2(PATIENT INDEX), V3(RESPONDER INDEX), V4(SAMPLE INDEX)
#----idx5Sum.1  MUST have V1 (CLUSTER INDEX), V2(TOTAL SAMPLE SIZE), V3(CELL LINE SIZE), V4(PATIENT RESP SIZE), V4(PATIENT NON-RESP SIZE) 
PCONLY = idx5Sum.1[ as.logical((!((idx5Sum.1$V4==0) * (idx5Sum.1$V5==0))) * (idx5Sum.1$V3 >0)), 1 ]  
idx4.PC = idx4.1[idx4.1$V1 %in% PCONLY ,];
idx4.C = idx4.PC[idx4.PC$V2==1,c("V1","V4")]; idx4.C<-idx4.C[order(idx4.C$V1,idx4.C$V4),]
idx4.P = idx4.PC[idx4.PC$V2==0,c("V1","V3","V4")]; idx4.P<-idx4.P[order(idx4.P$V1,idx4.P$V4),]
Adj<-merge(idx4.C,idx4.P,by="V1",suffixes=c(".C",".P"))
Rnk<- cbind(Rank=1:length(PCONLY),V1=PCONLY)
Adj1<-merge(Rnk,Adj,by="V1");Adj1<-Adj1[order(Adj1$Rank,Adj1$V4.C,Adj1$V4.P),]
#RespTag<-c("Non-responder","Responder")
#df = data.frame(from = AllNames[Adj1$V4.C],
#                to = RespTag[Adj1$V3+1],
#                value = 1,
#                stringsAsFactors = FALSE)
Adj2<-aggregate(V4.P~V1+Rank+V4.C+V3,data=Adj1,length);Adj2<-Adj2[order(Adj2$Rank,Adj2$V4.C,Adj2$V3),]
df1 = data.frame(from = AllNames[Adj2$V4.C],
                 to = RespTag[Adj2$V3+2],
                 value = Adj2$V4.P,
                 stringsAsFactors = FALSE)
ClSize<-idx5Sum.1[as.character(PCONLY),"V3"]
circos.par(start.degree = 180,gap.after = c(rep(0,ClSize[1]-1),5,rep(0,ClSize[2]-1),5,rep(0,ClSize[3]-1),5,rep(0,ClSize[4]-1),5,
                                            rep(0,ClSize[5]-1),5,rep(0,ClSize[6]-1),5,rep(0,ClSize[7]-1),5,5,5,5),clock.wise = FALSE)
chordDiagram(df1, grid.col=c(rep("grey",194),rep("purple",3)),col = Adj2$Rank,link.border= Adj2$Rank, annotationTrack = c("grid",axis),transparency=0.2, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df1))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if((xplot[2] - xplot[1]) < -20) {
    circos.text(mean(xlim), mean(ylim), sector.name, facing = "outside",
                niceFacing = TRUE,  col = "black")
    print((xplot[2] - xplot[1]))
  }else{
    if((xplot[2] - xplot[1]) < -2) {
      circos.text(mean(xlim), mean(ylim), sector.name, facing = "clockwise",
                  niceFacing = TRUE,  col = "black")
      print((xplot[2] - xplot[1]))
    }
  } 
}, bg.border = NA)
circos.clear()
}



