#-------------------------Code for data extraction and processing for iRx model---by  Abhisek Saha---------------------
#=============================================================================================================== 
# Note: part of the code follows paul geeleher's website (http://genemed.uchospitals.edu/$\sim$pgeeleher/natureMed.html), 
# which has been subsequently modified to extract required data and process it to suit our needs.
#========================================================================================================
#----Fix the data folder
#----Fix the data folder to the current one:
#setwd("../Data_Processing/")
library("ridge")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")
library("sva")
library("car")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library("preprocessCore")
library("ROCR")
#BiocManager::install("GEOquery", version = "3.8")
biocLite("GEOquery")
library("GEOquery") 
# Adding the  scripts required .
scriptsDir <- "../scripts"
source(file.path(scriptsDir, "R_scripts.R")) 

#===========================================Data extraction and validation==========================
#Download the expression data and sensitivity of multiple myeloma patients from GEO. The data are already downloaded and strored in "PreprocData/bortezomibData/". To manually download, please uncomment the following code.
load("../PreprocessedData/bortezomibData/bortGeo.RData") 
# bortezomib_mas5 <- getGEO("GSE9782") 
exprDataU133a <- cbind(exprs(bortezomib_mas5[[1]]), exprs(bortezomib_mas5[[2]]))
bortIndex <- c(which(pData(phenoData(bortezomib_mas5[[1]]))[,"characteristics_ch1.1"] == "treatment = PS341"), 255 + which(pData(phenoData(bortezomib_mas5[[2]]))[,"characteristics_ch1.1"] == "treatment = PS341"))
dexIndex <- c(which(pData(phenoData(bortezomib_mas5[[1]]))[,"characteristics_ch1.1"] == "treatment = Dex"), 255 + which(pData(phenoData(bortezomib_mas5[[2]]))[,"characteristics_ch1.1"] == "treatment = Dex"))
studyIndex <- c(as.character(pData(bortezomib_mas5[[1]])[, "characteristics_ch1"]), as.character(pData(bortezomib_mas5[[2]])[, "characteristics_ch1"]))
# The following code should be used instead of the above if getGEO does not produce samples in batches of 255
# exprDataU133a <- exprs(bortezomib_mas5[[1]])
# bortIndex <- which(pData(phenoData(bortezomib_mas5[[1]]))[,"characteristics_ch1.1"] == "treatment = PS341")
# dexIndex <- which(pData(phenoData(bortezomib_mas5[[1]]))[,"characteristics_ch1.1"] == "treatment = Dex")

#U133A array analysis: Map the rownames, which are currently probeset ids
# BiocManager::install("hgu133a.db", version = "2.8")
#source("https://bioconductor.org/biocLite.R")
#biocLite("hgu133a.db")
library("hgu133a.db") # version 2.8.0
x <- hgu133aSYMBOL
mapped_probes <- mappedkeys(x)
names(mapped_probes) <- as.character(x[mapped_probes])
affy2sym <- as.character(x[mapped_probes])
sym2affy <- affy2sym
names(sym2affy) <- names(affy2sym)
rownames(exprDataU133a) <- sym2affy[rownames(exprDataU133a)]

#Load cell lines expression data, which have been already pre-processed using rma function and BrainArray remapped CDF files.
# For detailed codes for pre-processing, please go the folder: "Code_FOR_PreProc_P_Geeleher".
load(file="../PreprocessedData/GdscProcData/gdsc_brainarray_syms.RData")
#Load the Bortezomib drug sensitivity data on cell lines from GDSC website by clicking the ``Sensitivity Data: http://www.cancerrxgene.org/translation/Drug/104
# The fie has been downloaded and stored in "../PreprocData/bortezomibData/" folder.
sensBortezomib <- read.csv("../PreprocessedData/bortezomibData/sensitivity_data_for_drug_104.csv", 
                           as.is=TRUE)
bortic50s <- sensBortezomib$"IC.50"
names(bortic50s) <- sensBortezomib$"Cell.Line.Name"
tissue <- sensBortezomib$"Tissue"
names(tissue) <- sensBortezomib$"Cell.Line.Name"


#Read the GDSC phenotypic data, which maps cell lines to CEL files. 
#Subset to only CEL files which map uniquely to a cell line and order and subset all of the data correctly.
pData <- read.delim("../PreprocessedData/GdscPdata/E-MTAB-783.sdrf.txt", as.is=TRUE)
pDataUnique <- pData[pData$Source.Name %in% names(which(table(pData$Source.Name) ==1)), ]
rownames(pDataUnique) <- pDataUnique$Source.Name
commonCellLines <- rownames(pDataUnique)[rownames(pDataUnique) %in% names(bortic50s)]
pDataUniqueOrd <- pDataUnique[commonCellLines, ]
bortic50sOrd <- bortic50s[commonCellLines]  # has sensitivity info
trainDataOrd <- gdsc_brainarray_syms[, pDataUniqueOrd$"Array.Data.File"] #has expression info

# Reading the respone phenotype
resp133a <- c(as.character(pData(bortezomib_mas5[[1]])[, "characteristics_ch1.8"]),
              as.character(pData(bortezomib_mas5[[2]])[, "characteristics_ch1.8"]))[bortIndex]
IndResp <- rep(1,length(resp133a));
IndResp[resp133a== "PGx_Responder = NR"]=0;
IndResp[resp133a== "PGx_Responder = IE"]=2;

#===================== Data cleaning and processing============================
testExprData= exprDataU133a[, bortIndex];  trainingExprData=trainDataOrd;  trainingPtype=bortic50sOrd;
batchCorrect="eb" # Stands for empirical Bayes option for COMBAT method
powerTransformPhenotype=TRUE # logical parameter indicating power transformation
removeLowVaryingGenes=.2    # proportion of low gene varying genes to be removed
#minNumSamples=10
selection= 1; # summarize duplicate gene ids by their mean  
printOutput=TRUE

# Get the homogenized data
homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
# Do variable selection if specified. By default we remove 20% of least varying genes.
# Otherwise, keep all genes.
if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
{
  keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
}else 
  keepRows <- seq(1:nrow(homData$train))

# PowerTranform phenotype if specified.
offset = 0
if(powerTransformPhenotype)
{
  if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
  {
    offset <- -min(trainingPtype) + 1
    trainingPtype <- trainingPtype + offset
  }
  
  transForm <- powerTransform(trainingPtype)[[6]]
  trainingPtype <- trainingPtype^transForm
}

#===================Saving Data as required========================================================
pval<- QuickPvalue( homData$train[keepRows,], trainingPtype)
Top1k <-  order(pval)[1:1000];
C <- homData$train[keepRows[Top1k],];
P <- homData$test[keepRows[Top1k],];
D <-  trainingPtype;
D_raw <- bortic50sOrd;
Top_P_geneIds <- keepRows[Top1k];
Allnames_CP <- c(names(trainingPtype),colnames(homData$test));
GeneNames<- rownames(homData$train)[Top_P_geneIds] ;
save(C,P, D, D_raw, offset, transForm, IndResp,Top_P_geneIds, keepRows, GeneNames,Allnames_CP, file='DataForFittingMyeloma.Rdata');
library(R.matlab);
writeMat("Input_Myeloma.mat", C=C, 
         P=P, D=D, IndResp=IndResp) 
#===================== End of the subroutine========================         
