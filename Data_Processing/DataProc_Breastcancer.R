#-------------------------Code for data extraction and processing for iRx model---by  Abhisek Saha---------------------
#=============================================================================================================== 
# Note: part of the code is  taken from paul geeleher's website (http://genemed.uchospitals.edu/$\sim$pgeeleher/natureMed.html), 
# which has been subsequently modified to extract required data and process it to suit our needs.
#========================================================================================================
# Please install all the libraries before they are executed.
#----Fix the data folder to the current one:
#setwd("../Data_Processing/")
library("ridge")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library("sva")
library("car")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
library("preprocessCore")
library("ROCR")

# Adding the  scripts required .
scriptsDir <- "../Scripts"
#source(file.path(scriptsDir, "compute_phenotype_function.R"))
#source(file.path(scriptsDir, "summarizeGenesByMean.R"))
#source(file.path(scriptsDir, "homogenize_data.R"))
#source(file.path(scriptsDir, "do_variable_selection.R")) 
source(file.path(scriptsDir, "R_scripts.R")) 

#===========================================Data extraction and validation==========================
#Extract the expression data for breast cancer patients from GEO. The data are already pre-processed after being downloaded from GEO website, using rma function and  and BrainArray remapped CDF files. 
#The data is then strored in "PreprocessedData/docetaxelData/". One can manually download and pre-process the data following codes as given in "/Codes_FOR_PreProc_P_Geeleher" due to Paul Geeleher.
load("../PreprocessedData/docetaxelData/doce_rma_syms_brainArray.RData") 

#Load cell lines expression data, which have been already pre-processed using rma function and BrainArray remapped CDF files.
# For detailed codes for pre-processing, please go the folder: ""../PreprocessedData/Code_FOR_PreProc_P_Geeleher/".
load(file="../PreprocessedData/GdscProcData/gdsc_brainarray_syms.RData")

#Load the Docetaxel drug sensitivity data on cell lines from GDSC website by clicking the ``Sensitivity Data: http://www.cancerrxgene.org/translation/Drug/1007
# The fie has been downloaded and stored in "../PreprocData/docetaxelData/" folder.
sensDoce <- read.csv("../PreprocessedData/docetaxelData/sensitivity_data_for_drug_1007.csv", 
                           as.is=TRUE)
doceic50s <- sensDoce$"IC.50"
names(doceic50s) <- sensDoce$"Cell.Line.Name"


#Read the GDSC phenotypic data, which maps cell lines to CEL files. 
#Subset to only CEL files which map uniquely to a cell line and order and subset all of the data correctly.
pData <- read.delim("../PreprocessedData/GdscPdata/E-MTAB-783.sdrf.txt", as.is=TRUE)
pDataUnique <- pData[pData$Source.Name %in% names(which(table(pData$Source.Name)  == 1)), ]
rownames(pDataUnique) <- pDataUnique$Source.Name
commonCellLines <- rownames(pDataUnique)[rownames(pDataUnique) %in% names(doceic50s)]
pDataUniqueOrd <- pDataUnique[commonCellLines, ]
doceic50sOrd <- doceic50s[commonCellLines]
trainDataOrd <- gdsc_brainarray_syms[, pDataUniqueOrd$"Array.Data.File"]

# Reading the respone phenotype
# Among the 24  patient samples collected, first 10 are responder (sensitive) and rest 14 are non-responder (resistant) 


#===================== Data cleaning and processing============================
testExprData= doceVivoNorm_syms;  trainingExprData=trainDataOrd;  trainingPtype=doceic50sOrd;
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

# In this study, first 10 people are responders:
IndResp <- rep(1,24);
IndResp[11:24]=0;
#===================Saving Data as required========================================================
pval<- QuickPvalue( homData$train[keepRows,], trainingPtype)
Top1k <-  order(pval)[1:1000];
C <- homData$train[keepRows[Top1k],];
P <- homData$test[keepRows[Top1k],];
D <-  trainingPtype;
D_raw <- doceic50sOrd;
Top_P_geneIds <- keepRows[Top1k];
Allnames_CP <- c(names(trainingPtype),colnames(homData$test));
GeneNames<- rownames(homData$train)[Top_P_geneIds] ;
save(C,P, D, D_raw, offset, transForm, IndResp,Top_P_geneIds, keepRows, GeneNames,Allnames_CP, file='DataForFittingBreastcancer.Rdata');
library(R.matlab);
writeMat("Input_Breastcancer.mat", C=C, 
         P=P, D=D, IndResp=IndResp) 
#===================== End of the subroutine========================         
