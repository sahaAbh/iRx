# iRx: A Bayesian Precision Medicine Framework for Calibrating Individualized Therapeutic Indices in Cancer

Please visit [BioRxiv link](https://www.biorxiv.org/content/10.1101/2021.08.09.455722v1) to have full access to the manuscript, and to know more details on the notations, and contexts the attached code refers to. We have demontrataed below how to reprodue the key figures. A slightly more detailed version can be found at [Dr. Veera's lab](https://github.com/bayesrx/iRx).

# REPRODUCTION OF THE ANALYSIS

## STEP 0: SOFTWARE ENVIRONMENT REQUIREMENT:
Both R and Matlab softwares are utilized for programming depending on the type of functionality required. R is used for data extraction and processing.  MATLAB is then used to implement iRx model, i.e., to generate parameter estimates, calibrated scores;  to perform clustering  on factors, and  to compute performance measures, such as AUC, EP score for iRx and NI methods.   R is again used for visualizing results from Matlab output as well as  for  implementation of PC_ind, PC_joint.  Thus, we break the process of implementing iRx  miodel into 3 steps .  1  Data Extraction And Processing    2. iRx Model Implementation   3.  Implementation of PC_ind and PC_joint + Visualizing  Outputs. Details of R and Matlab codes are listed in the table attached (Table_of_functions.docx).   

## STEP 1: DATA EXTRACTION AND PROCESSING: 
FOR MULTIPLE MYELOMA
1.  Run  “./Data_And_Codes/ Data_Processing / DataProc_Myeloma.R”   for Multiple Myeloma (MM) data. 
2. Obtain the processed data in the same folder:  (a) Input_Myeloma.mat, (b)  DataForFittingMyeloma.Rdata. A copy of them are provided in “./Data_And_Codes /MetaData/”

FOR BREAST CANCER
1.  Run  “./Data_And_Codes/ Data_Processing / DataProc_Breastcancer.R” for Breast cancer(BC) data.
2. Obtain the processed data in the same folder:  (a) Input_Breastcancer.mat, (b)  DataForFittingBreastcancer.Rdata. A copy of them are provided in “./Data_And_Codes /MetaData/”

## STEP 2: iRx IMPLEMENTATION:  
1. Run  “./Data_And_Codes/iRx_main_real_implementation.m”   
2. Obtain two “.mat” files. (a) MatOut_Myeloma  (b) MatOut_Breastcancer.  A copy of them are provided in “./Data_And_Codes /MetaData/”.
 
## STEP 3:   Implementation of PC_ind, PC_joint + VISUALIZING OUTPUTS 
1. Run  “./Data_And_Codes/Output_Visualization/Myeloma.R”   for Myeloma study.
2. Run  “./Data_And_Codes/Output_Visualization/Breastcancer.R” for Breast cancer study

The aforementioned programs in this folder may use auxiliary subroutines and pre-downloaded data-sets from other folders, as specified in the “Requirements and Outputs” column in the table below. Each program is heavily commented for being self-explanatory. To execute them smoothly, please load whole folder of Data_And_Codes. 
    

## SIMULATION RESULTS AND PLOTS
We briefly give an outline of the process to obtain simulation results. The data generation process is straightforward (See section 4 of the main manuscript) and any statistical software can be used. We have used MATLAB for this purpose and have supplied a function, “./Scripts/Data_like_simulation.m”  which can be used to generate data replicates for various settings. These function are heavily commented and inputs and outputs are self-explanatory from the context.  Then, one can fit iRx and NI using “iRx_main.m” as used in real data settings earlier.  To obtain predicted values of PC_ind and PC_joint one can make use two subroutines, “pcr_pred” (for ind) and “pcr_pred2” from  “./Scripts/PCR_ind_joint.R” , After that computing the performance measure is straight-forward. We have also supplied “./Scripts/rand_index.m” to compute RI score.  
