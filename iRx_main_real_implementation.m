seed=0;
addpath './MetaData';
load('Input_Myeloma.mat');
nrun=100; burn=50; thin=5;

% Uncomment and run this code if MCMC samples are required. 
%[iRx,Ni,iRx_std,Ni_std, irx_P,NI_P, AUC_NI, AUC_irx, EP_irx, RP, NRP, idx, L_mh,Lambda_MCMC,tau_MCMC,lambda_MCMC,beta_MCMC, nofout_MCMC,time]= iRx_main(D, C,P,IndResp, nrun,burn,thin,seed);
% save('MatOut_Myeloma_full.mat','iRx','Ni','iRx_std','Ni_std','irx_P','NI_P', 'AUC_NI', 'AUC_irx', 'EP_irx','RP', 'NRP','IndResp', 'RP','NRP','idx', 'L_mh','Lambda_MCMC','tau_MCMC','lambda_MCMC','beta_MCMC','nofout_MCMC','time');
[iRx,Ni,iRx_std,Ni_std, irx_P,NI_P, AUC_NI, AUC_irx, EP_irx, RP, NRP, idx, L_mh,~,~,~,~, ~,time]= iRx_main(D, C,P,IndResp, nrun,burn,thin,seed);
save('MatOut_Myeloma.mat','iRx','Ni','iRx_std','Ni_std','irx_P','NI_P', 'AUC_NI', 'AUC_irx', 'EP_irx','RP', 'NRP','IndResp', 'RP','NRP','idx', 'L_mh','time');

load('Input_BreastCancer.mat')
nrun=100; burn=50; thin=5;
warning('off','all');
% Uncomment and run this code if MCMC samples are required. 
%[iRx,Ni,iRx_std,Ni_std, irx_P,NI_P, AUC_NI, AUC_irx, EP_irx, RP, NRP, idx, L_mh,Lambda_MCMC,tau_MCMC,lambda_MCMC,beta_MCMC, nofout_MCMC,time]= iRx_main1(D, C,P,IndResp, nrun,burn,thin,seed);
% save('MatOut_Breastcancer_full.mat','iRx','Ni','iRx_std','Ni_std','irx_P','NI_P', 'AUC_NI', 'AUC_irx', 'EP_irx','RP', 'NRP','IndResp', 'RP','NRP','idx', 'L_mh','Lambda_MCMC','tau_MCMC','lambda_MCMC','beta_MCMC','nofout_MCMC','time');
[iRx,Ni,iRx_std,Ni_std, irx_P,NI_P, AUC_NI, AUC_irx, EP_irx, RP, NRP, idx, L_mh,~,~,~,~, ~,time]= iRx_main1(D, C,P,IndResp, nrun,burn,thin,seed);
save('MatOut_Breastcancer.mat','iRx','Ni','iRx_std','Ni_std','irx_P','NI_P', 'AUC_NI', 'AUC_irx', 'EP_irx','RP', 'NRP','IndResp', 'RP','NRP','idx', 'L_mh','time');
