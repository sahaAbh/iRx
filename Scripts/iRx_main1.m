function [iRx,Ni,iRx_std,Ni_std, irx_P,NI_P, AUC_NI, AUC_irx, EP_irx, RP, NRP, idx, L_mh,Lambda_MCMC,tau_MCMC,lambda_MCMC,beta_MCMC, nofout_MCMC,time]= iRx_main(D, C,P,IndResp, nrun,burn,thin,seed)
% -- Abhisek Saha -- %
% This is the same function as that of iRx_main except for the fact that it
% uses iRx_spfact1 which considers the adjustment for badly scaled matrix.
addpath './Scripts';

tic;
[tau_MCMC,lambda_MCMC,beta_MCMC,tauout, betaout,lambdaout]= iRx_Bayes_Ridge(D,C,P,nrun,burn,thin, seed);
[etacout,etapout,Lambdaout,Lambda_MCMC, nof1out,nofout_MCMC,Omegapureout,Omegapout,Omegacout,StdTrain, StdTest, VY]= ...
    iRx_spfact1(C,P,nrun,burn,thin,seed); 
K_Avg = round(mean(nof1out));
clear nof1out;
[p,~]=  size(Lambdaout);
omegap_Avg = reshape(Omegapout,p,p); 
omegac_Avg = reshape(Omegacout,p,p);
AAprimeAvg = Lambdaout(:,1: K_Avg)*Lambdaout(:,1: K_Avg)';
ShrMat =  AAprimeAvg /omegap_Avg; 

predictp_ind = betaout' * StdTest';
predictp_dep = betaout' * ShrMat * StdTest';

Vscur= omegac_Avg -  ShrMat*AAprimeAvg;
sigma_DgivenP_est = sqrt(tauout^(-1)+ betaout' * Vscur * betaout);
clear Omegapout AAprimeAvg omegap_Avg;

IndNotIE= (IndResp <2 ); 
%----------Making Boxplots
%=============== Z score and boxplots and t test
iRx = predictp_dep;
Ni = predictp_ind;
iRx_std= (predictp_dep - mean(predictp_dep))/(std(predictp_dep));
Ni_std =  (predictp_ind -  mean(predictp_ind))/(std(predictp_ind));


boxplot(iRx_std(IndNotIE), IndResp(IndNotIE), 'Labels', {'Non-responder','Responder'})
title('Prediction with iRx model');
saveas(gcf,'iRx_std_boxplot','jpg')
[~,~,~,AUC_irx] = perfcurve(IndResp(IndNotIE),iRx_std(IndNotIE),0);

boxplot(Ni_std(IndNotIE), IndResp(IndNotIE),  'Labels', {'Non-responder','Responder'})
title('Prediction with NI model');
saveas(gcf,'NI_std_boxplot','jpg')
[~,~,~,AUC_NI] = perfcurve(IndResp(IndNotIE),Ni_std(IndNotIE),0);

%================= T-TEST==============================
R= find(IndResp==1);
NR= find(IndResp==0);
irx_R = iRx_std(R);
irx_NR = iRx_std(NR);
NI_R= Ni_std(R);
NI_NR = Ni_std(NR);

[~,irx_P]=ttest2(irx_R,irx_NR,'tail','left','vartype','unequal');
[~,NI_P]=ttest2(NI_R,NI_NR,'tail','left','vartype','unequal');

%===============  Clustering and gene identifications  =======
[idx, L_mh]= iRx_clusters(K_Avg, etacout,etapout, Lambdaout,IndResp,0);
%================  Computing EP scores ==================
N=10000;seed=1000;
[EP_irx, ~, RP, NRP]=  EP(iRx_std,IndResp, 0.5, N ,seed);
time=toc;





