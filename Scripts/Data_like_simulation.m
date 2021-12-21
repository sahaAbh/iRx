function[dat_c,dat_p,D,D2,seed,trc]=Data_like_simulation(seed,psi_1t,psi_2t,mp,p,k,nc,np)
% Data_like_simulation(seed,psi_1t,psi_2t,mp,p,k,nc,np) generates data for
% fixed seed with source specific precisions, psi_1 and psi_2. mp controls
% the percentage of total variance explained by common source. 
%--OTHER (OPTIONAL) INPUTS 
% p = number of genes ( defaults to 1000)
% k = number of factors (defaults to 15)
% nc = number of celllines (defaults to 300)
% np = number of patient-samples (defaults to 200)
%--written by Abhisek Saha (abhisek.saha@uconn.edu)----

%save('data_like_simulation_2_1.mat')
%sim_data = 'data_like_simulation.mat';
% data like simuation
%load('Mat3_25k_Combined.mat','Lambdaout','nof1out');
%K_Avg = round(mean(nof1out));
% Lambdat= Lambdaout(:,1: K_Avg)
%[p,~]=  size(Lambdaout);
% t1 = clustergram(AAprimeAvg); takes some time to plot
% Simulation from known Lambda

%---------------------------Setting Defaults INPUTS------------------------
if nargin < 7,  nc = 300; np= 200; end
if nargin < 6,  k = 15; end
if nargin < 5,  p = 1000; end
if nargin < 4,  mp = 1; end
if nargin < 2,  psi_1t = 1.5;psi_2t = 1; end
%--------------------------------------------------------------------------
sig_mean = sqrt((psi_1t^(-1) + psi_2t^(-1))/2);
sparsity_prob= 1:(-0.05): (1-0.05*(k-1));
rng(seed,'twister');
Lambdat = sqrt(mp*sig_mean) * genLambda(p,k,sparsity_prob);

AAprimeAvg = Lambdat*Lambdat';
%AAprimeAvg(1:3,1:3   %print check
rep = 1;  Nc = nc*rep; Np =np*rep;

Otc = AAprimeAvg + (psi_1t^(-1))* eye(p);
Otp = AAprimeAvg + (psi_2t^(-1))* eye(p);
%Otc(1:3,1:3) %print check
%Otp(1:3,1:3)  %print check

%---generate P,C------------
mu = zeros(1,p);
dat_c =  mvnrnd(mu,Otc,Nc); 
dat_p =  mvnrnd(mu,Otp,Np);

%-------generate beta from old analysis----------
%load('posterior_drug.mat', 'betaout');
%beta2= betaout;
%  hist(betaout)

%-------generate beta from normal----------
beta1= normrnd(0,1,p,1);
P= dat_p';C= dat_c';
Tmat = [C, P] ;
M= mean(Tmat,2);VY=var(Tmat');
Tmat = bsxfun(@minus,Tmat,M);                     % center the training data
Tmat = bsxfun(@times,Tmat,1./sqrt(VY'));  
Yc= Tmat(:,1:nc)'; Yp=Tmat(:,(nc+1):(nc+np))';
%StdTrain = Yc ; StdTest = Yp;
clear Tmat  C P;
 Otc1 = Otc.*(1./sqrt(VY'*VY));                % true dispersion matrix of the transformed data
 Otp1 = Otp.*(1./sqrt(VY'*VY));
AAprimeAvg1 = AAprimeAvg.*(1./sqrt(VY'*VY));
%AAprimeAvg1(1:3,1:3)
% Otc1(1:3,1:3)
mu_c= Yc * beta1;
%mu_c= Yc * beta2;
% real sigma 1.8
D= mu_c+ normrnd(0,2,nc,1);
%---------true sensitivity data for patients------
%----Defining V: variation due to cell line model
Vc = Otc1-AAprimeAvg1*(Otp1\AAprimeAvg1);  
sigma_p = sqrt(4+ beta1'*Vc*beta1); % approximately 32.82
mu_p= Yp * (Otp1\AAprimeAvg1) * beta1;
D2= mu_p+ normrnd(0,sigma_p,np,1);
a1 =trace(AAprimeAvg1)/trace(Otp1);
a2= trace(AAprimeAvg1)/trace(Otc1);
a3= trace(AAprimeAvg1);
a4= trace(Otp1);
a5=trace(Otc1);
trc = [a1,a2,a3,a4,a5];



%save('data_like_simulation_2_1_temp.mat')
%==========Myeloma specific simulation==================
%========== Sparsity of Lambda================
%sum(Lambdaout(:,1:11)>0.001)/10

