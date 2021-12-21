function [idx4, L_mh]= iRx_clusters(K_Avg, etacout,etapout, Lambdaout,IndResp,seed)
rng(seed, 'twister');
%load('current_all','Lambdaout','etacout','etapout');
[n1,~]= size(etacout); [n2,~]=size(etapout);n=n1+n2;
X= [etacout(:,1:K_Avg) ; etapout(:,1:K_Avg)];
[idx,~] = kmeans(X,K_Avg);
%----Cluster wise distribution
idx2= [idx,(1:n)'];
Ind= [ones(n1,1);zeros(n2,1)];
% 1 responder 0 non responder 2 others
IndResp2= 3*Ind;
IndResp2(IndResp2==0)=IndResp;
idx3= [idx,Ind,IndResp2, (1:n)'];
idx4=sortrows(idx3,1);

clear nof1out;
%[p,~]=  size(Lambdaout);
Loadings = Lambdaout(:,1: K_Avg);
Load2= Loadings.^2;
Scale= sum(Load2,1);
L_mh = bsxfun(@times,Load2,1./Scale);
