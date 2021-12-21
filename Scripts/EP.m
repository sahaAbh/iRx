function [Prob1, Prob2, RP, NRP]=  EP(centered_score,TruClu, sig, N ,seed)
% Helps compute enrichment probability and other percentages 
rng(seed);
IndNotIE2= (TruClu <2);
TruClu= TruClu(IndNotIE2);
centered_score=centered_score(IndNotIE2);
centered_score= centered_score- mean(centered_score);

theta0_param  = [sum(centered_score(logical(TruClu)) > sig)+1, sum(centered_score(logical(~TruClu)) > sig)+1];%prediced nonresponder
theta1_param  = [sum(centered_score(logical(TruClu)) < -1*sig)+1, sum(centered_score(logical(~TruClu)) < -1*sig)+1]; %predicted responder
theta0 = drchrnd( theta0_param , N);
theta1 = drchrnd( theta1_param , N);

Predicted=2*ones(numel(centered_score),1);
Predicted(centered_score <  -1*sig)=1;
Predicted(centered_score >  sig)=0;
RP = mean(TruClu(Predicted==1));
NRP = 1- mean(TruClu(Predicted==0));

Prob1=  mean(theta1(:,1) > theta0(:,1));
Prob2=  mean(theta0(:,2) > theta1(:,2));
