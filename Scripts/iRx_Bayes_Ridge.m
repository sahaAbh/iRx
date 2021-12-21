function [tau1,lambda1,beta1,tauout, betaout,lambdaout]= iRx_Bayes_Ridge(D,C,P,nrun,burn,thin,seed)
rng(seed,'twister');
%------ Few data characteristics
[p,nc]= size(C);
Tmat = [C, P] ;
M= mean(Tmat,2);VY=var(Tmat');
Tmat = bsxfun(@minus,Tmat,M);                     % center the training data
Tmat = bsxfun(@times,Tmat,1./sqrt(VY'));  
Yc= Tmat(:,1:nc)';  
clear Tmat Train Test C P;


%-----------------------------------
sp =(nrun - burn)/thin;
 betaout = zeros(p,1);
  tauout = 0;
    lambdaout = 0;
tau1 = zeros(sp,1) ;
lambda1= zeros(sp,1);
beta1=   zeros(sp,10);
as=1;bs=0.3;
ad3=1;bd3=0.3;
lambda= gamrnd(ad3,1/bd3);  
 
    %------start gibbs sampling-----%

    for i = 1:nrun
    %-----Update Drug model
    %params----------------------------------------------
     % -- Updatetau -- % tau %--F in our case
        Veta1 = eye(nc) +(lambda*(Yc*Yc'));
        T = cholcov(Veta1); [~,R] = qr(T);
        S = inv(R); Veta = S*S';                   % Veta = inv(Veta1)
        tau = gamrnd(nc/2 +as , 1/((D'*Veta*D)/2+bs));
         %----    Update beta % p x 1 
        Veta1 = lambda*eye(p) +(Yc'*Yc);
        T = cholcov(Veta1); [~,R] = qr(T);
        S = inv(R); Veta = S*S';
        Mbeta = Veta*Yc'*D; 
        beta = Mbeta + (tau^(-0.5)*S)*normrnd(0,1,[p,1]);
         %-----  Update lambda 
        lambda = gamrnd(p/2+ad3, 1/( tau*(beta'*beta)/2+ bd3));
        if mod(i,thin)==0 && i > burn
            betaout = betaout + beta./sp;
            tauout = tauout + tau/sp;
            lambdaout= lambdaout + lambda/sp;
            tau1((i-burn)/thin)= tau;
            lambda1((i-burn)/thin)=lambda;
            beta1((i-burn)/thin,:)= beta([3 12 56 77 123 323 567 611 741 897],1)';
        end

        if mod(i,1000) == 0
            disp(i);
        end
    end

end
