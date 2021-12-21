function [etacout,etapout,Lambdaout,Lambda11, nof1out,nofout,Omegapureout,Omegapout,Omegacout,StdTrain, StdTest, VY]= iRx_spfact(C,P,nrun,burn,thin,seed)
% iRx_spfact runs addaptive gibbs sampler and computes the sprase loading
% matrix, and precision matrices. 
rng(seed,'twister');
[p,nc]= size(C);[~,np]= size(P);
Tmat = [C, P] ;
M= mean(Tmat,2);VY=var(Tmat');
Tmat = bsxfun(@minus,Tmat,M);                     % center the training data
Tmat = bsxfun(@times,Tmat,1./sqrt(VY'));  
E1=eig(Tmat*Tmat');
m1=max(E1);m2=min(E1);
Diff= E1(2:end)-E1(1:(end-1));
Jump=(2*(m1-m2))/p;
k1=find(Diff > Jump,1);
Yc= Tmat(:,1:nc)'; Yp=Tmat(:,(nc+1):(nc+np))';
StdTrain = Yc ; StdTest = Yp;
clear Tmat Train Test C P;
% --- define global constants --- %
 
sp =(nrun - burn)/thin; % number of posterior samples

kinit = (p-k1+1)+10;                 % number of factors to start with
b0 = 1; b1 = 0.0005;
epsilon = 1e-3;                                        % threshold limit
prop = 1.00;                                           % proportion of redundant elements within columns


%---- define output files across replicates-----%
  
%    Otc1 = Otc.*(1./sqrt(VY'*VY));                % true dispersion matrix of the transformed data
%    Otp1 = Otp.*(1./sqrt(VY'*VY));
%    Ot_pure = Lamtr*Lamtr';
%    Ot_pure1 = Ot_pure.*(1./sqrt(VY'*VY));                % true dispersion matrix of the transformed data
    
    num = 0;
    k=kinit; 
    % ------end read data--------%

    % --- Define hyperparameter values --- %
    
    as = 1;bs = 0.3;                           % gamma hyperparameters for residual precision
    df = 3;                                    % gamma hyperparameters for t_{ij}
    ad1 = 2.1;bd1 = 1;                         % gamma hyperparameters for delta_1
    ad2 = 3.1;bd2 = 1;                         % gamma hyperparameters delta_h, h >= 2
    adf = 1; bdf = 1;                          % gamma hyperparameters for ad1 and ad2 or df

    % --- Initial values --- %
    ps1= gamrnd(as,1/bs,p,1); Sigmac=diag(1./ps1); 
    ps2= gamrnd(as,1/bs,p,1); Sigmap=diag(1./ps2); 
   
    % Sigma = diagonal residual covariance
    Lambda = zeros(p,k);                                
    psijh = gamrnd(df/2,2/df,[p,k]);                            % local shrinkage coefficients
    delta = ...
        [gamrnd(ad1,bd1);gamrnd(ad2,bd2,[k-1,1])];              % gobal shrinkage coefficients multilpliers
    tauh = cumprod(delta);                                      % global shrinkage coefficients
    Plam = bsxfun(@times,psijh,tauh');                          % precision of loadings rows

    % --- Define output files specific to replicate --- %
    nofout = zeros(nrun+1,1);                  % number of factors across iteartions
    nofout(1) = k;
    nof1out = zeros(sp,1);
     Omegapout = zeros(p^2,1);
  %  Omegap1out = zeros(p^2,1);
    Omegacout = zeros(p^2,1);
  %  Omegac1out = zeros(p^2,1);
     Omegapureout = zeros(p^2,1);
  %    Omegapure1out = zeros(p^2,1);
      etapout = zeros(np, kinit); 
      etacout= zeros(nc, kinit);
      Lambdaout= zeros(p,kinit);
      Lambda11= zeros(p,kinit,sp);   % difference spotted
    %------start gibbs sampling-----%

    for i = 1:nrun
        % -- Update eta -- %  %--F in our case
        Lmsg1 = bsxfun(@times,Lambda,ps1);
        Lmsg2 = bsxfun(@times,Lambda,ps2);
        Veta1 = eye(k) + Lmsg1'*Lambda;
        Veta2 = eye(k) + Lmsg2'*Lambda;
        T = cholcov(Veta1); [Q,R] = qr(T);
        S = inv(R); Veta = S*S';                   % Veta = inv(Veta1)
        Metac = Yc*Lmsg1*Veta;                        % n x k 
        etac = Metac + normrnd(0,1,[nc,k])*S';        % update eta in a block
        T = cholcov(Veta2); [Q,R] = qr(T);
        S = inv(R); Veta = S*S';   
        Metap = Yp*Lmsg2*Veta; 
        etap = Metap + normrnd(0,1,[np,k])*S';
        %eta = [etac ; etap];
        etactemp = etac;etaptemp = etap;
        % -- update Lambda (rue & held) -- %
        eta2c = etac'*etac;
        eta2p = etap'*etap;
        for j = 1:p
            Qlam = diag(Plam(j,:)) + ps1(j)*eta2c+ps2(j)*eta2p;
            blam = ps1(j)*(etac'*Yc(:,j))+ps2(j)*(etap'*Yp(:,j));
            Llam = chol(Qlam,'lower'); zlam = normrnd(0,1,k,1);
         
            if(rcond(Llam) < 1e-12)
               Llam= Llam+eye(k)*0.01;
            end
            vlam = Llam\blam; mlam = Llam'\vlam; ylam = Llam'\zlam;
            Lambda(j,:) = (ylam + mlam)';
        end
        Lambdatemp = Lambda;
        %------Update psi_{jh}'s------%
        psijh = gamrnd(df/2 + 0.5,1./(df/2 + bsxfun(@times,Lambda.^2,tauh')));

        %------Update delta & tauh------%
        mat = bsxfun(@times,psijh,Lambda.^2);
        ad = ad1 + 0.5*p*k; bd = bd1 + 0.5*(1/delta(1))*sum(tauh.*sum(mat)');
        delta(1) = gamrnd(ad,1/bd);
        tauh = cumprod(delta);

        for h = 2:k
            ad = ad2 + 0.5*p*(k-h+1); bd = bd2 + 0.5*(1/delta(h))*sum(tauh(h:end).*sum(mat(:,h:end))');
            delta(h) = gamrnd(ad,1/bd); tauh = cumprod(delta);
        end

        % -- Update Sigma -- %
        Ytilc = Yc - etac*Lambda';
        Ytilp = Yp - etap*Lambda';
        ps1=gamrnd(as + 0.5*nc,1./(bs+0.5*sum(Ytilc.^2)))';
        ps2=gamrnd(as + 0.5*np,1./(bs+0.5*sum(Ytilp.^2)))';
        Sigmac=diag(1./ps1);
        Sigmap=diag(1./ps2);

        %---update precision parameters----%
        Plam = bsxfun(@times,psijh,tauh');

        % ----- make adaptations ----%
        prob = 1/exp(b0 + b1*i);                % probability of adapting
        uu = rand;
        lind = sum(abs(Lambda) < epsilon)/p;    % proportion of elements in each column less than eps in magnitude
        vec = lind >=prop;num = sum(vec);       % number of redundant columns

        if uu < prob
            if  i > 20 && num == 0 && all(lind < 0.995)
                k = k + 1;
                Lambda(:,k) = zeros(p,1);
                etac(:,k) = normrnd(0,1,[nc,1]);
                etap(:,k) = normrnd(0,1,[np,1]);
                psijh(:,k) = gamrnd(df/2,2/df,[p,1]);
                delta(k) = gamrnd(ad2,1/bd2);
                tauh = cumprod(delta);
                Plam = bsxfun(@times,psijh,tauh');
            elseif num > 0
                nonred = setdiff(1:k,find(vec)); % non-redundant loadings columns
                k = max(k - num,1);
                Lambda = Lambda(:,nonred);
                psijh = psijh(:,nonred);
                etac = etac(:,nonred);
                etap = etap(:,nonred);
                delta = delta(nonred);
                tauh = cumprod(delta);
                Plam = bsxfun(@times,psijh,tauh');
            end
        end
        nofout(i+1) = k;

        % -- save sampled values (after thinning) -- %
        if mod(i,thin)==0 && i > burn
            Omegap = Lambda*Lambda' + Sigmap;
        %    Omegap1 = Omegap .* sqrt(VY'*VY);
            Omegac = Lambda*Lambda' + Sigmac;
        %    Omegac1 = Omegac .* sqrt(VY'*VY);
            Omegapure = Lambda*Lambda';
        %    Omegapure1 = Omegapure .* sqrt(VY'*VY);
        %-----------------storing the posterior mean----------    
            Omegapout = Omegapout + Omegap(:)/sp;
        %    Omegap1out = Omegap1out + Omegap1(:)/sp;
            Omegacout = Omegacout + Omegac(:)/sp;
         %   Omegac1out = Omegac1out + Omegac1(:)/sp;
            Omegapureout= Omegapureout + Omegapure(:)/sp;
          %  Omegapure1out= Omegapure1out + Omegapure1(:)/sp; 
            nof1out((i-burn)/thin) = (nofout(i) - num);
        etacout(:,1:nofout(i))= etacout(:,1:nofout(i))+ (etactemp./sp);
        etapout(:,1:nofout(i))= etapout(:,1:nofout(i))+ (etaptemp./sp);
        Lambdaout(:,1:nofout(i))=Lambdaout(:,1:nofout(i))+(Lambdatemp./sp);
        Lambda11(:,1:nofout(i),(i-burn)/thin )= Lambdatemp;

        if mod(i,1000) == 0
            disp(i);
        end
    end

end