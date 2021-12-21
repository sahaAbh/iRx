
function  Lambdat = genLambda(p, k, sparsity_prob, lambda_true_sigma)
% DOCUMENTATION:  REQUIRED INPUTS
%   genLambda(p,k,sparsity_prob, lambda_true_sigma) computes the Lambda matrix (also called the loading
%   matrix) with p many genes and k many latent variables. The
%   sparisity_prob vector specified the the vector of proportions of non-zero's
%   starting from the leftmost columns
%   vector 
% OPTIONAL INPUT:
%   lambda_true_sigma = sigma parameter for the normal (default = 1)
% Written by Abhisek Saha (abhisek.saha@uconn.edu)    
%--------------------------------------------------------------------------
if nargin < 2, k  = 10; end
if nargin < 3,  sparsity_prob= 1:(-0.05):(1-(k-1)*0.05);  end
if nargin < 4, lambda_true_sigma =1; end
%----Some helps----------------
% make sure to include the following for compiling
%  example--------------------------------------------------------------------------
%if ischar(reps),          eval(sprintf('reps = %s;',reps)),           end
%   
% see also mcc2, generate_swarm_kmeans_image, par_kmeans_imag

if length(sparsity_prob) ~= k, error("ERROR: \n sparsity_prob must of length %d ",k); end   

% numeff = # non-zero entries in each column of Lambda
% numeff varies 

Lambdat = zeros(p,k);   % A in my case, 
randperm(p);

%sparsity_prob= 1:(-0.05):.30;
numeff= floor(sparsity_prob*p);

% generate loadings
for h = 1:k
    temp = randsample(p,numeff(h));
    Lambdat(temp,h) = normrnd(0,lambda_true_sigma,[numeff(h),1]);
end
end

 