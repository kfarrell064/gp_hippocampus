function pp = fitGPregress_MaxEv(xvals,yvals,fdprs)
% Maximize log evidence for GP regression model with Gaussian noise using
% fourier-domain GP 
%
% pp = fitGPregress_MaxEv(xvals,yvals,fdprs)
%
% INPUT:
% -------
%    xvals [nt x d] - stimulus matrix (each row contains spatial stimulus at single time)
%    yvals [nt x 1] - response vector
%    fdprs (struct) - struct governing Fourier domain representation
%         .minlen [1 x 1] or [m x 1] - minimum length scale for each dimension (can be scalar)
%         .circinterval [2 x d] - circular support in each stimulus dimension
%         .condthresh   [1 x 1] - condition number for thresholding for small eigenvalues 
%
% OUTPUT:
% --------
%   pp - parameter struct with fitted GP model

%% ========= Parse inputs and initialize hyperparams =====================

[nsamps,nd] = size(xvals); % number of samples and input dimensions

% set Fourier-Domain params (fdprs) if necessary
if nargin < 3
    fdprs = initFDprs(xvals);
end

%% ==== Compute sufficient statistics in Fourier domain ==============
[fBasis,wwsq] = mkRBFfourierBasis(xvals,fdprs.minlen,fdprs);
dd.xx = fBasis'*fBasis; % regressors autocovariance
dd.xy = fBasis'*yvals; % cross-covariance of regressors with reponses
dd.yy = yvals'*yvals;
dd.nsamps = nsamps;
wwnrmtot = sum(wwsq,2); % total squared frequencies for each coeff


%% ==== initialize range for hyperparameter grid search =================
% This is just for a  grid search to initialize numerical search. Ranges
% are based on crude heuristics, so definitely room for improvement here.

% lengthscale range
lrange = [min(fdprs.minlen),max(max(fdprs.minlen)*nd,min(range(xvals))/2)];
loglrange = log(lrange);

% Rho range
yvar = dd.yy/nsamps; % variance of output
rhomax = yvar;             % max rho to consider
rhomin = min(1,yvar*.05);  % min rho to consider
rhorange = [rhomin,rhomax];
% Change of variables to tilde rho (which separates rho and length scale)
logtrhorange = log(transformRho(rhorange,lrange,nd,1));

% noise variance sigma_n^2
lognsevarmax = log(yvar); % marginal variance of y
lognsevarmin = log(min(1,yvar*.05));  % minimum to explore
lognsevarrange = [lognsevarmin, lognsevarmax];


%% ========= Grid search for initial hyperparameters =====================

% Make handle for loss function
lfun = @(prs)neglogEv_GPregress(prs,dd,wwnrmtot,fdprs.condthresh);  % loss function

ngrid = 4;  % search a 4 x 4 x 4 grid for initial value of hyperparams
rnges = [loglrange;logtrhorange;lognsevarrange];
[nllvals,gridpts] = grideval(ngrid,rnges,lfun);
[kprs0,~] = argmin(nllvals,gridpts(:,1),gridpts(:,2),gridpts(:,3)); % find minimum


%% ========== Optimize evidence ==========================================

%LB = [log(min(fdp.minlen));-inf;-inf];  % lower bounds
%UB = inf(3,1); % upper bounds
%fminopts = optimset('gradobj','on','Hessian','on','display','iter','algorithm','trust-region-reflective');
% HessCheck(lfun,kprs0);  % check gradient and Hessian numerically, if desired
%kprshat = fmincon(lfun,kprs0,[],[],[],[],LB,UB,[],fminopts); % run optimization

fminopts = optimoptions('fminunc','algorithm', 'trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','display', 'iter');
kprshat = fminunc(lfun,kprs0,fminopts); % run optimization

% Check if length scale is at minimum allowed range
if exp(kprshat(1)) <= min(fdprs.minlen)+.001
    fprintf('\n---------------------------------------------------\n');
    warning(['Solution below minimum length scale; ',...
        'consider reducing ''minlen'' and re-running']);
    fprintf('---------------------------------------------------\n');
end

%% ========= Compute posterior mean and covariance at maximizer of hyperparams ===

% Compute Hessian and Fourier-domain support at posterior mode
[neglogEv,~,H,muFFT,LpostFFT,ii] = lfun(kprshat);

pp.wmean = muFFT;
pp.wcov = LpostFFT;

% Report rank of prior at termination
fprintf('GP evidence optimization terminated with Fourier basis rank = %d\n',sum(ii));

% transform hyperparameters back to original units
[kprshat_untransformed,Hessian] = transformHprs(kprshat,nd,-1,H); % length scale

pp.kprs.len = kprshat_untransformed(1);      % length scale hyperparameter
pp.kprs.rho = kprshat_untransformed(end-1);  % rho hyperparameter
pp.kprs.nsevar = kprshat_untransformed(end); % noise variance
pp.kprspost.Hessian = Hessian; % Hessian of hyperparameters
pp.kprspost.CI = sqrt(diag(inv(Hessian))); % 1SD posterior CI for hyperparams
pp.neglogEv = neglogEv; % negative log-evidence at solution

% Set Fourier-domain params
pp.fdprs = fdprs;

