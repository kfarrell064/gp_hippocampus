function check_LogEvidenceHessianAtMode
% Check calculation of log-evidence and Hessian of negative log-evidence w.r.t. the hyperparameters

if contains(pwd,'/unittests')
    cd ..; setpaths; % move up a directory and call setpaths
    addpath unittests;
end
clf;

%% 1. Sample data from a multi-D GP
% =================================

nd = 1; % dimensionality of inputs

% Set hyperparameters for covariance function 
rhox = 11.2; % marginal variance
lenx = 2.2; % length scale
signse = 1.5; % stdev of observation noise
nsamp = 1e3; % number of observations

% Sample points for collecting data
xrnge = [0,20]; % support
xsamp = rand(nsamp,nd)*diff(xrnge)+xrnge(1); % x locations

% Set parameters governing Fourier-domain representation
nSTD = 5; % number of length scales to extend circular interval
fdprs.minlen = 1.5;  % minimum length scale to consider
fdprs.circinterval = [xrnge(1);xrnge(2)+nSTD*lenx]*ones(1,nd); % circular interval for function
fdprs.condthresh = 1e8; 

% Sample function using Fourier representation
[cdiag,Bsamp] = Krbf_fourier(xsamp,lenx,rhox,fdprs);
nfreq = length(cdiag);
fsamp = Bsamp*(sqrt(cdiag).*randn(nfreq,1));
ysamp = fsamp+signse*randn(nsamp,1);

% Plot samples for 1st 2 dimensions
plot(xsamp,fsamp,'o',xsamp, ysamp, '.'); 
title('samples from GP and noisy observations');
xlabel('x'); ylabel('y');


%% 2. Optimize log-evidence to infer hyperparameters
% =================================

pp = fitGPregress_MaxEv(xsamp,ysamp,fdprs);

% fprintf('\nLearned hyperparameters (+/- 2SD):\n----------------------------------\n');
% fprintf('len:    %.2f  (+/- %.2f) [true=%2.2f]\n', pp.kprs.len, 2*pp.kprspost.CI(1),len);
% fprintf('rho:    %.2f  (+/- %.2f) [true=%2.2f]\n', pp.kprs.rho, 2*pp.kprspost.CI(2),rho);
% fprintf('nsevar: %.3f (+/- %.2f) [true=%2.3f]\n', pp.kprs.nsevar, 2*pp.kprspost.CI(3),signse.^2);


%% 3. Check evaluation of log-evidence at mode
% ============================================

% This is only checking that the function compLogEv_GP computes evidence in
% the same way as the function neglogEv_GPregress.m, which is used during
% optimization.  Both use Fourier domain formula

logev1 = -pp.neglogEv; 
logev2 = compLogEv_GP(xsamp,ysamp,pp.kprs.len,pp.kprs.rho,pp.kprs.nsevar,pp.fdprs);

fprintf('\nValidating computation of log evidence:\n=======================================\n');
fprintf('optimization code: %.4f\n',logev1);
fprintf('computed directly: %.4f\n',logev2);
fprintf('             Diff: %.6f\n\n', logev1-logev2);

if  abs(logev1-logev2)<1e-6
    fprintf('unit test 1: PASSED\n');
else
    warning('unit test 1: FAILED');
end


%% Check Hessian from optimization code with Hessian computed numerically 
% =======================================================================
%
% NOTE: relies on function 'compHess.m' from ncclabcode

fprintf('\nValidating Hessian of neg-log-evidence w.r.t. hyperparams (len, rho, nsevar)\n');
fprintf('==============================================================================\n');

% Extract from code
Hessian1 = pp.kprspost.Hessian

% Compute numerically
lfun = @(x)(-compLogEv_GP(xsamp,ysamp,x(1),x(2),x(3),pp.fdprs));
kprsMAP = [pp.kprs.len;pp.kprs.rho;pp.kprs.nsevar];
Hessian2 = compHess(lfun,kprsMAP,1e-4)

Diff = Hessian1-Hessian2

fprintf('max abs diff = %.4f\n\n', max(abs(Hessian1(:)-Hessian2(:))));

if  max(abs(Hessian1(:)-Hessian2(:)))< 0.25
    fprintf('unit test 2: PASSED\n');
else
    warning('unit test 2: FAILED');
end
