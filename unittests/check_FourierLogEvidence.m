function check_FourierLogEvidence
% Check calculation of log-evidence using Fourier domain compared to classic formula

if contains(pwd,'/unittests')
    cd ..; setpaths; % move up a directory and call setpaths
    addpath unittests;
end
clf;

% 1. Sample data from a multi-D GP
% =================================

nd = 1; % dimensionality of inputs

% Set hyperparameters for covariance function 
rhox = 11.2; % marginal variance
lenx = 1.2; % length scale
signse = 1.5; % stdev of observation noise
nsamp = 100; % number of observations

% Sample points for collecting data
xrnge = [0,20]; % support
xsamp = rand(nsamp,nd)*diff(xrnge)+xrnge(1); % x locations

% Set parameters governing Fourier-domain representation
nSTD = 10; % number of length scales to extend circular interval
fdprs.minlen = lenx;  % minimum length scale to consider
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


% 2. Compute log-evidence
% =======================

rhotest = rand*rhox*2+0.01;
lentest = rand*lenx*2+0.5;
nsevartest = rand*signse.^2*2+.01;

% compute using Fourier domain representation
logev_fd = compLogEv_GP(xsamp,ysamp,lentest,rhotest,nsevartest,fdprs)

% compute using classic formula

Lcov = Krbf(xsamp,lentest,rhotest)+nsevartest*eye(nsamp); % covariance
logev_classic = logmvnpdf(ysamp',zeros(1,nsamp),Lcov)

fprintf('\nValidating computation of log evidence:\n=======================================\n');
fprintf('Fourier:  %.4f\n',logev_fd);
fprintf('classic:  %.4f\n',logev_classic);
fprintf('   Diff: %.6f\n\n', logev_fd-logev_classic);

if  abs(logev_fd-logev_classic)<1e-3
    fprintf('unit test 1: PASSED\n');
else
    warning('unit test 1: FAILED');
end


