% demo3_GPregress1Dcirc.m
%
% Fourier-domain GP regression with 100K 1D inputs on a circular domain.
%
% NOTE: UNFINISHED (Just started, July 2022).

if contains(pwd,'/demos')
    cd ..; setpaths; % move up a directory and call setpaths
end

% 1. Sample data from a GP

% Set hyperparameters for covariance function 
kprs.rho = 10; % marginal variance
kprs.len = 3; % length scale
signse = 2; % additive noise

% Sample points for collecting data
xrnge = [-50,50]; % support
nsamp = 1e5; % number of observations
xsamp = sort(rand(nsamp,1)*diff(xrnge)+xrnge(1)); % x locations

% Set parameters governing Fourier-domain representation
nSTD = 5; % number of length scales to extend circular interval
fdprs.circinterval = [xrnge(:,1),xrnge(:,2)+nSTD*kprs.len]'; % circular interval 
fdprs.condthresh = 1e8;  % threshold for cutting off small singular values
fdprs.minlen = 2; % minimum length scale to consider (set higher for increased speed)

% Sample function using Fourier representation
[cdiag,Bsamp] = Krbf_fourier(xsamp,kprs.len,kprs.rho,fdprs);
nfreq = length(cdiag);
fsamp = Bsamp*(sqrt(cdiag).*randn(nfreq,1));
ysamp = fsamp+signse*randn(nsamp,1);

% Make grid of x points (for visualizing function)
ngrid = 200; % number of points
xgrid = linspace(xrnge(1),xrnge(2),ngrid)'; % grid points
ftrue = interp1(xsamp,fsamp,xgrid,'linear'); % interpolated true func

% Plot samples
subplot(211);
plot(xsamp,ysamp, '.');
hold on; plot(xgrid,ftrue, 'k--', 'linewidth', 2); hold off;
title('samples');
xlabel('x'); ylabel('y');
box off;

%% Optimize log-evidence to infer hyperparameters

pp = fitGPregress_MaxEv(xsamp,ysamp,fdprs);

fprintf('\nLearned hyperparameters (+/- 2SD):\n----------------------------------\n');
fprintf('len:    %.2f  (+/- %.2f) [true=%2.2f]\n', pp.kprs.len, 2*pp.kprspost.CI(1),kprs.len);
fprintf('rho:    %.2f  (+/- %.2f) [true=%2.2f]\n', pp.kprs.rho, 2*pp.kprspost.CI(2),kprs.rho);
fprintf('nsevar: %.3f (+/- %.2f) [true=%2.3f]\n', pp.kprs.nsevar, 2*pp.kprspost.CI(3),signse.^2);


%% Compute predictive distribution on a grid

[fgrid,fci] = compPredDist_GP(xgrid,pp);

subplot(212);
p1 = plot(xgrid, fgrid, 'linewidth', 3);
feb20sd = 20*fci;
ebregionplot(xgrid,fgrid,feb20sd,feb20sd);
hold on; p2 = plot(xgrid,ftrue, 'k--', 'linewidth', 2); hold off;
title('posterior mean (+/- 20SD)');
box off;
legend([p2,p1], 'true function', 'posterior mean');
xlabel('x'); ylabel('f(x)');
