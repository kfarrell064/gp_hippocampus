% demo2_GPregress2D.m
%
% Full example of fourier-domain GP regression with 100K 2D inputs, including
% learning of hyperparameters via evidence maximization

if contains(pwd,'/demos')
    cd ..; setpaths; % move up a directory and call setpaths
end

% 1. Sample data from a 2D GP

% Set hyperparameters for covariance function 
kprs.rho = 100; % marginal variance
kprs.len = 2.5; % length scale
signse = 5; % additive noise

% Sample points for collecting data
xrnge = [-10,10]; % support
nsamp = 1e5; % number of observations
nd = 2; % dimensionality of inputs
xsamp = rand(nsamp,nd)*diff(xrnge)+xrnge(1); % x locations

% Optional: set parameters governing Fourier-domain representation
fdprs.circinterval = [xrnge(1);xrnge(2)+5*kprs.len]*ones(1,nd); % circular interval for function
fdprs.condthresh = 1e6;  % threshold for cutting off small singular values
fdprs.minlen = 2; % minimum length scale to consider (set higher for increased speed)

% Sample function using Fourier representation
[cdiag,Bsamp] = Krbf_fourier(xsamp,kprs.len,kprs.rho,fdprs);
nfreq = length(cdiag);
fsamp = Bsamp*(sqrt(cdiag).*randn(nfreq,1));
ysamp = fsamp+signse*randn(nsamp,1);

% Plot samples
subplot(121);
plot3(xsamp(:,1),xsamp(:,2),ysamp, '.');
title('samples');
xlabel('x'); ylabel('y');

%% Optimize log-evidence to infer hyperparameters

pp = fitGPregress_MaxEv(xsamp,ysamp,fdprs);

fprintf('\nLearned hyperparameters (+/- 2SD):\n----------------------------------\n');
fprintf('len:    %.2f  (+/- %.2f) [true=%2.1f]\n', pp.kprs.len, 2*pp.kprspost.CI(1),kprs.len);
fprintf('rho:    %.2f  (+/- %.2f) [true=%2.1f]\n', pp.kprs.rho, 2*pp.kprspost.CI(2),kprs.rho);
fprintf('nsevar: %.3f (+/- %.2f) [true=%2.3f]\n', pp.kprs.nsevar, 2*pp.kprspost.CI(3),signse.^2);


%% Compute predictive distribution on a grid

% Make grid of x points (for visualizing function)
ngrid = 40; % number of grid points
dx = diff(xrnge)/ngrid; % grid spacing
xgrid0 = (xrnge(1)+dx/2:dx:xrnge(2))'; % grid points for 1 axis
[xgrid,ygrid]=meshgrid(xgrid0); % x and y grids
xxgrid = [xgrid(:),ygrid(:)]; % stored as paired [x,y] values
nngrid = size(xxgrid,1);

[fgrid,fstd] = compPredDist_GP(xxgrid,pp);

subplot(122);
mesh(xgrid0,xgrid0,reshape(fgrid,ngrid,ngrid));
title('posterior mean');