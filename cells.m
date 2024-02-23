celldata = load("celldata.mat");
setpaths;
tc = celldata.tc(:,1:41756);
pos = celldata.y(1:41756);
[numcells, timepoints] = size(tc);

celln = 61;
% Set hyperparameters for covariance function 
kprs.rho = 100; % marginal variance
kprs.len = 2.5; % length scale
signse = 5; % additive noise

% Sample points for collecting data
xrnge = [0,400];
nsamp = timepoints; % number of observations
nd = 1; % dimensionality of inputs
xsamp = pos; % x locations

% Optional: set parameters governing Fourier-domain representation
fdprs.circinterval = [xrnge(1);xrnge(2)+5*kprs.len]*ones(1,nd);
fdprs.condthresh = 1e6;  % threshold for cutting off small singular values
fdprs.minlen = 2; % minimum length scale to consider (set higher for increased speed)

% Sample function using Fourier representation
% [cdiag,Bsamp] = Krbf_fourier(xsamp,kprs.len,kprs.rho,fdprs);
% nfreq = length(cdiag);
% fsamp = Bsamp*(sqrt(cdiag).*randn(nfreq,1));
ysamp = tc(celln,:)';

% Plot samples
subplot(121);
plot(xsamp,ysamp, '.');
title('Cell #61');
xlabel('Position (cm)'); ylabel('Response');

%% Optimize log-evidence to infer hyperparameters

pp = fitGPregress_MaxEv(xsamp,ysamp,fdprs);

fprintf('\nLearned hyperparameters (+/- 2SD):\n----------------------------------\n');
fprintf('len:    %.2f  (+/- %.2f) [true=%2.1f]\n', pp.kprs.len, 2*pp.kprspost.CI(1),kprs.len);
fprintf('rho:    %.2f  (+/- %.2f) [true=%2.1f]\n', pp.kprs.rho, 2*pp.kprspost.CI(2),kprs.rho);
fprintf('nsevar: %.3f (+/- %.2f) [true=%2.3f]\n', pp.kprs.nsevar, 2*pp.kprspost.CI(3),signse.^2);


%% Compute predictive distribution on a grid
% error bars
% information per cell
% git repository
% convolve ysamp to smooth

% Make grid of x points (for visualizing function)
ngrid = 400; % number of grid points
dx = diff(xrnge)/ngrid; % grid spacing
xgrid = (xrnge(1)+dx/2:dx:xrnge(2))'; % grid points
ngrid = size(xgrid);
[fgrid,fstd] = compPredDist_GP(xgrid,pp);

subplot(122);
% 
% plot(xgrid,fgrid);
% title('posterior mean');

p1 = plot(xgrid, fgrid, 'linewidth', 3);
fstd20 = 20*fstd;
ebregionplot(xgrid,fgrid,fstd20,fstd20);
title('posterior mean (+/- 20SD)');
box off;

trial_ms = size(pos,1);
avg_activity = sum(tc,2)/trial_ms;
bins = 50;
binlen = max(pos)/bins;
bins_ms = zeros(bins,1);

for i = 1:size(pos)
    if pos(i) >= 0
        binnum = ceil(pos(i) / binlen);
        if binnum == 0
            binnum = 1
        end
        bins_ms(binnum) = bins_ms(binnum)+1;
    end
end

% for cell = 1:numcells
%      gp = 
