setpaths;  % set some necessary path variabls

% Load the data: 236 cells x 6261 positions
celldata = load("fentonlab.mat").output;

nsamp = 6261;  % the number of pos samples
tc = celldata.C(:,1:nsamp);  % matrix of single-cell calcium responses
pos = [1:6261];  % for now, assume each sampled index represents a unique position
ncells = 236; % number of cells in dataset

nd = 1; % dimensionality of inputs
xrnge = [min(pos),max(pos)];  % this is the minimum and maximum x location

% Set parameters governing Fourier-domain representation
fdprs.circinterval = [xrnge(1);xrnge(2)+100];
fdprs.condthresh = 1e8;  % threshold for cutting off small singular values
fdprs.minlen = 10; % minimum length scale to consider (set higher for increased speed)
% use min length scale of 20?

%% Specify the cell to analyze and plot its raw response data

for celln = 1:ncells % set the cell number
    
    % Process input data
    xx = pos; % x locations
    yy = tc(celln,:)'; % data for a single cell
    
    % Optionally blur out the xx and yy values by adding some Gaussian noise
    % xx = xx+randn(size(xx))*20;
    % yy = yy+randn(size(yy))*.2;
    
    % Plot samples
    subplot(211);
    plot(xx,yy, '.');
    title(sprintf('Cell #%d',celln));
    xlabel('Position (cm)'); ylabel('Response');
    set(gca,'xlim',xrnge);
    
    %% Fit GP model & hyperparamters to the data from a single cell
    
    % Optimize log-evidence to infer hyperparameters
    pp = fitGPregress_MaxEv(xx,yy,fdprs);
    
    fprintf('\nLearned hyperparameters (+/- 2SD):\n----------------------------------\n');
    fprintf('len:    %.2f  (+/- %.2f)\n', pp.kprs.len, 2*pp.kprspost.CI(1));
    fprintf('rho:    %.2f  (+/- %.2f)\n', pp.kprs.rho, 2*pp.kprspost.CI(2));
    fprintf('nsevar: %.3f (+/- %.2f)\n', pp.kprs.nsevar, 2*pp.kprspost.CI(3));
    
    
    %% Compute fitted GP mean on a grid
    
    % Make grid of x points (for visualizing function)
    ngrid = 400; % number of grid points
    dx = diff(xrnge)/ngrid; % grid spacing
    xgrid = (xrnge(1)+dx/2:dx:xrnge(2))'; % grid points
    
    % Set lengthscale to minimum length scale (if optimized value is below it)
    if (pp.kprs.len < fdprs.minlen)
        pp.kprs.len = fdprs.minlen;
    end
    
    % Compute the posterior mean on our grid
    [fgrid,fstd] = compPredDist_GP(xgrid,pp);
    fstd20 = 20*fstd;  % 20 times the standard deviation (for plotting)
    
    % Make plot of GP fit
    subplot(212);
    plot(xx,yy,'.'); hold on;
    p1 = plot(xgrid, fgrid, 'linewidth', 3);
    ebregionplot(xgrid,fgrid,fstd20,fstd20);
    set(gca,'xlim',xrnge);
    xlabel('Position (cm)'); ylabel('Response');
    box off; hold off;
    h = get(gca,'children');
    legend(h([3,1,2]), 'raw data', 'GP mean','+/- 20 SD error bars');
    
    pause; % pause so we can inspect each fit.

end