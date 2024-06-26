setpaths;  % set some necessary path variabls

% Load the data
celldata = load("celldata.mat");

nsamp = 41756;  % the number of samples to keep from raw data (why?)
tc = celldata.tc(:,1:nsamp);  % matrix of single-cell calcium responses
pos = celldata.y(1:nsamp);  % animal's sampled 1D spatial position over time
nd = 1; % dimensionality of inputs

%% Fit GP model & hyperparamters to the data from a single cell

celln = 1; % set the cell number

% Process input data
xrnge = [min(pos),max(pos)];
xx = pos; % x locations
yy = tc(celln,:)'; % data for a single cell

% Set parameters governing Fourier-domain representation
fdprs.circinterval = [xrnge(1);xrnge(2)+50];
fdprs.condthresh = 1e8;  % threshold for cutting off small singular values
fdprs.minlen = 50; % minimum length scale to consider (set higher for increased speed)

% Plot samples
subplot(221);
plot(xx,yy, '.');
title('Cell #61');
xlabel('Position (cm)'); ylabel('Response');

%% Optimize log-evidence to infer hyperparameters

pp = fitGPregress_MaxEv(xx,yy,fdprs);

fprintf('\nLearned hyperparameters (+/- 2SD):\n----------------------------------\n');
fprintf('len:    %.2f  (+/- %.2f)\n', pp.kprs.len, 2*pp.kprspost.CI(1));
fprintf('rho:    %.2f  (+/- %.2f)\n', pp.kprs.rho, 2*pp.kprspost.CI(2));
fprintf('nsevar: %.3f (+/- %.2f)\n', pp.kprs.nsevar, 2*pp.kprspost.CI(3));


%% Compute predictive distribution on a grid

% Make grid of x points (for visualizing function)
ngrid = 400; % number of grid points
dx = diff(xrnge)/ngrid; % grid spacing
xgrid = (xrnge(1)+dx/2:dx:xrnge(2))'; % grid points
ngrid = size(xgrid);

if (pp.kprs.len < fdprs.minlen)
    pp.kprs.len = fdprs.minlen;
end
[fgrid,fstd] = compPredDist_GP(xgrid,pp);

subplot(122);

p1 = plot(xgrid, fgrid, 'linewidth', 3);
fstd20 = 20*fstd;
ebregionplot(xgrid,fgrid,fstd20,fstd20);
title('posterior mean (+/- 20SD)');
box off;

%% Compute spatial information I without GPR - don't binarize
bins = 20;
total_ms = size(pos,1);
avg_activity = sum(tc,2)/total_ms;
bin_activity = zeros(size(tc,1), bins);
binlen = max(pos)/bins;
bins_ms = zeros(1, bins);
% on_threshold = 3.5*std(tc,0,2);

for i = 1:size(pos)
    if pos(i) >= 0
        binnum = ceil(pos(i) / binlen);
        if binnum == 0
            binnum = 1;
        end
        bins_ms(1,binnum) = bins_ms(1,binnum)+1;
        bin_activity(:,binnum) = bin_activity(:,binnum) + tc(:,i);   
    end
end

bin_avgs = bin_activity./bins_ms;
info_original = zeros(size(tc,1),1);

for i = 1:size(info_original)
    r = bin_avgs(i,:)/avg_activity(i);
    % sometimes a cell has 0 activity in a given bin, can't take log
    r(r==0) = 1;
    info_original(i) = sum(r.*log2(r));
end
figure
scatter(1:1325,info_original, 10, "filled");
title('Information per cell (unshuffled, no GPR)');
xlabel("Cell number")
ylabel("Information (bits)")

%% Non GP shuffling test
% plot information per cell on the x, y axes to compare them
n_shuffles = 500;
shuffled_distributions = zeros(size(tc,1),n_shuffles);

for s=1:n_shuffles
    bins_ms_shuff = zeros(1, bins);
    bin_activity_shuff = zeros(size(tc,1), bins);
    tc_shuff = circshift(tc, randi([3,size(pos,1)]), 2);
    for x = 1:size(pos)
        if pos(x) >= 0
            binnum = ceil(pos(x) / binlen);
            if binnum == 0
                binnum = 1;
            end
            bins_ms_shuff(1,binnum) = bins_ms_shuff(1,binnum)+1;
            bin_activity_shuff(:,binnum) = ...
                bin_activity_shuff(:,binnum) + tc_shuff(:,x);  
        end
    end
    bin_avgs_shuff = bin_activity_shuff./bins_ms_shuff;
    info_shuff = zeros(size(tc_shuff,1),1);

    for i = 1:size(info_shuff)
        r = bin_avgs_shuff(i,:)/avg_activity(i);
        r(r==0) = 1;
        info_shuff(i) = sum(r.*log2(r));
    end
    shuffled_distributions(:,s) = info_shuff;
end
%% Plot non GP shuffled information
sdists = load("shuffled_nogp.mat").shuffled_distributions;
figure
scatter((1:500),sdists(celln,:), 10, "red", "filled");
yline(info_original(celln))
title("Information for Cell " +celln+" over 500 shuffles (no GPR)");
xlabel("Shuffle number");
ylabel("Information (bits)");
%% Compute spatial information I with GPR  
% fit before looping, store gridded posterior/hyperparameters
% average each larger spatial bin 
% turn zeros into 1 / nansum
% if there is high information try shuffling by trial

celldata = load("celldata.mat");
setpaths;
tc = celldata.tc(:,1:41756);
pos = celldata.y(1:41756);
[numcells, timepoints] = size(tc);
cell_gridpts = zeros(size(tc,1), 400);

for celln = 1:size(tc,1)
    % Set hyperparameters for covariance function 
    kprs.rho = 100; % marginal variance
    kprs.len = 2.5; % length scale
    signse = 5; % additive noise
    
    % Sample points for collecting data
    xrnge = [0,400];
    nsamp = timepoints; % number of observations
    nd = 1; % dimensionality of inputs
    xx = pos; % x locations
    
    % Optional: set parameters governing Fourier-domain representation
    fdprs.circinterval = [xrnge(1);xrnge(2)+5*kprs.len]*ones(1,nd);
    fdprs.condthresh = 1e6;  % threshold for cutting off small singular values
    fdprs.minlen = 2; % minimum length scale to consider (set higher for increased speed)
    
    yy = tc(celln,:)';
    
    pp = fitGPregress_MaxEv(xx,yy,fdprs);
    
    
    % Make grid of x points (for visualizing function)
    ngrid = 400; % number of grid points
    dx = diff(xrnge)/ngrid; % grid spacing
    xgrid = (xrnge(1)+dx/2:dx:xrnge(2))'; % grid points
    ngrid = size(xgrid);
    [fgrid,fstd] = compPredDist_GP(xgrid,pp);
    cell_gridpts(celln,:) = cell_gridpts(celln) + fgrid;
    
end
%% GP shuffling test

%% Notes

% ignore pksdis
% compare no GP 20 bin representation to information captured with the GP
% plot on the same axes
% error bars
% information per cell
% git repository
% convolve ysamp to smooth
% use bins without binarizing
% make function to compute information per cell
