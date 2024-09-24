% Load data
addpath(genpath(pwd)); warning off
datasetname = 'celldata_lvm.mat';  % name of dataset
load(datasetname);
xx = lvmdata.latentVariable;
yy = lvmdata.spikes;
ff = lvmdata.spikeRates;

% Get sizes and spike counts
[nt,nneur] = size(yy); % nt: number of time points; nneur: number of neurons
nf = size(xx,2); % number of latent dimensions

%% == 1. Compute baseline estimates ====

% Initialize the log of spike rates with the square root of spike counts.
ffmat = sqrt(yy);

% % Compute LLE
xlle = lle(ffmat,nf,20);
xllemat = align_xtrue(xlle,xx); % align the estimate with the true latent variable.

% Compute Poisson Linear Dynamic System (PLDS)
xplds = run_plds(yy,nf)';
xpldsmat = align_xtrue(xplds,xx); % align the estimate with the true latent variable.

%% == 2. Compute P-GPLVM ====

% Set up options
setopt.sigma2_init = 2; % initial noise variance
setopt.lr = 0.95; % learning rate
setopt.latentTYPE = 1; % kernel for the latent, 1. AR1, 2. SE
setopt.ffTYPE = 2; % kernel for the tuning curve, 1. AR1, 2. SE
setopt.initTYPE = 1; % initialize latent: 1. use PLDS init; 2. use random init; 3. true xx
setopt.la_flag = 1; % 1. no la; 2. standard la; 3. decoupled la
setopt.rhoxx = 10; % rho for Kxx
setopt.lenxx = 50; % len for Kxx
setopt.rhoff = 10; % rho for Kff
setopt.lenff = 50; % len for Kff
setopt.hypid = [1,2,3,4]; % 1. rho for Kxx; 2. len for Kxx; 3. rho for Kff; 4. len for Kff; 5. sigma2 (annealing it instead of optimizing it)
% setopt.xpldsmat = xllemat; % for plotting purpose
% setopt.xplds = xlle; % for initialization purpose
setopt.xpldsmat = xpldsmat; % for plotting purpose
setopt.xplds = xplds; % for initialization purpose
setopt.niter = 50; % number of iterations

% Compute P-GPLVM with Laplace Approximation
result_la = pgplvm_la(yy,xx,ff,setopt);

%% == 3. Plot latent variables and tuning curves ====

xxsampmat = align_xtrue(result_la.xxsamp,simdata.latentVariable);
subplot(211); plot(1:nt,xx,'b-',1:nt,xllemat,'r.-',1:nt,xpldsmat,'m.-',1:nt,xxsampmat,'k:','linewidth',2); legend('true x','LLE x','PLDS x','P-GPLVM x');
xlabel('time bin'); drawnow;

xgrid = gen_grid([min(xxsampmat(:,1)) max(xxsampmat(:,1))],50,nf); % x grid for plotting tuning curves
fftc = exp(get_tc(xxsampmat,result_la.ffmat,xgrid,result_la.rhoff,result_la.lenff));

neuronlist = 1:nneur;
ii = randperm(nneur);
neuronlist = neuronlist(ii(1:4));
for ii=1:4
    neuronid = neuronlist(ii);
    subplot(4,2,4+ii), cla
    hold on,
    plot(simdata.xgrid,simdata.tuningCurve(:,neuronid),'b-')
    plot(xgrid,fftc(:,neuronid),'r-')
    legend('true tc','estimated tc')
    title(['neuron ' num2str(neuronid)])
    hold off
end