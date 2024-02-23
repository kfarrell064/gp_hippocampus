% setpaths
%
% Add necessary paths for fourier-domain GPs (fdgp)

if ~exist('Krbf','file')
    addpath kernelfuns
end

if ~exist('neglogEv_GPregress','file')
    addpath inference
end

if ~exist('demo1_GPregress1D','file')
    addpath demos
end


if ~exist('kronmulttrp','file')
    addpath utils
end

addpath tools_misc