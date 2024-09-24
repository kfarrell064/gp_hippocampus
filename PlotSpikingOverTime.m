setpaths;  % set some necessary path variabls

% Load the data
celldata = load("fentonlab.mat").output.S;

nsamp = 6261;  
ncells = size(celldata,1)
celln = 1
%% Specify the cell to analyze and plot its raw response data

if celln <= ncells
    % Process input data
    xx = [1:1:nsamp]; % x axis
    yy = celldata(celln,:)'; % data for a single cell
    
    % Plot samples
    subplot(211);
    plot(xx,yy, '.');
    xlim([0 nsamp])
    hold off
    title(sprintf('Cell #%d',celln));
    xlabel('Timestamp'); ylabel('Response');
    celln = celln + 1


end