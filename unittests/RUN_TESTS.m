% RUN_TESTS

% Runs all unit tests contained in this directory
% (NOTE: must be called from the unittests directory!)
cd ~/Dropbox/Docs/gitcode/fdGP_regression/unittests/  % EDIT THIS LINE

% Get list of files in current directory
D=dir;

% Set to 1 to pause between unit tests (eg, to examine plots)
PAUSEON = 1;

% Run all scripts except this one
for jj=3:length(D)
    nm = D(jj).name;
    if ~strcmp(nm(1),'.') && ~strcmp(nm(1:3),'RUN')
        fprintf('\n\nRUNNING UNIT TEST: %s\n', nm);
        eval(nm(1:end-2))
        if PAUSEON
            fprintf('press any key to run next test...\n');
            pause;
        end
    end
end