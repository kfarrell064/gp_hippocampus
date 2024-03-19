function check_transformHprs
% check_transformHyperprs
%
% Check that function transformHyperprs performs as expected

% Change to parent directory if needed
if contains(pwd,'/unittests')
    cd ..; setpaths; % move up a directory and call setpaths
    addpath unittests;
end

% 1. Check that transformHprs is invertible
fprintf('\nChecking transformHprs\n');
fprintf('======================\n')

kprs0 = 2*rand(3,1);  % generate random values for [len, rho, sig^2];
err1 = zeros(3,1);
for nd=1:3
    kprs1 = transformHprs(kprs0,nd,1); % transform to "trho" coords used for fitting
    kprs2 = transformHprs(kprs1,nd,-1); % transform back to standard "len", "rho", "sig^2".
    
    % Report difference
    err1(nd) = sum(abs(kprs2-kprs0));
    fprintf('Difference (nd=%d): %.6f\n',nd,err1(nd));
end

if  max(abs(err1)) < 1e-12
    fprintf('\nunit test 1: PASSED\n');
else
    warning('unit test 1: FAILED');
end

%% 2. Check Jacobian 
fprintf('\nChecking transformHprs Jacobian\n');
fprintf('===============================\n');

y1 = randn(3,1); % point in transformed parameter space
dx = 1e-5;

err2 = zeros(3,1);
for nd= 1:3
    [x1,~,J] = transformHprs(y1,nd,-1,eye(3));
    
    Jtest = zeros(3,3);
    for jj = 1:3
        dxvec = zeros(3,1);
        dxvec(jj) = dx;
        Jtest(jj,:) = (transformHprs(x1+dxvec,nd,1)-transformHprs(x1-dxvec,nd,1))'/(2*dx);
    end

    err2(nd) = max(abs(J(:)-Jtest(:))/max(J(:))); % compute error
    fprintf('Summed Abs Difference (nd=%d): %.6f\n',nd,err2(nd));
end

if  max(abs(err2)) < dx*max(J(:))
    fprintf('\nunit test 2: PASSED\n');
else
    warning('unit test 2: FAILED');
end


%% 3. Check that transformRho.m works the same as transformHpers

fprintf('\nChecking agreement of transformHprs and transformRho\n');
fprintf('====================================================\n')

kprs0 = 2*rand(3,1);  % generate random values for [len, rho, sig^2];
err3 = zeros(3,2);    % errors
for nd=1:3
    kprs1 = transformHprs(kprs0,nd,1); % transform to "trho" coords used for fitting
    trho = transformRho(kprs0(2),kprs0(1),nd,1);
    rho = transformRho(trho,kprs0(1),nd,-1);
    err3(nd,1) = exp(kprs1(2))-trho;
    err3(nd,2) = kprs0(2)-rho;
end

fprintf('Summed Abs Difference: %.6f\n',sum(abs(err3(:))));

if  max(abs(err3)) < 1e-12
    fprintf('\nunit test 3: PASSED\n');
else
    warning('unit test 3: FAILED');
end
