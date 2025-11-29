clear all
close all
clc

%%%%%%%%%%%%%%
% PARAMETERS %
%%%%%%%%%%%%%%
N               = 64;                       % --- The image size is NxN
maxIter         = 200;                      % --- Maximum number of CG iterations
tol             = 1e-6;                     % --- CG convergence tolerance

%%%%%%%%%%%%%%%%%%%%%
% DATA CONSTRUCTION %
%%%%%%%%%%%%%%%%%%%%%
phant           = phantom(N);               % --- Reference phantom
theta           = 0 : 2 : 177;              % --- Projection angles
[sinogram, xp]  = radon(phant, theta);
m               = numel(sinogram);          % --- Number of data
n               = N * N;                    % --- Number of unknowns
b               = sinogram(:);              % --- Known term of the linear system to be solved

%%%%%%%%%%%%%%%%%%%
% SYSTEM MATRIX A %
%%%%%%%%%%%%%%%%%%%
fprintf('Building system matrix A (m=%d, n=%d)...\n', m, n);
A = sparse(m, n);
for j = 1:n
    % Create canonical basis image
    img = zeros(N, N);
    img(j) = 1;  % j is linear index
    
    proj = radon(img, theta);
    A(:, j) = proj(:);
    
    if mod(j, 500) == 0, fprintf('.'); end
end
fprintf('\n...done\n');

%%%%%%%%%%%%%%%%%%%
% NORMAL EQUATION %
%%%%%%%%%%%%%%%%%%%
ATA             = A' * A;
ATb             = A' * b;

fprintf('Running unpreconditioned CG...\n');
[x_cg, flag1, relres1, iter1, resvec1]  = pcg(ATA, ATb, tol, maxIter);
x_cg_img                                = reshape(x_cg, N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% JACOBI PRECONDITIONING %
%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing Jacobi preconditioner...\n');
diagA                   = diag(ATA);
diag_min                = max(diagA) * 1e-10;  
diagA(diagA < diag_min) = diag_min;
M_jac                   = spdiags(1 ./ diagA, 0, n, n);

fprintf('Running Jacobi-preconditioned CG...\n');
[x_jac, flag2, relres2, iter2, resvec2] = pcg(ATA, ATb, tol, maxIter, M_jac);
x_jac_img                               = reshape(x_jac, N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INCOMPLETE CHOLESKY PRECONDITIONING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing incomplete Cholesky preconditioner...\n');
opts.droptol            = 1e-2;       
opts.type               = 'nofill';

% Method 1: Try standard ichol first
success                 = false;
tau_values              = [1e-4, 1e-3, 1e-2, 1e-1];  % Try different regularization levels
for i = 1:length(tau_values)
    tau = tau_values(i) * mean(diagA);
    try
        L = ichol(ATA + tau*speye(n), opts);
        fprintf('Ichol successful with tau = %.3e\n', full(tau));
        success = true;
        break;
    catch
        fprintf('  Failed with tau = %.3e\n', full(tau));
    end
end

% Method 2: If standard ichol fails, use modified ichol for ill-conditioned matrices
if ~success
    fprintf('Trying modified ichol approach...\n');
    try
        % Use a different dropping strategy
        opts.type = 'ict';
        opts.droptol = 0.1;  % More aggressive dropping
        L = ichol(ATA + 0.1 * speye(n), opts);  % Strong regularization
        fprintf('Modified ichol successful\n');
        success = true;
    catch
        fprintf('All ichol attempts failed, using Jacobi preconditioner\n');
        L = spdiags(sqrt(1 ./ diagA), 0, n, n);  % Fallback to Jacobi
    end
end

fprintf('Running preconditioned pcg...\n');
[x_ichol, flag3, relres3, iter3, resvec3]   = pcg(ATA, ATb, tol, maxIter, L, L');
x_ichol_img                                 = reshape(x_ichol, N, N);

% --- Plot reconstructions
figure('Name','Reconstructions','Position',[100 100 1000 350]);
subplot(1,4,1); imagesc(phant); axis image off; title('True phantom');
subplot(1,4,2); imagesc(x_cg_img); axis image off; title(sprintf('CG (%d it)', iter1));
subplot(1,4,3); imagesc(x_jac_img); axis image off; title(sprintf('Jacobi PCG (%d it)', iter2));
subplot(1,4,4); imagesc(x_ichol_img); axis image off; title(sprintf('Ichol PCG (%d it)', iter3));
colormap gray;

% --- Plot residual histories (normalized)
figure('Name','Residual histories');
semilogy(resvec1 / resvec1(1),'LineWidth',2); hold on;
semilogy(resvec2 / resvec2(1),'LineWidth',2);
semilogy(resvec3 / resvec3(1),'LineWidth',2);
xlabel('Iteration'); ylabel('Relative residual');
legend('CG','Jacobi PCG','Ichol PCG'); grid on;
title('Convergence Comparison');
