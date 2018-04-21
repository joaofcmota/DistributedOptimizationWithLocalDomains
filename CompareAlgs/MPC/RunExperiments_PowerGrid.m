
% Filename for saving data
FILE_SAVE = 'Results/PowerGrid_Stable.mat';

% =========================================================================
% Execution parameters of the algorithms:

max_communications = 1000;
eps_opt = 1e-4;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}

% We will select the rho that is best for 1e-4 of accuracy
%rhos = [20, 25, 30, 35];                           
rho_DADMM = 25;                  % All rhos selected with precision +/-5
rho_Kekatos = 30;
rho_Boyd = 25;
% =========================================================================


% =========================================================================
% Directories and filenames

addpath('../../GenerateData/Networks/MPC/ProcessedNetwork/');    % Networks
addpath('../../GenerateData/Data/MPC/Data/');                    % Data

% Algorithms to be compared
addpath('../../');
addpath('../../KekatosADMM');
addpath('../../SpecialPurposeAlgs/BoydADMM');
addpath('../../SpecialPurposeAlgs/NesterovMethod');

% Solvers and auxiliary functions
addpath('../../UsageExamples/MPC_PowerGrid');
addpath('../../SpecialPurposeAlgs/NesterovMethod/MPC');

file_networks = 'Network_PowerGrid.mat'; 
file_data     = 'MPC_PGrid_P_4941_STABLE.mat';
% =========================================================================

% =========================================================================
% Extract Data

% Load files
load(file_networks);
load(file_data);

pos_iter_errors = -log10(eps_opt);  %  Corresponding row in iter_for_errors

% Network data
P = Network.P;
Adj = Network.Adj;
Num_Colors = Network.Num_Colors;
Partition = Network.Partition;
Neighbors = Network.Neighbors;
Degrees = Network.Degrees;

centers = kron(1:P,ones(1,T*m_p));
% =========================================================================


% =========================================================================
% vars_prob for all algorithms

% **********************************************************
% D-ADMM, Kekatos, and BoydADMM

vars_prob = struct('handler', {@SolverMPC}, ...
    'components', {components}, ...
    'centers', {centers}, ...                    % For Boyd's algorithm
    'E_p', {E_p}, ...
    'w_p', {w_p} ...
    );


vars_network = struct('P', {P}, ...
    'neighbors', {Neighbors}, ...
    'partition_colors', {Partition} ...
    );
% **********************************************************

% **********************************************************
% Gradient and Nesterov Methods

vars_prob_Grad = struct('gradients', {@GradientsMPC}, ...
    'proj_constraints', {@projConsMPC}, ...
    'components', {components}, ...
    'centers', {centers}, ...
    'E_p', {E_p}, ...
    'w_p', {w_p} ...
    );

% =========================================================================


% =========================================================================
% Execute Algorithms


%******************************************************************
% DADMM_Partial

% Optional input
ops = struct('rho', {rho_DADMM}, ...
    'max_iter', {max_communications}, ...
    'x_opt', {solution}, ...
    'eps_opt', {eps_opt} ...
    );

fprintf('DADMM_Partial: start\n');
[X_DADMM_Partial, vars_prob_DADMM_Partial, ops_out] = ...
    DADMM_Partial(length(solution), vars_prob, vars_network, ops);
fprintf('DADMM_Partial: finish\n');

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
for i_g = 1 : 6
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end

errors_DADMM_Partial = error_iterations;
best_rhos_DADMM_Partial = rho_DADMM;
iter_for_errors_DADMM_Partial = iter_for_errors(pos_iter_errors,2);
%******************************************************************


%******************************************************************
% Kekatos ADMM

% Optional input
ops = struct('rho', {rho_Kekatos}, ...
    'max_iter', {max_communications}, ...
    'x_opt', {solution}, ...
    'eps_opt', {eps_opt} ...
    );

fprintf('KekatosADMM: start\n');
[X_KekatosADMM, vars_prob_KekatosADMM, ops_out] = ...
    KekatosADMM(length(solution), vars_prob, vars_network, ops);
fprintf('KekatosADMM: finish\n');

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
for i_g = 1 : 6
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end

errors_KekatosADMM = error_iterations;
best_rhos_KekatosADMM = rho_Kekatos;
iter_for_errors_KekatosADMM = iter_for_errors(pos_iter_errors,2);
%******************************************************************


%******************************************************************
% BoydADMM
    
% Optional input
ops = struct('rho', {rho_Boyd}, ...
    'max_iter', {max_communications}, ...
    'x_opt', {solution}, ...
    'eps_opt', {eps_opt} ...
    );

fprintf('BoydADMM: start\n');
[X_BoydADMM, vars_prob_BoydADMM, ops_out] = ...
    BoydADMM(length(solution), vars_prob, vars_network, ops);
fprintf('BoydADMM: finish\n');

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
for i_g = 1 : 6
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end

errors_BoydADMM = error_iterations;
best_rhos_BoydADMM = rho_Boyd;
iter_for_errors_BoydADMM = iter_for_errors(pos_iter_errors,2);
%******************************************************************


%******************************************************************
% Nesterov Method

% Optional input
ops = struct('max_iter', {max_communications}, ...
    'x_opt', {solution}, ...
    'eps_opt', {eps_opt} ...
    );

fprintf('Nesterov: start\n');
[X_Nest, vars_prob_Grad, ops_out] = NesterovMethod(length(solution), Lipschitz, ...
    vars_prob_Grad, vars_network, ops);
fprintf('Nesterov: finish\n');

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
for i_g = 1 : 6
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end

errors_NesterovMethod = error_iterations;
iter_for_errors_NesterovMethod = iter_for_errors(pos_iter_errors,2);
%******************************************************************

% =========================================================================



% =========================================================================
% Save data
save(FILE_SAVE, 'rho_DADMM', 'rho_Kekatos', 'rho_Boyd', 'max_communications', 'eps_opt', ...
    'errors_DADMM_Partial', 'best_rhos_DADMM_Partial', 'iter_for_errors_DADMM_Partial', ...
    'errors_KekatosADMM', 'best_rhos_KekatosADMM', 'iter_for_errors_KekatosADMM', ...
    'errors_BoydADMM', 'best_rhos_BoydADMM', 'iter_for_errors_BoydADMM', ...
    'errors_NesterovMethod', 'iter_for_errors_NesterovMethod' ...
    );
% =========================================================================


figure(1);clf;
semilogy(errors_DADMM_Partial, 'b-');
hold on;
semilogy(errors_KekatosADMM, 'r-');
semilogy(errors_BoydADMM, 'k-');
semilogy(errors_NesterovMethod, 'c-');
legend('DADMMp', 'Kekatos', 'Boyd', 'Nesterov')
drawnow;









