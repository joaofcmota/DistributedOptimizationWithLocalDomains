
% Filename for saving data
FILE_SAVE = 'Results/Barab100_NonStar_Stable.mat';

% =========================================================================
% Execution parameters of the algorithms:

max_communications = 1000;
eps_opt = 1e-4;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}

% We will select the rho that is best for 1e-4 of accuracy
rhos = 40;             % Precision +/-5
% =========================================================================


% =========================================================================
% Directories

addpath('../../GenerateData/Networks/MPC/ProcessedNetwork/');   % Networks
addpath('../../GenerateData/Data/MPC/Data/');                   % Data

% Algorithms to be compared
addpath('../../');
addpath('../../KekatosADMM');

% Solvers and auxiliary functions
addpath('../../UsageExamples/MPC_PowerGrid');

file_networks = 'Barabasi_P_100.mat'; 
file_data     = 'Barab_NonStar_P_100_STABLE.mat';
% =========================================================================


% =========================================================================
% Extract Data

% Load files
load(file_networks);
load(file_data);

pos_iter_errors = -log10(eps_opt);  %  Corresponding row in iter_for_errors

len_rhos = length(rhos);

% Network data
P = Network.P;
Adj = Network.Adj;
Num_Colors = Network.Num_Colors;
Partition = Network.Partition;
Neighbors = Network.Neighbors;
Degrees = Network.Degrees;
% =========================================================================


% =========================================================================
% vars_prob and vars_network

vars_prob = struct('handler', {@SolverMPC}, ...
    'components', {components}, ...
    'E_p', {E_p}, ...
    'w_p', {w_p} ...
    );

vars_network = struct('P', {P}, ...
    'neighbors', {Neighbors}, ...
    'partition_colors', {Partition} ...
    );
% =========================================================================


% =========================================================================
% Execute Algorithms


%******************************************************************
% DADMM_Partial
errors_DADMM_Partial = 0;
best_rhos_DADMM_Partial = 0;
iter_for_errors_DADMM_Partial = 0;

best_iter_DADMM_Partial = Inf;
for i_rhos = 1 : len_rhos
    
    % Optional input
    ops = struct('rho', {rhos(i_rhos)}, ...
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
    
    if iter_for_errors(pos_iter_errors,2) < best_iter_DADMM_Partial
        best_iter_DADMM_Partial = iter_for_errors(pos_iter_errors,2);
        errors_DADMM_Partial = error_iterations;
        best_rhos_DADMM_Partial = rhos(i_rhos);
        iter_for_errors_DADMM_Partial = iter_for_errors(pos_iter_errors,2);
    end
end
%******************************************************************


%******************************************************************
% Kekatos ADMM
errors_KekatosADMM = 0;
best_rhos_KekatosADMM = 0;
iter_for_errors_KekatosADMM = 0;

best_iter_KekatosADMM = Inf;
for i_rhos = 1 : len_rhos
    
    % Optional input
    ops = struct('rho', {rhos(i_rhos)}, ...
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
    
    if iter_for_errors(pos_iter_errors,2) < best_iter_KekatosADMM
        best_iter_KekatosADMM = iter_for_errors(pos_iter_errors,2);
        errors_KekatosADMM = error_iterations;
        best_rhos_KekatosADMM = rhos(i_rhos);
        iter_for_errors_KekatosADMM = iter_for_errors(pos_iter_errors,2);
    end
end
%******************************************************************


% =========================================================================


% =========================================================================
% Save data
save(FILE_SAVE, 'rhos', 'max_communications', 'eps_opt', ...
    'errors_DADMM_Partial', 'best_rhos_DADMM_Partial', 'iter_for_errors_DADMM_Partial', ...
    'errors_KekatosADMM', 'best_rhos_KekatosADMM', 'iter_for_errors_KekatosADMM' ...
    );
% =========================================================================


figure(1);clf;
semilogy(errors_DADMM_Partial, 'b-');
hold on;
semilogy(errors_KekatosADMM, 'r-');
legend('DADMMp', 'Kekatos')
drawnow;



