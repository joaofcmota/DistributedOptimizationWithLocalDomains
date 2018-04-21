
% Filename for saving data
FILE_SAVE = 'Results/NetFlow_Barab100.mat';

% =========================================================================
% Execution parameters of the algorithms:

max_communications = 2000;
eps_opt = 1e-4;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}

% We will select the rho that is best for 1e-4 of accuracy
rhos_DADMM = [0.05 0.10 0.15];
% Same for the Lipschitz constant
Lips = [8400 8450 8500];
% =========================================================================


% =========================================================================
% Directories and filenames

addpath('../../GenerateData/Networks/MPC/ProcessedNetwork/');  % Networks
addpath('../../GenerateData/Data/NetworkFlows/Results/');      % Data

% Algorithms to be compared
addpath('../../');
addpath('../../KekatosADMM');
addpath('../../SpecialPurposeAlgs/BoydADMM');
addpath('../../SpecialPurposeAlgs/NesterovMethod');

% Solvers and auxiliary functions
addpath('../../UsageExamples/NetworkFlow/');
addpath('../../SpecialPurposeAlgs/NesterovMethod/NetworkFlow');

file_networks = 'Barabasi_P_100.mat'; 
file_data     = 'NFData_Barabasi_P_100.mat';
% =========================================================================


% =========================================================================
% Extract Data

% Load files
load(file_networks);
load(file_data);

pos_iter_errors = -log10(eps_opt);  %  Corresponding row in iter_for_errors

len_rhos = length(rhos_DADMM);
len_Lips = length(Lips);

Adj = Network.Adj;                             % Adjacency matrix
partition_colors = Network.Partition;          % Color partition of network
P = length(Adj);                               % Number of nodes
neighbors = Network.Neighbors;
Dp = Network.Degrees;

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );


% Create vars_prob

B = incidence_matrix;
num_edges = size(incidence_matrix,2);
d = flows;

components = cell(P,1);

for p = 1 : P
    components{p} = find(B(p,:) ~= 0);
end
% =========================================================================


% =========================================================================
% vars_prob for all algorithms

centers_primal = zeros(num_edges,1);
for edge = 1 : num_edges
    centers_primal(edge) = find(B(:,edge) == -1);
end

% For D-ADMM, Kekatos
vars_prob = struct('handler', {@NetFlowSolver}, ...
    'function_eval'     , {@Delays_f}, ...
    'gradient_eval'     , {@Delays_grad}, ...
    'projection'        , {@Delays_proj}, ...
    'NetFlowSolver_proj', {@NetFlowSolver_proj}, ...
    'spg_bb'            , {@SPG_BB}, ...
    'components', {components}, ...
    'B', {B}, ...
    'd', {d}, ...
    'capacities', {capacities}, ...
    'centers', {centers_primal}, ...
    'Dp', {Dp} ...
    );


% For NesterovMethod
components_Dual = cell(P,1);
for p = 1 : P
    components_Dual{p} = sort([neighbors{p}, p]);
end
centers_Dual = (1:P)';
x_estimate = zeros(num_edges,1);    % Will have the current estimate for x
   
vars_prob_Nest = struct('components', {components_Dual}, ...
    'centers', {centers_Dual}, ...
    'gradients', {@grads_NetFlowNest}, ...
    'proj_constraints', {@proj_NetFlowNest}, ...    
    'B', {B}, ...
    'd', {d}, ...
    'capacities', {capacities}, ...
    'x_estimate', {x_estimate}  ...
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
    ops = struct('rho', {rhos_DADMM(i_rhos)}, ...
        'max_iter', {max_communications}, ...        
        'x_opt', {solution}, ...
        'eps_opt', {eps_opt} ...
        );    
    
    fprintf('DADMM_Partial: start\n');    
    [X_DADMM_Partial, vars_prob_DADMM_Partial, ops_out] = ...
        DADMM_Partial(num_edges, vars_prob, vars_network, ops);    
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
        best_rhos_DADMM_Partial = rhos_DADMM(i_rhos);
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
    ops = struct('rho', {rhos_DADMM(i_rhos)}, ...
        'max_iter', {max_communications}, ...        
        'x_opt', {solution}, ...
        'eps_opt', {eps_opt} ...
        );    
    
    fprintf('KekatosADMM: start\n');    
    [X_KekatosADMM, vars_prob_KekatosADMM, ops_out] = ...
        KekatosADMM(num_edges, vars_prob, vars_network, ops);
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
        best_rhos_KekatosADMM = rhos_DADMM(i_rhos);
        iter_for_errors_KekatosADMM = iter_for_errors(pos_iter_errors,2);
    end
end
%******************************************************************


%******************************************************************
% BoydADMM

errors_BoydADMM = 0;
best_rhos_BoydADMM = 0;
iter_for_errors_BoydADMM = 0;

best_iter_BoydADMM = Inf;
for i_rhos = 1 : len_rhos
    
    % Optional input
    ops = struct('rho', {rhos_DADMM(i_rhos)}, ...
        'max_iter', {max_communications}, ...        
        'x_opt', {solution}, ...
        'eps_opt', {eps_opt} ...
        );    
    
    fprintf('BoydADMM: start\n');    
    [X_BoydADMM, vars_prob_BoydADMM, ops_out] = ...
        BoydADMM(num_edges, vars_prob, vars_network, ops);    
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
    
    if iter_for_errors(pos_iter_errors,2) < best_iter_BoydADMM
        best_iter_BoydADMM = iter_for_errors(pos_iter_errors,2);
        errors_BoydADMM = error_iterations;
        best_rhos_BoydADMM = rhos_DADMM(i_rhos);
        iter_for_errors_BoydADMM = iter_for_errors(pos_iter_errors,2);
    end
end
%******************************************************************


%******************************************************************
% Nesterov Method

errors_NesterovMethod = 0;
best_Lips_NesterovMethod = 0;
iter_for_errors_NesterovMethod = 0;

best_iter_Nesterov = Inf;
for i_Lip = 1 : len_Lips
    
    Lipschitz = Lips(i_Lip);
    
    % Optional input    
    ops = struct('max_iter', {max_communications}, ...
        'x_opt', {solution}, ...
        'error_fun', {@Error_gradNetFlow}, ...
        'eps_opt', {eps_opt} ...
        );    
    
    
    fprintf('Nesterov: start\n');    
    [X_Nest, vars_prob_Nest, ops_out] = NesterovMethod(P, Lipschitz, ...
        vars_prob_Nest, vars_network, ops);    
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
    
    if iter_for_errors(pos_iter_errors,2) < best_iter_Nesterov
        best_iter_Nesterov = iter_for_errors(pos_iter_errors,2);
        errors_NesterovMethod = error_iterations;
        best_Lips_NesterovMethod = Lips(i_Lip);
        iter_for_errors_NesterovMethod = iter_for_errors(pos_iter_errors,2);
    end
end
%******************************************************************


% =========================================================================


% =========================================================================
% Save data
save(FILE_SAVE, 'rhos_DADMM', 'Lips', 'max_communications', 'eps_opt', ...
    'errors_DADMM_Partial', 'best_rhos_DADMM_Partial', 'iter_for_errors_DADMM_Partial', ...
    'errors_KekatosADMM', 'best_rhos_KekatosADMM', 'iter_for_errors_KekatosADMM', ...
    'errors_BoydADMM', 'best_rhos_BoydADMM', 'iter_for_errors_BoydADMM', ...
    'errors_NesterovMethod', 'best_Lips_NesterovMethod', 'iter_for_errors_NesterovMethod' ...    
    );
% =========================================================================


figure(1);clf;
semilogy(errors_DADMM_Partial, 'b-');
hold on;
semilogy(errors_KekatosADMM, 'r-');
semilogy(errors_BoydADMM, 'k-');
semilogy(errors_NesterovMethod(2:end), 'c-');
legend('DADMMp', 'Kekatos', 'Boyd', 'Nesterov')
drawnow;







