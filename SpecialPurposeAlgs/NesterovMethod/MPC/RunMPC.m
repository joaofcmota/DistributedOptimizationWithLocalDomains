% Solves the MPC problem with NesterovMethod. Given a time-horizon T, it 
% solves
%
%     minimize     sum_{t=0}^{T-1} [x(t)'*Q*x(t) + u(t)'*R*u(t)]
%                                                          + x(T)'*Q_f*x(T)
%
%     subject to   x_p(t+1) = A_pp x_p(t) + sum_{j \in Omega} B_pj*u_j[t], 
%                                                              p=1,...,P,
%
% where Q and R are diagonal matrices and each state has a dynamics that
% only depends on its own state but it is controlled by the actuators to
% which the sensor is connected.
%
% See the file '../../../GenerateData/Data/MPC/Readme/MPC.pdf' 
% for a more comprehensive explanation of the problem.
% -------------------------------------------------------------------------
% Distributed Optimization With Local Domains: Applications in MPC and
% Network Flows
% Copyright (C) 2014 Joao Mota
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
%
% =========================================================================
% Please send any questions, comments, or bug reports to j.mota@ucl.ac.uk
% =========================================================================


% =========================================================================
% Directories and filenames

addpath('../../../GenerateData/Networks/MPC/ProcessedNetwork/'); % Networks
addpath('../../../GenerateData/Data/MPC/Data/');                 % Data
addpath('../');                                              % Algorithm

file_networks = 'Barabasi_P_100.mat';
file_data     = 'Barab_P_100_STABLE.mat'; 
%file_networks = 'Barabasi_P_2000.mat';
%file_data     = 'Barab_P_2000_STABLE.mat'; 
%file_networks = 'Network_PowerGrid.mat';
%file_data     = 'MPCDataPowerGrid.mat'; 
% =========================================================================


% =========================================================================
% Parameters

max_iter = 1000;         % Maximum number of iterations
eps_opt = 1e-3;          % Tolerance
% =========================================================================


% =========================================================================
% Extract data and networks from files

% Load files
load(file_networks);
load(file_data);


% Network data
P = Network.P;
Adj = Network.Adj;
Num_Colors = Network.Num_Colors;
Partition = Network.Partition;
Neighbors = Network.Neighbors;
Degrees = Network.Degrees;

centers = kron(1:P,ones(1,T*m_p));

vars_prob = struct('gradients', {@GradientsMPC}, ...
    'proj_constraints', {@projConsMPC}, ...
    'components', {components}, ...
    'centers', {centers}, ...
    'E_p', {E_p}, ...
    'w_p', {w_p} ...
    );

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {Neighbors} ...
    );
% =========================================================================


% =========================================================================
% Execute NesterovMethod

% Optional input
 ops = struct('max_iter', {max_iter}, ...
     'x_opt', {solution}, ...
     'eps_opt', {eps_opt} ...     
 );

[X, vars_prob, ops_out] = NesterovMethod(length(solution), Lipschitz, ...
    vars_prob, vars_network, ops);
% =========================================================================


% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
fprintf('iter_for_errors = \n');
num_rows = size(iter_for_errors, 1);
for i_g = 1 : num_rows
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end

figure(1);clf;
semilogy(1:iterations,error_iterations(1:iterations), 'b');
title('error\_{iterations}');
% =========================================================================





