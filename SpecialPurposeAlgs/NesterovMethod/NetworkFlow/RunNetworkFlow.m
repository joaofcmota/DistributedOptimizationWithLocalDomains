% Solves a network flow problem with NesterovMethod. The problem is
%
%        minimize     sum_{(i,j) in E}  phi_{ij}(x_{ij})
%        {x_{ij}}
%        subject to   B*x = d
%                     0 <= x <= c
%
% where the variable x_{ij} is defined over the edge ij of the network E, B
% is the node-arc incidence matrix, and d is the vector of sources (sum(d) 
% = 0). phi_{ij} is a function associated to the edge ij and it only 
% depends on the flow on that edge, x_{ij}.
%
% We will use for each function 
%
%    phi_{ij}(x_{ij}) = x_{ij}/(c_{ij} - x_{ij}) ,
%
% from a multi-commodity flow problem modeling delays at the links.
%
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

addpath('../../../GenerateData/Networks/MPC/ProcessedNetwork/');% Networks
addpath('../../../GenerateData/Data/NetworkFlows/Results/');    % Data
addpath('../');                                                 % Algorithm

%file_networks = 'Barabasi_P_2000.mat';
%file_data     = 'NFData_Barabasi_P_2000.mat'; 
file_networks = 'Barabasi_P_100.mat';
file_data     = 'NFData_Barabasi_P_100.mat'; 
% =========================================================================


% =========================================================================
% Parameters
Lipschitz = 15000;
max_iter = 5000;          % Maximum number of iterations
eps_opt = 1e-4;          % Tolerance
% =========================================================================


% =========================================================================
% Extract data and networks from files

% Load files
load(file_networks);
load(file_data);

Adj = Network.Adj;                             % Adjacency matrix
partition_colors = Network.Partition;          % Color partition of network
P = length(Adj);                               % Number of nodes
neighbors = Network.Neighbors;
Dp = Network.Degrees;

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors} ...    
    );


% Create vars_prob

B = incidence_matrix;
d = flows;

num_edges = size(incidence_matrix,2);

components = cell(P,1);

for p = 1 : P
    components{p} = sort([neighbors{p}, p]);
end

centers = (1:P)';

x_estimate = zeros(num_edges,1);    % Will have the current estimate for x
   
vars_prob = struct('components', {components}, ...
    'centers', {centers}, ...
    'gradients', {@grads_NetFlowNest}, ...
    'proj_constraints', {@proj_NetFlowNest}, ...    
    'B', {B}, ...
    'd', {d}, ...
    'capacities', {capacities}, ...
    'x_estimate', {x_estimate}  ...
    );
% =========================================================================


% =========================================================================
% Execute NesterovMethod

% Optional input
 ops = struct('max_iter', {max_iter}, ...
     'x_opt', {solution}, ...
     'error_fun', {@Error_gradNetFlow}, ...
     'eps_opt', {eps_opt} ...     
 );

[X, vars_prob, ops_out] = NesterovMethod(P, Lipschitz, vars_prob, ...
    vars_network, ops);
% =========================================================================


% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('Lipschitz = %d\n', Lipschitz);
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




