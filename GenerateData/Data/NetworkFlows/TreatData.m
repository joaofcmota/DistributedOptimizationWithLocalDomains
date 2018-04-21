% Treats the data for the Network Flow problem:
%
%      minimize     \sum_{(i,j) in edges}  x_ij/(c_ij - x_ij)      (1)
%     x = {x_ij}
%      subject to   B*x = d
%                   0 <= x <= c,
%
% where B is the node-arc incidence matrix of the graph, c is the vector of
% capacities, and d is the vector of flows.
%
% B, c, d, the network, and the solution of (1) were already computed with 
% the Sage script
%
%    GenerateNetworkFlows.sage
%
% and the results are stored in the folder Results.
%
% Here, we just treat that data.

% =========================================================================
% Filenames (NOTE: THE INPUT FILE WILL BE REWRITTEN)

FILENAME_DATA_INPUT  = 'Results/NFData_Barabasi_P_6.mat';
FILENAME_DATA_OUTPUT = FILENAME_DATA_INPUT;
% =========================================================================


load(FILENAME_DATA_INPUT);

capacities = double(capacities);
flows = double(flows);
incidence_matrix = sparse(double(incidence_matrix));
number_of_commodities = length(solution_mcfp);
node_out_node_in_intensity = double(node_out_node_in_intensity);

[num_nodes, num_edges] = size(incidence_matrix);


save(FILENAME_DATA_OUTPUT, 'capacities', 'flows', 'incidence_matrix', ...
    'num_nodes', 'num_edges', 'number_of_commodities', ...
    'node_out_node_in_intensity', 'solution_mcfp', 'solution')