function [] = CreateVarScript(num_nodes, n, depth)

% Uses the function 'CreateVarPartialConnected.m' to create a factor graph
% for each network with 'num_nodes'. The networks are located in the folder
% /GlobalVar/GenerateData/Networks.
%
% For the inputs 'n' and 'depth' see the function 
% 'CreateVarPartialConnected.m'

% =========================================================================
% Directories and filenames

% path to networks
dir_nets = '../../../../GlobalVar/GenerateData/Networks/Nets_';

% filename of networks file
file_nets = [dir_nets, num2str(num_nodes), '_nodes.mat'];

% Output filename
filename_output = ['GeneratedFGs/FG_Nets_', num2str(num_nodes), '_nodes.mat'];
% =========================================================================


% =========================================================================
% Code

load(file_nets);

Num_Nets = length(Networks);

FactorGraphs = cell(Num_Nets,1);

for i = 1 : Num_Nets
   
    P = Networks{i}.P;
    Neighbors = Networks{i}.Neighbors;
           
    FactorGraphs{i} = CreateVarPartialConnected(n, P, Neighbors, depth);            
end
% =========================================================================

save(filename_output, 'FactorGraphs');