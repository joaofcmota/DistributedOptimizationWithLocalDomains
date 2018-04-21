function [] = generate_series_bipartite()

% =========================================================================
% Generates a series of bipartite networks and factor graphs

% Number of nodes in each set
C = [2  5   9  10  20  25  25  40];
P = [5  10 10  20  25  28  40  50];

FILENAME = 'Bipartite_Nets.mat';     % Where the networks will be stored

prob = 0.25;   % connection probability in Erdos-Renyi model

MAX_TRIALS = 100;   % Will try to generate a connected network up to 
                    % MAX_TRIALS                   
                   
python_filename = 'gen_bipartite_erdos.py';  % Python script that generates
                                             % network with Networkx
% =========================================================================

filename_temp = 'single_netw_tmp.txt';


num_networks = length(C);

if length(P) ~= num_networks
    error('C and P do not have the same number of elements');
end

Networks = cell(num_networks,1);
FactorGraphs = cell(num_networks,1);

for net = 1 : num_networks
    
    for trial = 1 : MAX_TRIALS
        
        system(['./', python_filename, ' ', num2str(C(net)), ' ', ...
            num2str(P(net)), ' ', num2str(prob), ' ', filename_temp ]);
        
        [conn, Adj, Diameter, colors] = read_Adjacency_and_connectivity(filename_temp);
        system(['rm', ' ', filename_temp]);
        
        if conn == 1
            break;
        end        
    end
    if conn == 0
        error('Was never able to generate a connected network for C x P = %d x %d', C(net),P(net));
    end
    
    part1 = find(colors);
    part2 = find(~colors);
    
    % We make the central nodes always less than the peripheral nodes
    if length(part1) < length(part2)
        central_nodes = part1;
        periphe_nodes = part2;
    else
        central_nodes = part2;
        periphe_nodes = part1;
    end
    
    num_nodes = C(net) + P(net);
    
    %  ====================================================================
    % Compute network information
    
    Parameters = prob;
    Num_Colors = 2;    
    Partition = {part1' , part2'};
    Neighbors = cell(num_nodes,1);
    Degrees = zeros(num_nodes,1);
    for p = 1 : num_nodes
        Neighbors{p} = find(Adj(p,:));
        Degrees(p) = length(Neighbors{p});
    end        
    Max_Degree = max(Degrees);
        
    Network = struct('P', {num_nodes}, ...
        'central_nodes', {central_nodes}, ...
        'periphe_nodes', {periphe_nodes}, ...
        'Adj', {Adj}, ...
        'Type', 'Erdos-Renyi_bipartite', ...
        'Parameters', {Parameters}, ...
        'Num_Colors', {Num_Colors}, ...
        'Partition', {Partition}, ...
        'Diameter', {Diameter}, ...
        'Neighbors', {Neighbors}, ...
        'Degrees', {Degrees}, ...
        'Max_Degree', {Max_Degree} ...
        );
    %  ====================================================================
    
    %  ====================================================================
    % Compute factor graph
    
    components = cell(num_nodes, 1);
    
    % Central nodes have only an associated variable component
    for ind_c = 1 : length(central_nodes)
        c = central_nodes(ind_c);
        components{c} = ind_c;
    end
        
    for ind_p = 1 : length(periphe_nodes)
        
        p = periphe_nodes(ind_p);
        neighbs_p = Neighbors{p};   % These should be central nodes
        Dp = Degrees(p);
        
        components_p = zeros(Dp,1);
        for i_Dp = 1 : Dp
            c = neighbs_p(i_Dp);
            components_p(i_Dp) = components{c};
        end
        
        components{p} = components_p;
    end
    
    FactorGraph = struct('components', {components}, ...
        'P', {num_nodes}, ...
        'neighbors', {Neighbors} ...
    );
    %  ====================================================================
    
    Networks{net} = Network;
    FactorGraphs{net} = FactorGraph;
    
    save(FILENAME, 'Networks', 'FactorGraphs');
end


end




function [conn, Adj, diameter, colors] = read_Adjacency_and_connectivity(filename)
% [conn, Adj, diameter, colors] = read_Adjacency_and connectivity(filename)
%
% It reads values from a file with 'filename' with the following format:
% 
% conn = val
% C = val
% P = val
% diameter = val
% Adj =
% val val val ... val
% val val val ... val
% ...
% val val val ... val
% colors = 
% val val val ... val
%
% where val is a number. The size of Adj matrix is (C+P) x (C+P). 



fid = fopen(filename, 'r');
data = textscan(fid, '%s');
fclose(fid);

conn = str2double(data{1}(3));

if conn == 0
    Adj = [];
    colors = [];
    diameter = [];
    return;
end

C = str2double(data{1}(6));
P = str2double(data{1}(9));
diameter = str2double(data{1}(12));

num_nodes = C + P;

Adj = zeros(num_nodes,num_nodes);

index = 14;

for i = 1 : num_nodes
    for j = 1 : num_nodes
        Adj(i,j) = str2double(data{1}(index + num_nodes*(i-1) + j));
    end
end

index = 14 + num_nodes^2 + 2;

colors = zeros(num_nodes,1);

for i = 1 : num_nodes
    colors(i) = str2double(data{1}(index + i));
end

end
