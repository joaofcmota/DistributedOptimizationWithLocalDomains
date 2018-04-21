% Reads the power grid network and stores the .mat file in
%
%  ../ProcessedNetwork/
%
% It uses the .sage files to create a tmp file with the network


% =========================================================================
% Filenames

% Sage script that reads the network, colors it, and gets other stats
SAGE_FILE = 'Barabasi.sage';
%SAGE_FILE = 'processPowerGrid.sage';

% Temporary file created by the Sage script; it will be erased aftewards
TMP_FILE = 'tmp.txt';     

% Where the final network will be stored (_node_num.mat will be added)
FILENAME_OUTPUT = '../ProcessedNetwork/Barabasi';
%FILENAME_OUTPUT = '../ProcessedNetwork/Network_PowerGrid.mat';
% =========================================================================


% =========================================================================
% Execute Sage script

% Name of temporary file for storing the statistics. It will be moved after
% knowing what the number of nodes is
STATS_TMP_FILE = 'stats_tmp.txt'; 

system(['sage ', SAGE_FILE, ' > ', STATS_TMP_FILE]);
% =========================================================================



% =========================================================================
% Read TMP_FILE
%
% The file has the following format:
%
% conn = num
% num_nodes = num
% num_edges = num
% diameter = num
% numColors = num
% nodePerColor =
% num
% num
%  ...
% num
% colors =
% num num num num num num num num num num num num num num
% num num num num
% num num num num num num num num
% Edges =
% num num
% num num
% num num
%  ...

fid = fopen(TMP_FILE, 'r');
data_aux = textscan(fid, '%s');
fclose(fid);

data = data_aux{1};

% Delete TMP file
system(['rm ', TMP_FILE]);

conn = str2double(data(3));

if conn == 0
    error('Network is not connected.')
end

ind = 6;
num_nodes = str2double(data(ind));
ind = ind + 3;
num_edges = str2double(data(ind));
ind = ind + 3;
diameter  = str2double(data(ind));
ind = ind + 3;
numColors = str2double(data(ind));

nodesPerColor = zeros(numColors,1);

% Read nodesPerColor
ind = ind + 3;
for i = 1 : numColors
    nodesPerColor(i) = str2double(data(ind));
    ind = ind +1;
end

% Read Colors
Colors = cell(numColors,1);

ind = ind + 2;
for i = 1 : numColors
    vec_aux = zeros(1,nodesPerColor(i));
    for nd = 1 : nodesPerColor(i)
        vec_aux(nd) = str2double(data(ind)) + 1;  % In Python vectors start
                                                  % in 0
        ind = ind + 1;
    end
    Colors{i} = vec_aux;
end

% Read Edges and build adjacency matrix
Adj = sparse(zeros(num_nodes,num_nodes));

ind = ind + 2;
for i = 1 : num_edges
    node1 = str2double(data(ind)) + 1;     % In python nodes start in 0
    ind = ind + 1;
    node2 = str2double(data(ind)) + 1;
    ind = ind + 1;
    Adj(node1,node2) = 1;
    Adj(node2,node1) = 1;
end
% =========================================================================


% =========================================================================
% Store Network

% First, get Neighbors, Degrees, and Max_Degree
Neighbors = cell(num_nodes,1);
Degrees = zeros(num_nodes,1);

for p = 1 : num_nodes
    Neighbors{p} = find(Adj(p,:));
    Degrees(p) = length(Neighbors{p});
end

Max_Degree = max(Degrees);

Network = struct('P', {num_nodes}, ...
    'Adj', {Adj}, ...
    'Type', 'Power Grid US', ...
    'Num_Colors', {numColors}, ...
    'Partition', {Colors}, ...
    'Diameter', {diameter}, ...
    'Neighbors', {Neighbors}, ...
    'Degrees', {Degrees}, ...
    'Max_Degree', {Max_Degree} ...
    );


% FILE_NAMES
ending_filename = [FILENAME_OUTPUT, '_P_', num2str(num_nodes)];

filename_data = [ending_filename, '.mat'];
filename_stat = [ending_filename, '.txt'];

save(filename_data, 'Network');
system(['mv ' STATS_TMP_FILE, ' ', filename_stat]);
% =========================================================================




