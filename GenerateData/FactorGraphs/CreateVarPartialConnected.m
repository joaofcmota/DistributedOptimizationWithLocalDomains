function [FactorGraph] = CreateVarPartialConnected(n, P, neighbors, depth)

% [FactorGraph] = CreateVarPartialConnected(n, P, neighbors, depth)
%
% Creates a random instance of a problem with a partial connected variable.
% Given a network with P nodes, it implements the following algorithm. For
% each component xl, it chooses one node at random; then, it selects one of
% the neighbors of this node, again, at random. Next, it selects one node
% from the set of neighbors of the first and the second nodes. This is
% performed a number 'depth' of times. The function at the nodes in the 
% resulting subgraph will depend on xl. This process is repeated for all
% components.
%
% Different mode (n == P)
%   When n == P, the algorithm is slightly different. Namely, there are P
%   variables and P nodes. So, the pth variable is assigned to node p. The
%   nodes that will depend on xp will be determined by the neighbors and 
%   depth in the same way described above.
%
% Inputs:
%   - n: dimension of the variable of the problem
%   - P: number of nodes in the network
%   - neighbors: Px1 cell. The pth entry has a vector with the neighbors of
%                node p.
%   - depth: a number between 2 and P-1. It is the number of nodes in each 
%            subgraph. It can also be a vector of size n, where each entry
%            specifies the number of nodes in each subgraph (each subgraph 
%            is generated n times).
%
% Output: Struct with the following fields:
%   - components: Px1 cell, where the pth entry contains a vector with the
%                 indices of the components node p depends on.
%   - P: number of nodes in the network
%   - neighbors: Px1 cell. The pth entry has a vector with the neighbors of
%                node p.
%   - depth: same as in the input 
%   - n: same as in the input
%  
% =========================================================================
% NOTE: If at the end of the algorithm there are nodes which were not
%       assigned any variable, the algorithm is repeated. The number of
%       times the algorithm is repeated is defined in the variable REPEAT,
%       in the code.
% =========================================================================

REPEAT = 30;    % Maximum number of times the algorithm can be repeated


% =========================================================================
% Check for errors in the input
if isscalar(depth)
    depth = depth*ones(n,1);
else
    if length(depth) ~= n
        error('Input vector depth should have n entries or be a scalar. Type ''help CreateVarPartialConnected''.');
    end
end
if sum( (depth > (P-1)) + (depth < 2) ) > 0
    error('Input ''depth'' should be between 2 and P-1. Type ''help CreateVarPartialConnected''.')
end
% =========================================================================

% Determine the mode
if n == P
    NORMAL_MODE = 0;
else
    NORMAL_MODE = 1;
end


components = cell(P,1);


% =========================================================================
% Algorithm

for rep = 1 : REPEAT
    
    for p = 1 : P
        components{p} = [];
    end
    
    for l = 1 : n

        fringe = [];
        nodes_in_subgraph = zeros(P,1);   % pth entry = 1 => node p is in the graph
                                          %           = 0 => otherwise        
        
        if NORMAL_MODE                                  
            p = round((P-1)*rand) + 1;    % node selected randomly
        else
            p = l;
        end
        
        components{p} = [components{p} ; l];
        
        fringe = alterfringe(fringe, neighbors{p}, 'add');
        
        nodes_in_subgraph(p) = 1;
        
        for dep = 1 : depth(l)-1
            
            len_fringe = length(fringe);
            k = round((len_fringe-1)*rand) + 1;  % node selected randomly 
                                                 % from the fringe
            p = fringe(k);
            
            components{p} = [components{p} ; l];
            
            nodes_in_subgraph(p) = 1;
            
            fringe = alterfringe(fringe, p, 'rem');
            fringe = alterfringe(fringe, neighbors{p}, 'add');
            fringe = alterfringe(fringe, find(nodes_in_subgraph), 'rem');
            
            if sum(nodes_in_subgraph) == P-1
                fprintf('Warning: component %d contains P-1 nodes. Premature stop.\n', l);
                break;
            end
        end
    end
    
    
    % Check if there are nodes without any component
    node_without_comp = zeros(P,1);
    for p = 1 : P
        if isempty(components{p})
            node_without_comp(p) = 1;
        end
    end
    num_nodes = sum(node_without_comp);
    if num_nodes > 0
        %fprintf('There are %d nodes without components. Repeating (repetition number = %d)\n', num_nodes, rep);
    else
        break;
    end
        
end
% =========================================================================


% =========================================================================
% Output

if rep == REPEAT
    error('The generated problem had always at least one node with no assigned variables. Please increase depth or n next time.');
else
    fprintf('Number of repetitions: %d \n', rep);
end

FactorGraph = struct('components', {components}, ...
    'P', {P}, ...
    'neighbors', {neighbors}, ...
    'depth', {depth}, ...
    'n', {n} ...
    );
% =========================================================================


end




function [fringe] = alterfringe(fringe, node, operation)

% Receives the vector 'fringe' and adds or removes 'node', which is either 
% a vector of nodes or a scalar, according to the value of 'operation'. The 
% input operation is a string with 'add' if we want to add 'node' to the 
% 'fringe', or with 'rem'if we want to remove 'node' from the fringe. 
    
switch operation
    case 'add'
        
        for i = 1 : length(node)
            
            if sum(fringe == node(i)) == 0  % Only add the node if it is not in the fringe
                fringe = [fringe ; node(i)]; %#ok<AGROW>
            end
            
        end
        
        
    case 'rem'
        
        for i = 1 : length(node)
        
            index = (fringe == node(i));
            fringe = fringe(not(index));
            
        end
        
        
    otherwise
        error('Operation not recognized in alterfringe, file CreateVarPartialConnected.');
end
    

end
