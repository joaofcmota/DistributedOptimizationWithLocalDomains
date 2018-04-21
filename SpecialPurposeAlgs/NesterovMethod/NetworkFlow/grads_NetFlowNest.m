function [grad, vars_prob] = grads_NetFlowNest(p, x, vars_prob)

% ASSUMPTION: Between any pair of nodes i and j, there is only one link at
% most, either from i to j, or from j to i


% From vars_prob
capacities = vars_prob.capacities;
B = vars_prob.B;
d = vars_prob.d;
components = vars_prob.components;
x_estimate = vars_prob.x_estimate;

num_edges = length(capacities);

comp_p = components{p};


% To avoid confusion with notation
lambda = x;

ind_lambda_p = find(comp_p == p);

lambda_p = lambda(ind_lambda_p);

n_p = length(lambda);


grad = zeros(n_p,1);

for i_n_p = 1 : n_p
   
    if i_n_p == ind_lambda_p       % contribution to the gradient of itself
        grad(i_n_p) = d(p);
    else        
        node = comp_p(i_n_p);      % node for which we want to compute the 
                                   % contribution of node p to its gradient
        lambda_j = lambda(i_n_p);
    
        % Find the index in B (and lambda) that corresponds to edge (p,node)
        for ed = 1 : num_edges
            if abs(B(p,ed)) == 1 && abs(B(node,ed)) == 1
                break;
            end
        end
    
        if ed == num_edges && abs(B(p,ed)) ~= 1 && abs(B(node,ed)) ~= 1
            error('Programming error')
        end
    
        if B(p,ed) == -1    % Link out node p
            grad(i_n_p) =  -compute_x_our(lambda_j - lambda_p, capacities(ed));
            x_estimate(ed) = -grad(i_n_p);
        else
            grad(i_n_p) = compute_x_our(lambda_p - lambda_j, capacities(ed));
        end
    end
end

vars_prob.x_estimate = x_estimate;

end


function [x_our] = compute_x_our(v, capacity)

if v >= 0
    x_our = 0;
    return;
end

x_our_quad = capacity - sqrt(-capacity/v);

if x_our_quad < 0
    x_our = 0;
elseif x_our_quad > capacity
    x_our = capacity;   
else
    x_our = x_our_quad;
end
end
