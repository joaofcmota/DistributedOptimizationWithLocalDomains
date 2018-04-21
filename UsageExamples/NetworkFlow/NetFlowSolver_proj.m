function [proj] = NetFlowSolver_proj(y, a, b, c)

% [proj] = NetFlowSolver_proj(y, a, b, c)
%
% Projects the point y onto the set {a : a'*x = b, 0 <= x <= c}. This boils
% down to solving
%
%                   minimize    0.5*||x - y||^2
%                       x
%                   subject to  a'*x = b
%                               0 <= x <= c
%
% The solution to this problem can be computed in a very simple way,
% without requiring any iterative solver.

% Check for input errors
if size(a) ~= size(y)
    error('Vectors a and y in ''NetFlowSolver_proj'' do not have the same dimensions.');
end
n = length(a);

if sum(abs(a)) ~= n
    error('Vector a in ''NetFlowSolver_proj'' should contain only +1 or -1 in its entries.');
end

ind_a_pos = (a >= 0);
ind_a_neg = (a <  0);

points_of_slop_change = [y(ind_a_pos) ; y(ind_a_pos) - c(ind_a_pos) ;...
    -y(ind_a_neg); c(ind_a_neg) - y(ind_a_neg)];

z = sort(points_of_slop_change);

% Evaluate the function g(lambda) = a'*max( y - lambda*a , 0 ) for all 
% lambda = z(i), i = 1, ..., len_z

len_z = length(z);
g_z = zeros(len_z,1);

for i = 1 : len_z    
    g_z(i) = Evaluate_g(y, a, z(i), c);
end

% Note: g_z should contain decreasing values

if b > g_z(1) || b < g_z(len_z)
    error('The projection problem passed to the function ''NetFlowSolver_proj'' is infeasible');
end

indices = (g_z <= b);
aux = find(indices);
l = aux(1)-1;

if l == 0
    lambda = z(1);
else
    lambda = z(l) + (z(l+1) - z(l))*(b - g_z(l))/(g_z(l+1) - g_z(l));
end


% Compute proj

proj = y - lambda*a;

for i = 1 : n        
    if proj(i) < 0
        proj(i) = 0;
    elseif proj(i) > c(i)
        proj(i) = c(i);
    end
end


end


function [g] = Evaluate_g(y, a, lambda, c)

n = length(y);

g = 0;

for i = 1 : n
   
    point = y(i) - lambda*a(i);
    
    if point < 0
        point = 0;
    elseif point > c(i)
        point = c(i);
    end
    
    g = g + a(i)*point;
end


end

