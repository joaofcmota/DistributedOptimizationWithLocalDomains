function [g, vars_BB] = Delays_grad(x, vars_BB)

v = vars_BB.v;
a = vars_BB.a;
capacities = vars_BB.capacities;

% If there are entries =c_i, the derivative is infinity. To avoid it,
ind_plus_c = (x ==  capacities);

x(ind_plus_c) =  capacities(ind_plus_c) - 1e-3;

g = capacities./(capacities - x).^2 + v + 2*a.*x;

