function [f, vars_BB] = Delays_f(x, vars_BB)

v = vars_BB.v;
a = vars_BB.a;
capacities = vars_BB.capacities;


f = sum(x./(capacities - x)) + v'*x + a'*(x.^2);



