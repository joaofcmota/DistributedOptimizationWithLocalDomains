function [xp, vars_prob] = SolverMPC(p, v, c, X, vars_prob)

% From vars_prob
E_p = vars_prob.E_p{p};
w_p = vars_prob.w_p{p};

M = E_p + diag(c);
vec = w_p + v;

xp = -0.5*( M \ vec);

