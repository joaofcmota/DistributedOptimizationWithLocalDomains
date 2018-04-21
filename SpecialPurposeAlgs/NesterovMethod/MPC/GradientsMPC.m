function [grad, vars_prob] = GradientsMPC(p, x, vars_prob)

% From vars_prob
E_p = vars_prob.E_p{p};
w_p = vars_prob.w_p{p};

grad = 2*E_p*x + w_p;

