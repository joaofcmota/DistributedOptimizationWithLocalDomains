function [error_p, vars_prob] = Error_gradNetFlow(X, x_opt, vars_prob)

x_estimate = vars_prob.x_estimate;

error_p = norm(x_estimate - x_opt, Inf)/norm(x_opt, Inf);

