function [xp, vars_prob] = NetFlowSolver(p, v, a, X, vars_prob)

% Solver for Network Flow problems using DADMM_Partial. We use the
% Barzilai-Borwein method (see function SPG_BB).

% Problem variables
B = vars_prob.B;
d = vars_prob.d;
Dp = vars_prob.Dp;
capacities = vars_prob.capacities;

Delays_f           = vars_prob.function_eval;
Delays_grad        = vars_prob.gradient_eval;
Delays_proj        = vars_prob.projection;
NetFlowSolver_proj =  vars_prob.NetFlowSolver_proj;
SPG_BB             = vars_prob.spg_bb;


dp = d(p);

indices_of_x = (B(p,:) ~= 0);

dim_var = sum(indices_of_x);

if Dp(p) ~= dim_var
    error('Dp(p) ~= sum(indices_of_x) in function NetFlowSolver');
end


% Solve the problem with SPG_BB

const = B(p, indices_of_x)';

vars_BB = struct('f_val', {Delays_f}, ...
    'grad', {Delays_grad}, ...
    'proj', {Delays_proj}, ...
    'NetFlowSolver_proj', {NetFlowSolver_proj}, ...
    'v', {v}, ...
    'a', {a}, ...
    'capacities', {capacities(indices_of_x)}, ...
    'const', {const}, ...
    'd', {dp} ...
    );


[xp, vars_BB] = SPG_BB(dim_var, X{p}, vars_BB);





