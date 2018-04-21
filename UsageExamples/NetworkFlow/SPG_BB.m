function [x_opt, vars_prob] = SPG_BB(n, x0, vars_prob)

% [x_opt, vars_prob] = SPG_BB(n, x0, vars_prob)
%
% Implements a Spectral Projected Gradient (SPG) Method, using a technique
% similar to that of Barzilai-Borwein. More specifically, it implements the
% algorithm SPG1 from
%
% E. Birgin, J. Martinez, and M. Raydan, "Nonmonotone Spectral Projected
% Gradient Methods on Convex Sets," SIAM J. Optim, Vol.10, No.4, 2000, pp.
% 1196-1211
%
% SPG_BB solves
%
%                         minimize   f(x)
%                         subject to x in X
%
% where f(x) is a convex, continuously differentiable function, and X a
% closed convex set. It is required from the user to provide three 
% functions in the struct vars_prob. In the field 'f_val' should be a 
% handler for
%
%        [f, vars_prob] = function_name(x, vars_prob)
%
% where x is a vector of size n, and f is the value of the function f at x,
% i.e., f(x). In the field 'grad' should be a handler for
%
%        [g, vars_prob] = function_name(x, vars_prob)
%
% where x is a vector of size n, and g is the gradient of f at the point x.
% In the field 'proj' should be a handler for
%
%        [p, vars_prob] = function_name(x, vars_prob)
%
% where x is has size n, and p is the projection of x onto the set X, i.e.,
% p = argmin{ ||z-x|| : z in X}.
%
% The remaining fields of the struct vars_prob can be used by the user to 
% store information useful for the implementation of the above functions.
%
% The SPG_BB inputs and outputs are:
%
% Inputs:
%
%   - n: size of the variable
%
%   - x0: initialization vector n x 1. If empty ([]), the algorithm is
%         initialized with the vector of zeros
%
%   - vars_prob: struct with the following compulsory fields (explained 
%                above):
%                  . f_val
%                  . grad
%                  . proj
%                The remaining struct can be used to optimize the
%                algorithm.
%
% Outputs:
%
%   - x_opt: optimal value for the variable x
%
%   - vars_prob: same struct as input
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Please send any questions, comments, or bug reports to joaomota@cmu.edu.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% =========================================================================
% Check for input errors

if ~isstruct(vars_prob)
    error('vars_prob (3rd argument of SPG_BB) should be a struct. Please type ''help SPG_BB''.');
end
if ~isfield(vars_prob, 'f_val')
    error('vars_prob (3rd argument of SPG_BB) should contain the field ''f_val''. Please type ''help SPG_BB''.');
end
if ~isfield(vars_prob, 'grad')
    error('vars_prob (3rd argument of SPG_BB) should contain the field ''grad''. Please type ''help SPG_BB''.');
end
if ~isfield(vars_prob, 'proj')
    error('vars_prob (3rd argument of SPG_BB) should contain the field ''proj''. Please type ''help SPG_BB''.');
end

if isempty(x0)
    x_init = zeros(n,1);
else
    x_init = x0;
    if length(x_init) ~= n
        error('Input vector x0 is not of size n. Please type ''help SPG_BB''.');
    end
end
% =========================================================================


% =========================================================================
% Parameters of the algorithm (we will use the ones suggested by the paper)

MAX_ITER = 1e3;     % Maximum number of iterations (warning if reached)

MAX_ITER_IN = 1e3;  % Maximum number of iterations for the inner loop

epsilon = 1e-10;     % Algorithm stops if ||P(x_k - g(x_k))-x_k||<= epsilon

M = 30;             % Memory of the nonmonotone algorithm

gamma = 1e-4;       % Sufficient decrease parameter

alpha_min = 1e-6;   % Maximum and minimum value for the stepsize
alpha_max = 1e6;

sigma_1 = 0.1;      % Safeguarding parameters
sigma_2 = 0.9;     
% =========================================================================


% =========================================================================
% Algorithm

% Function handlers
f = vars_prob.f_val;
g = vars_prob.grad;
proj = vars_prob.proj;

% Vector with the value of f for the past M iterations
past_f = Inf*ones(M,1);

% Initializations
x_k = x_init;
[g_k, vars_prob] = g(x_k, vars_prob);
[proj_xg, vars_prob] = proj(x_k - g_k, vars_prob);
alpha = min(alpha_min, 1/max(abs(proj_xg - x_k)));         % recommended by the paper's authors 

for iter = 1 : MAX_ITER
   
    lambda = alpha;
    
    for inner_iter = 1 : MAX_ITER_IN
        [x_new, vars_prob] = proj(x_k - lambda*g_k, vars_prob);
        
        gamma_inner_prod = gamma*(x_new - x_k)'*g_k;
        gamma_inner_prod_vec = ones(M,1)*gamma_inner_prod;
        
        max_val = max(past_f + gamma_inner_prod_vec);
        [f_new, vars_prob] = f(x_new, vars_prob);
        
        if f_new <= max_val
            x_prev = x_k;
            g_prev = g_k;
            x_k = x_new;
            s_k = x_k - x_prev;
            [g_k, vars_prob] = g(x_k, vars_prob);                        
            y_k = g_k - g_prev;
            break;
        else
            lambda_new = lambda/2;
            lambda_new = max(lambda_new, sigma_1*lambda);
            lambda = min(lambda_new, sigma_2*lambda);
        end
    end
    
    if inner_iter == MAX_ITER_IN
        fprintf('Warning: Maximum number of iterations reached in the inner loop of SPG_BB\n');
    end
    
    % Store the value of f in past_f
    [f_k, vars_prob] = f(x_k, vars_prob);
    index = mod(iter, M) + 1;
    past_f(index) = f_k;
    
    b_k = s_k'*y_k;
    
    if b_k <= 0
        alpha = alpha_max;
    else
        a_k = norm(s_k)^2;
        alpha = min( alpha_max , max(alpha_min, a_k/b_k) );
    end
            
    % Compute f and g at point x_k
    [g_k, vars_prob] = g(x_k, vars_prob);
    [proj_xg, vars_prob] = proj(x_k - g_k, vars_prob);
    
    % Stopping criterion
    if max(abs(proj_xg - x_k)) <= epsilon
        break;
    end    
end
% =========================================================================

%iter
if iter == MAX_ITER
    fprintf('Warning: Maximum number of iterations reached in SPG_BB\n');
end

x_opt = x_k;










