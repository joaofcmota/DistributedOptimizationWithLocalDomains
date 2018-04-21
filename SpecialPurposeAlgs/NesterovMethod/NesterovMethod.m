function [X, vars_prob, varargout] = NesterovMethod(n, Lipschitz, ...
    vars_prob, vars_network, varargin) 

% [X, vars_prob, varargout] = NesterovMethod(n, Lipschitz, ...
%   vars_prob, vars_network, varargin)
%
% Solves the optimization problem
%
%       minimize  f1(x_S1) + f2(x_S2) + ... + fP(x_SP)         (1)
%          x
% 
% where the variable x has n components: x = (x1, x2, ..., xn). It assumes
% that the subgraph induced by each component is a star. This is an
% implementation of Nesterov's algorithm directly applied to (1). Each
% function fp depends only on a subset Sp of components of x.
%
% Requirements on each fp:
%     - fp is convex, real-valued, differentiable
%     - the gradient of f1 + ... + fP is Lipschitz continuous, where the
%       Lipschitz constant is passed by the user as input
%     - constraints can be imposed by the nodes that are star centers;
%       these constraints have to be imposed componentwise, not on blocks
%       of components
%
% The user has to provide two functions:
%  1)  one that returns the gradients for each function fp
%  2)  a function that projects a point onto the set of constraints, for
%      each component
%
% The user gives the handler of these functions the input vars_prob. This 
% input is a struct that can be used to store information/variables 
% relevant to solve (1), but it should contain the fields 'components', 
% 'centers', 'gradients', and 'proj_constraints', as explained next.
%
% Complusory fields of the input struct vars_prob:
%
% components:
%   a cell of size Px1, where P is the number of nodes in the network. The 
%   pth entry contains a vector with the components node p depends on. For
%   example, if Sp = {3, 20, 22}, then the pth entry of components should 
%   have components{p} = [3, 20, 22]. We assume each vector to be ordered 
%   and with no repetitions.
%
% centers:
%   a vector of size nx1, where the lth entry contains the node that is the
%   star for component x_l.
%
% gradients:
%   a handler of a function (written by the user) that returns the gradient
%   of fp. The header of the function should be
%    
%       [grad, vars_prob] = function_name(p, x, vars_prob)
%
%   where p in the node number and x the point at which we want the 
%   gradient to be computed. The ordering of the components of x is the 
%   same as in 'components'. The output grad should have the same order.
%
% proj_constraints:
%   a handler of a function (written by the user) that returns the
%   projection of a point x_l onto a set X_l of constraints. The header is
%
%       [proj_point, vars_prob] = function_name(l, point, vars_prob) 
%
%   where l is the component number (between 1 and n), point is a scalar,
%   and proj_point is the projection of point onto the set of constraints.
%
% 
% The input vars_prob should then contain vars_prob.components, 
% vars_prob.centers, vars.gradients, and vars.proj_constraints. The 
% remaining fields can be used to store problem variables.
%
% The remaining inputs and outputs of NesterovMethod are now explained.
%
% Inputs:
%     - n: is the size of the variable x.
%
%     - Lipschitz: is a positive number: the Lipschitz constant of the
%                  global cost function of (1).
% 
%     - vars_prob: is a struct that contains at least four fields, 
%                  'components', 'centers', 'gradients', and 
%                  'proj_constraints' as explained above. The rest of 
%                  vars_prob can be used by the user to store information 
%                  (so that it does not get erased when Matlab leaves that 
%                  function).
%
%     - vars_network: is a struct containing information about the
%                     network: 
%           . vars_network.P is the number of nodes P.                      
%           . vars_network.neighbors is a cell array of the size of the 
%             number of nodes P, and each entry i contains the neighbors of
%             node i.           
%
%     - varargin: the optional argument is a struct which contains the
%                 following fields.
%           . max_iter: maximum number of iterations. The default is 1000.
%           . eps: the algorithm stops either because of max_iter or 
%                  because for all nodes ||X{p} - X_prev{p}|| < eps. It 
%                  defaults to 1e-4.
%           . init: all components of the variable are initialized with
%                   init. It defaults to 0.
%           . x_opt: the optimal solution of (1). If the optimal solution 
%                    is passed, several (optional) outputs are activated.
%           . error_fun: only valid if x_opt exists. It is a function
%                        handler to assess the error, e.g., if the user 
%                        wants to see the error in the dual variable. The
%                        header should be:
%
%             [error, vars_prob] = name_of_function(X, x_opt, vars_prob)
%
%                        where X is a cell of size P where X{p} is the
%                        estimate of node p to problem (1), x_opt is the
%                        optimal solution provided by the user, and
%                        vars_prob is the struct also provided by the user.
%                        This function is called at the end of each
%                        iteration (all nodes have updated their estimates)
%
%           . eps_opt: only valid in case x_opt exists. This turns off
%                      the eps stopping criteria (based on two consecutive
%                      iterations), and replaces it by
%                         || x_est - x_opt ||_Inf / ||x_opt||_Inf,
%                      where each component of x_est is the worst estimate
%                      over all nodes that estimate that component. eps_opt
%                      has to be > 0.
%           . turn_off_eps: if this field is present (no matter what
%                           value), the eps or the eps_opt stopping 
%                           criteria are turned off: NesterovMethod terminates 
%                           only when it reaches the maximum number of 
%                           iterations.
%
%
% Outputs:
%
%     - X: is a P x 1 cell array where the pth entry contains the solution
%          x, solving the optimization problem, for node p. X{p} has size
%          n_p = |S_p| and the components are by increasing order.
%
%     - vars_prob: is the struct that the user provides as input. Since it
%                  may have relevant information, we return it to the user.
%
%     - varargout: (optional) it is a struct with the following fields.
%           . iterations: total number of iterations
%           . stop_crit: string with 'MAX ITERATIONS REACHED', if the 
%                        maximum number of iterations was reached, or
%                        'EPS REACHED', otherwise.
%           . error_iterations: this output is activated when x_opt is 
%                               passed as input. It is a vector with the
%                               size of the number of iterations, and
%                               contains the relative error in each entry 
%                               (corresponding to each iteration)
%           . iter_for_errors: activated when x_opt is passed. It is a
%                              matrix 2x10:
%                                          [  1e-1  ,  iter1 ]
%                                          [  1e-2  ,  iter2 ]
%                                          [       ...       ]
%                                          [  1e-10 , iter10 ]
%                              that contains the number of iterations to
%                              achieve a "canonical" relative error of
%                              1e-1, 1e-2, ..., 1e-10.
%
%
% *************************************************************************
% NOTE: error_iterations and iter_for_errors contain the errors relative to
%       the variable x. If the user provides its own error function (in 
%       error_fun), then these outputs will be in respect to that error.
% *************************************************************************
%
% -------------------------------------------------------------------------
% Distributed Optimization With Local Domains: Applications in MPC and
% Network Flows
% Copyright (C) 2014 Joao Mota
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
%
% =========================================================================
% Please send any questions, comments, or bug reports to j.mota@ucl.ac.uk
% =========================================================================


OPT_OUTPUT = 0;    % 0 if there is no optional output; 1 otherwise

% =========================================================================
% Check for input and output errors
optargs = nargin - 4;
if optargs > 1
    error('Number of input arguments exceeds the maximum value. Please type ''help NesterovMethod''.');
end
if optargs == 1
    if ~isstruct(varargin{1})
        error('Optional input should be a struct. Please type ''help NesterovMethod''.');
    end
end
if Lipschitz <= 0
    error('Lipschitz constant (2nd argument of NesterovMethod) is not positive. Please type ''help NesterovMethod''.');
end
if ~isstruct(vars_prob)
    error('vars_prob (3rd argument of NesterovMethod) should be a struct. Please type ''help NesterovMethod''.');
end
if sum(isfield(vars_prob, {'components', 'centers', 'gradients', 'proj_constraints'})) ~= 4
    error('vars_prob (3rd argument of NesterovMethod) should contain the fields ''components'', ''centers'', ''gradients'', and ''proj_constraints''. Please type ''help NesterovMethod''.');
end
if ~isstruct(vars_network)
    error('vars_network 4th argument of NesterovMethod) should be a struct. Please type ''help NesterovMethod''.');
end
if sum(isfield(vars_network, {'P', 'neighbors'})) ~= 2
    error('Fieldnames of the struct vars_network are not correct. Please type ''help NesterovMethod''.');
end
nout = max(nargout,1)-2;
if nout > 1
    error('Number of outputs is greater than the maximum. Please type ''help NesterovMethod''.'); 
end
if nout == 1
    OPT_OUTPUT = 1; 
end
% =========================================================================


% =========================================================================
% Take care of optional input

% Predefined and default variables
MAX_ITER = 1000;
EPS = 1e-4;
INIT_X = 0;

% Optional input
EXISTS_X_OPT = 0;     % 0 if there is no x_opt; 1 otherwise
TURN_OFF_EPS = 0;     % 1 if stopping criterion is maximum number of 
                      % iterations only; 0 otherwise
                      
EPS_OPT = 0;          % Equal to 0 if turned off. Oth., it has the eps_opt value 

if ~isempty(varargin)
    opt_input = varargin{1};
    
    if isfield(opt_input, 'max_iter')
        MAX_ITER = opt_input.max_iter;
    end
    if isfield(opt_input, 'eps')
        EPS = opt_input.eps;
    end
    if isfield(opt_input, 'init')
        INIT_X = opt_input.init;
    end
    if isfield(opt_input, 'x_opt')
        x_opt = opt_input.x_opt;
        EXISTS_X_OPT = 1;
    end
    
    ERROR_FUN = 0;        % 1 if user provides error function
    
    if EXISTS_X_OPT == 1
        if isfield(opt_input, 'error_fun')
            ERROR_FUN = 1;
        end
        if isfield(opt_input, 'eps_opt')
            EPS_OPT = opt_input.eps_opt;
        end
    end
    if isfield(opt_input, 'turn_off_eps')
        TURN_OFF_EPS = 1;
    end
end

% Network variables
P = vars_network.P;
neighbors = vars_network.neighbors;

% Problem variables
components = vars_prob.components;
centers = vars_prob.centers;
gradients = vars_prob.gradients;
proj_constraints = vars_prob.proj_constraints;


if length(components) ~= P
    error('Size of components is not P. Please type ''help NesterovMethod''.');
end

if length(centers) ~= n
    error('Size of centers is not n. PLease type ''helop NesterovMethod''.')
end
% =========================================================================


% =========================================================================
% Initializations

Sp = zeros(P,1);                % Number of components each node depends on
                                % node p depends on at the pth entry
Grads = cell(P,1);              % Cell (one/node) with the gradients
X = cell(P,1);                  % Cell (one/node) with current estimates
Y = cell(P,1);                  % Nesterov variable
X_prev = cell(P,1);             % Cell (one/node) with previous estimates

for p = 1 : P
    Sp(p) = length(components{p});
    X{p} = INIT_X*ones(Sp(p),1);
    Y{p} = X{p};
    X_prev{p} = X{p};
    Grads{p} = zeros(Sp(p),1);
end


% n x 1 cell where the lth entry contains V_l, i.e., the set of nodes that
% depends on x_l.
induced_subgraphs = cell(n,1);

for l = 1 : n
    
    induced_subgraph_x_l = [];    
    for p = 1 : P
        if sum(components{p} == l) > 0
            induced_subgraph_x_l = [induced_subgraph_x_l, p]; %#ok<AGROW>
        end
    end
    induced_subgraphs{l} = induced_subgraph_x_l;
    
    % Check if each induced subgraph is a star    
    center_of_star = centers(l);
    neighbs_cs = [neighbors{center_of_star} , center_of_star];
    
    len_V_l = length(induced_subgraph_x_l);
    for i = 1 : len_V_l
        if sum(induced_subgraph_x_l(i) == neighbs_cs) == 0
            error('The subgraph induced by x_%d is not a star.', l);
        end
    end
end


Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated

                                      
% Struct with the internal variables of NesterovMethod. 
vars = struct( 'Sp', {Sp}, ...
    'induced_subgraphs', {induced_subgraphs}, ...
    'Grads', {Grads}, ...
    'X', {X},...
    'Y', {Y}, ...
    'X_prev', {X_prev},...    
    'EPS', {EPS},...
    'MAX_ITER', {MAX_ITER},...
    'Stop', {Stop},...
    'TURN_OFF_EPS', {TURN_OFF_EPS},...    
    'EPS_OPT', {EPS_OPT},...
    'EXISTS_X_OPT', {EXISTS_X_OPT},...
    'x_opt', {0},...                    % To be filled if EXISTS_X_OPT == 1
    'error_fun', {0},...                % To be filled if EXISTS_X_OPT == 1
    'error_fun_handler', {0},...        % To be filled if EXISTS_X_OPT == 1
    'error_iterations', {[]},...        % To be filled if EXISTS_X_OPT == 1
    'iter_for_errors', {0}, ...         % To be filled if EXISTS_X_OPT == 1
    'iterations', {0}...                % To be filled if EXISTS_X_OPT == 1
    );                                      


if EXISTS_X_OPT == 1
    if ERROR_FUN == 1
        error_fun_handler = opt_input.error_fun;
        vars.error_fun = 1;
        vars.error_fun_handler = error_fun_handler;
        [error_iterations, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
    else
        error_vec = zeros(P,1);
        for p = 1 : P
            error_vec(p) = norm(X{p} - x_opt(components{p}), Inf);
        end
        error_iterations = norm(error_vec, Inf)/norm(x_opt, Inf);
%         error_iterations = 0;
%         for p = 1 : P
%             error_iterations = error_iterations + norm(X{p} - x_opt(components{p}))^2/Sp(p);
%         end
%         error_iterations = sqrt(error_iterations/norm(x_opt)^2);
    end
    iter_for_errors = zeros(10,2);
    iter_for_errors(:,1) = 10.^(-(1:10))';    % The first column has the errors
    iter_for_errors(:,2) = Inf*ones(10,1);    % The second column the number of iterations
    iterations = 0;
    
    vars.x_opt = x_opt;
    vars.error_iterations = error_iterations;
    vars.iter_for_errors = iter_for_errors;
    vars.iterations = iterations;
end
% =========================================================================
                                        
                          
% =========================================================================
% Algorithm

alpha = 1/Lipschitz;

for k = 1 : MAX_ITER    
    
    % Compute the gradient in all nodes   
    for p = 1 : P
        
        if ~Stop(p)              % Only compute if node has not terminated
            [Grads{p}, vars_prob] = gradients(p, Y{p}, vars_prob);
        end
        
    end    
        
    % Update X and Y
    X_prev = X;
    vars.X_prev = X_prev;
    
    for l = 1 : n
        
        induced_subgraph_x_l = induced_subgraphs{l};
        len_V_l = length(induced_subgraph_x_l);
        center_x_l = centers(l);
        ind_x_l_in_center = (components{center_x_l} == l);
        
        sum_grad_y_l = 0;
        for i = 1 : len_V_l
            node = induced_subgraph_x_l(i);
            ind_x_l_in_node = (components{node} == l);
            sum_grad_y_l = sum_grad_y_l + Grads{node}(ind_x_l_in_node);
        end
        
        
        point_to_proj = Y{center_x_l}(ind_x_l_in_center) - alpha*sum_grad_y_l;
        [new_x_l, vars_prob] = proj_constraints(l, point_to_proj, vars_prob);
                
        new_y_l = new_x_l + (k-1)/(k+2)*(new_x_l - X_prev{center_x_l}(ind_x_l_in_center));
        
        % Change x_l and y_l in all nodes in the induced subgraph
        for i = 1 : len_V_l
            node = induced_subgraph_x_l(i);
            ind_x_l_in_node = (components{node} == l);
            X{node}(ind_x_l_in_node) = new_x_l;
            Y{node}(ind_x_l_in_node) = new_y_l;
        end
    end    
              
    vars.X = X;
    vars.Y = Y;
    vars.Grads = Grads;    
    
    
    if EXISTS_X_OPT == 1
        if ERROR_FUN == 1
            error_fun_handler = vars.error_fun_handler;
            [new_error, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
        else
            error_vec = zeros(P,1);
            for p = 1 : P
                error_vec(p) = norm(X{p} - x_opt(components{p}), Inf);
            end
            new_error = norm(error_vec, Inf)/norm(x_opt, Inf);
            %         new_error = 0;
            %         for p = 1 : P
            %             new_error = new_error + norm(X{p} - x_opt(components{p}))^2/Sp(p);
            %         end
            %         new_error = sqrt(new_error/norm(x_opt)^2);
        end
        error_iterations = [error_iterations, new_error];%#ok<AGROW> %figure(10);clf;semilogy(error_iterations);drawnow;
        ind_nonfilled = (iter_for_errors(:,2) == Inf);
        ind_lesserror = (new_error < iter_for_errors(:,1));
        intersect = ind_nonfilled & ind_lesserror;
        iter_for_errors(intersect,2) = k;
        
        vars.error_iterations = error_iterations;
        vars.iter_for_errors = iter_for_errors;
        vars.iterations = k;
    end
    
    
    % Stopping criterion for the outer loop
    Stop = NesterovMethod_stopping_criterion(Stop, vars_network, vars, vars_prob, k);
    
    vars.Stop = Stop;
    
    % If all nodes have converged, stop the algorithm
    if sum(Stop) == P
        break;
    end
       
    
end
% =========================================================================


% =========================================================================
% Optional output
if OPT_OUTPUT == 1
    if EXISTS_X_OPT == 1
        error_iterations_out = vars.error_iterations;
        iter_for_errors_out  = vars.iter_for_errors;
    else
        error_iterations_out = [];
        iter_for_errors_out = [];
    end
    
    if k == MAX_ITER
        stop_crit = 'MAX ITERATIONS REACHED';
    else
        stop_crit = 'EPS REACHED';
    end
    
    varargout{1} = struct('iterations', k,...
        'stop_crit', stop_crit,...
        'error_iterations', error_iterations_out,...
        'iter_for_errors', iter_for_errors_out);
end
% =========================================================================


end


function [Stop] = NesterovMethod_stopping_criterion(Stop, vars_network, vars, vars_prob, k)


if k <= 2 || vars.TURN_OFF_EPS == 1
    return
end

% Network variables
P = vars_network.P;

% Problem variables
components = vars_prob.components;

% NesterovMethod variables
% Sp = vars.Sp;
X = vars.X;
X_prev = vars.X_prev;
EPS = vars.EPS;
EPS_OPT = vars.EPS_OPT;

EXISTS_X_OPT = vars.EXISTS_X_OPT;
ERROR_FUN = vars.error_fun;

if EXISTS_X_OPT == 1
    x_opt = vars.x_opt;
    
    
    if ~isempty(EPS_OPT)
        if ERROR_FUN == 1
            error_fun_handler = vars.error_fun_handler;
            [new_error, vars_prob] = error_fun_handler(X, x_opt, vars_prob); %#ok<NASGU>
        else
            error_vec = zeros(P,1);
            for p = 1 : P
                error_vec(p) = norm(X{p} - x_opt(components{p}), Inf);
            end
            new_error = norm(error_vec, Inf)/norm(x_opt, Inf);
%             new_error = 0;
%             for p = 1 : P
%                 new_error = new_error + norm(X{p} - x_opt(components{p}))^2/Sp(p);
%             end
%             new_error = sqrt(new_error/norm(x_opt)^2);
        end
        if new_error <= EPS_OPT
            Stop = ones(P,1);
        end
        
        return;
    end
end

for p = 1 : P
    if norm(X{p} - X_prev{p}) <= EPS
        Stop(p) = 1;
    end
    
end

end