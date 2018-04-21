function [X, vars_prob, varargout] = DADMM_Partial(n, vars_prob, vars_network,...
    varargin)

% [X, vars_prob, varargout] = DADMM_Partial(n, vars_prob, vars_network,...
%  varargin)
%
% Solves the optimization problem
%
%       minimize  f1(x_S1) + f2(x_S2) + ... + fP(x_SP)         (1)
%          x
% 
% where the variable x has n components: x = (x1, x2, ..., xn). Each
% function fp in (1) is associated to the node p of a network with P nodes. 
% fp depends only on some components of x, and this is indicated by the set 
% Sp: the indices of x node p depends on. This implementation simulates a 
% distributed scenario for solving (1) with the Extended Alternating Method 
% of Multipliers (ADMM).
%
% Requirements on each fp:
%     - fp is closed, proper, and convex on R^{np}, where np = |Sp|.
%     - fp can have constraints (see below how to handle them)
%
% The user has to provide a Matlab function that solves, for each p,
%
%           minimize    fp(x_Sp) + v'*x_Sp + x_Sp'*C*x_Sp      (2)
%             x_Sp
%
% where the vector v (of size np) and the diagonal matrix C are given 
% inputs, provided by DADMM_Partial. The user gives the handler of the 
% function solving (2) through the input vars_prob. This input is a struct 
% that can be used to store information/variables relevant to solve (1), 
% but it should contain the fields 'handler' and 'components', explained 
% next.
%
% Compulsory fields of the input struct vars_prob:
%
% handler: 
%   a handler of a function (written by the user) to solve (2). The header
%   of that function should be
%
%           [xp, vars_prob] = function_name(p, v, c, X, vars_prob)
%
%   where p is the node number, v is as in (2), c is the diagonal of the 
%   matrix C, and X is a cell array such that X{p} is the last solution 
%   returned by this function, in xp, for node p. (Having X as an input 
%   reduces the memory usage whenever the solver of (2) requires 
%   warm-starts) The output xp is the solution of (2), and vars_prob is 
%   also returned because it might have changed (that depends on the user's 
%   implementation).
%
% components:
%   a cell of size Px1, where P is the number of nodes in the network. The 
%   pth entry contains a vector with the components node p depends on. For
%   example, if Sp = {3, 20, 22}, then the pth entry of components should 
%   have components{p} = [3, 20, 22]. We assume each vector to be ordered 
%   and with no repetitions.
%
% The input vars_prob should then contain vars_prob.handler and
% vars_prob.components. The remaining fields can be used to store problem
% variables.
%
% The remaining inputs and outputs of DADMM_Partial are now explained.
%
% Inputs:
%     - n: is the size of the variable x.
% 
%     - vars_prob: is a struct that contains at least two fields, 'handler'
%                  and 'components', as explained above. The rest of 
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
%           . vars_network.partition_colors is a cell array of the size of
%             the number of (proper) colors of the network graph. Each 
%             entry contains a vector with the nodes that have that color.
%             For example, for a network with 6 nodes, 
%             {[1,2], [3,4,5], 6} would mean that we have three colors: 
%             nodes 1 and 2 get the first color, 3, 4, and 5 the second, 
%             and node 6 gets the third color.
%
%
%     - varargin: the optional argument is a struct which contains the
%                 following fields.
%           . rho:  a positive constant associated to the augmented 
%                   Lagrangian of f(x). The default value is 1.
%           . max_iter: maximum number of iterations. The default is 100.
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
%                           criteria are turned off: DADMM_Partial 
%                           terminates only when it reaches the maximum 
%                           number of iterations.
%           
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
optargs = nargin - 3;
if optargs > 1
    error('Number of input arguments exceeds the maximum value. Please type ''help DADMM_Partial''.');
end
if optargs == 1
    if ~isstruct(varargin{1})
        error('Optional input should be a struct. Please type ''help DADMM_Partial''.');
    end
end
if ~isstruct(vars_prob)
    error('vars_prob (2nd argument of DADMM_Partial) should be a struct. Please type ''help DADMM_Partial''.');
end
if sum(isfield(vars_prob, {'handler', 'components'})) ~= 2
    error('vars_prob (2nd argument of DADMM_Partial) should contain the field ''handler'' and the field ''components''. Please type ''help DADMM_Partial''.');
end
if ~isstruct(vars_network)
    error('vars_network (3rd argument of DADMM_Partial) should be a struct. Please type ''help DADMM_Partial''.');
end
if sum(isfield(vars_network, {'P', 'neighbors', 'partition_colors'})) ~= 3
    error('Fieldnames of the struct vars_network are not correct. Please type ''help DADMM_Partial''.');
end
if ~iscell(vars_network.partition_colors)
    error('vars_network.partition_colors is not a cell. Please type ''help DADMM_Partial''.');
end
nout = max(nargout,1)-2;
if nout > 1
    error('Number of outputs is greater than the maximum. Please type ''help DADMM_Partial''.'); 
end
if nout == 1
    OPT_OUTPUT = 1; 
end
% =========================================================================


% =========================================================================
% Take care of optional input

% Predefined and default variables
rho = 1;
MAX_ITER = 100;
EPS = 1e-4;
INIT_X = 0;

% Optional input
EXISTS_X_OPT = 0;     % 0 if there is no x_opt; 1 otherwise
TURN_OFF_EPS = 0;     % 1 if stopping criterion is maximum number of 
                      % iterations only; 0 otherwise
                      
EPS_OPT = 0;         % Equal to 0 if turned off. Oth., it has the eps_opt value 


if ~isempty(varargin)
    opt_input = varargin{1};
    
    if isfield(opt_input, 'rho')
        rho = opt_input.rho;
    end
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

if length(components) ~= P
    error('Size of components is not P. Please type ''help DADMM_Partial''.');
end
% =========================================================================


% =========================================================================
% Initializations

Sp = zeros(P,1);                % Number of components each node depends on
                                % node p depends on at the pth entry
U = cell(P,1);                  % Cell (one/node) with the dual variable
Diff = cell(P,1);               % Cell (one/node) with 
                                %    Dpl*Xl{p} - sum_{j in neighbs}Xl{j}
                                % for each component l node p depends on
X = cell(P,1);                  % Cell (one/node) with current estimates
X_prev = cell(P,1);             % Cell (one/node) with previous estimates


% P x 1 cell where each entry is a np x 1 cell, where np = Sp(p). The entry
% p x l of neighbs_components is a vector with the indices of the neighbors
% of node p that depend on the component l.
neighbs_components = cell(P,1);


for p = 1 : P
    Sp(p) = length(components{p});
    X{p} = INIT_X*ones(Sp(p),1);
    X_prev{p} = X{p};
    U{p} = zeros(Sp(p),1);
    Diff{p} = zeros(Sp(p),1);
    
    neighbs_components{p} = cell(Sp(p),1);
    neighbs = neighbors{p};
    comp_p = components{p};
    for l = 1 : Sp(p)
        vec_p_l = [];        
        for j = 1 : length(neighbs)
            if sum(components{neighbs(j)}' == comp_p(l)) > 0
                vec_p_l = [vec_p_l ; neighbs(j)]; %#ok<AGROW>
            end
        end
        neighbs_components{p}{l} = vec_p_l;
    end    
end

Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated

                                      
% Struct with the internal variables of DADMM_Partial. 
vars = struct( 'Sp', {Sp}, ...
    'neighbs_components', {neighbs_components}, ...
    'U', {U}, ...
    'Diff', {Diff}, ...
    'X', {X},...
    'X_prev', {X_prev},...
    'rho', {rho}, ...
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
                                
for k = 1 : MAX_ITER    
                        
    vars.X_prev = X;
    [vars, vars_prob] = DADMM_compute_gradient_lagr(vars, vars_prob, vars_network);
    Diff = vars.Diff;
    X = vars.X;
                                                                        
    for p = 1 : P
        % Only iterate if the nodes are still active
        if Stop(p) == 0  
            U{p} = U{p} + rho*Diff{p};
        end
    end
    vars.U = U;
                                              
    % Stopping criterion for the outer loop
    Stop = DADMM_stopping_criterion(Stop, vars_network, vars, vars_prob, k);
    
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




function [vars, vars_prob] = DADMM_compute_gradient_lagr(vars, vars_prob, vars_network)

% This function updates the primal variables according to the Extended ADMM 
% algorithm. It uses the graph coloring to select which nodes to update
% sequentially.

% DADMM_Partial variables
Sp = vars.Sp;
neighbs_components = vars.neighbs_components;
U = vars.U;
Diff = vars.Diff;
X = vars.X;
rho = vars.rho;
Stop = vars.Stop;
EXISTS_X_OPT = vars.EXISTS_X_OPT;
ERROR_FUN = vars.error_fun;

if EXISTS_X_OPT == 1
    x_opt = vars.x_opt;
    error_iterations = vars.error_iterations;
    iter_for_errors = vars.iter_for_errors;
    iterations = vars.iterations;
end


% Network variables
P = vars_network.P;
partition_colors = vars_network.partition_colors;    


% Problem variables (function handler and components)
user_solver = vars_prob.handler;
components = vars_prob.components;


% =========================================================================
% Algorithm

num_colors = length(partition_colors);

for color = 1 : num_colors
    
    X_aux = X;
    for p = partition_colors{color}
        
        if ~Stop(p)              % Only compute if not has not terminated                        
            comp_p = components{p};
            
            v_aux = zeros(Sp(p),1);
            c_aux = zeros(Sp(p),1);
            
            % Determine the sum of the X's of the neighbors for each
            % component l
            for ind_l = 1 : Sp(p)
                l = comp_p(ind_l);
                neighbs_Vl = neighbs_components{p}{ind_l};
                sum_neighbs_l = 0;
                for ind_n = 1 : length(neighbs_Vl)
                    neighbor = neighbs_Vl(ind_n);
                    index_nei = (components{neighbor} == l);
                    if sum(index_nei) > 1
                        error('There were double indices in components.');
                    end
                    sum_neighbs_l = sum_neighbs_l + X{neighbor}(index_nei);
                end
                v_aux(ind_l) = sum_neighbs_l;
                c_aux(ind_l) = length(neighbs_Vl);
            end
            
            
            % Solving each problem in the alternating direction minimization
            
            v = U{p} - rho*v_aux;
            c = (rho/2)*c_aux;
            [X_aux{p}, vars_prob] = user_solver(p, v, c, X, vars_prob);
        end
    end
    X = X_aux;
    
end
% =========================================================================

% =========================================================================
% Output

for p = 1 : P
    
    comp_p = components{p};    
    Diff_p = zeros(Sp(p),1);
    
    for ind_l = 1 : Sp(p)
        l = comp_p(ind_l);
        neighbs_Vl = neighbs_components{p}{ind_l};
        index_l_in_p = (components{p} == l);
        for ind_n = 1 : length(neighbs_Vl)
            neighbor = neighbs_Vl(ind_n);
            index_nei = (components{neighbor} == l);
            Diff_p(ind_l) = Diff_p(ind_l) + ( X{p}(index_l_in_p) - X{neighbor}(index_nei) );
        end
    end
    
    Diff{p} = Diff_p;
end

vars.Diff = Diff;
vars.X = X;

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
    error_iterations = [error_iterations, new_error];          %figure(10);clf;semilogy(error_iterations);drawnow;
    ind_nonfilled = (iter_for_errors(:,2) == Inf);
    ind_lesserror = (new_error < iter_for_errors(:,1));
    intersect = ind_nonfilled & ind_lesserror;
    iter_for_errors(intersect,2) = iterations + 1;
        
    vars.error_iterations = error_iterations;
    vars.iter_for_errors = iter_for_errors;
    vars.iterations = iterations + 1;
end
% =========================================================================

end




function [Stop] = DADMM_stopping_criterion(Stop, vars_network, vars, vars_prob, k)


if k <= 2 || vars.TURN_OFF_EPS == 1
    return
end

% Network variables
P = vars_network.P;

% Problem variables
components = vars_prob.components;

% D-ADMM variables
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
