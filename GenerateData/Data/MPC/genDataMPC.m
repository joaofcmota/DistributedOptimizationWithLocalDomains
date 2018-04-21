% Generates data for the MPC problem, given
%
%   T:   horizon of the MPC problem
%
%   n_p: dimension of the state at each subsystem
%
%   m_p: dimension of the input at each subsystem
%
%
% The variable has the following format:
%
% var = [ u_1 ;     % each with dim = m_p*T
%         u_2 ;
%          ...
%         u_P  ]
%
%
% This script also creates the factor graph for the problem.

% =========================================================================
% Filenames

% Network  (Power grid)
PATH_NETWORK = '../../Networks/MPC/ProcessedNetwork/';

FILENAME_NETWORK = 'Barabasi_P_100.mat';
%FILENAME_NETWORK = 'Barabasi_P_2000.mat';
%FILENAME_NETWORK = 'Network_PowerGrid.mat';


% Base to file for storing the data (.mat will be added)
FILENAME_OUTPUT_BASE = 'Data/Barab';            
%FILENAME_OUTPUT_BASE = 'Data/MPC_PGrid';
% =========================================================================


% =========================================================================
% Parameters

STABLE = 0;      % If True, each subsystem will be designed to be stable.

T = 5;           % Time horizon

n_p = 3;         % Dimension of each state

m_p = 1;         % Dimension of each control
% =========================================================================


% =========================================================================
% Load Network

filename_network = [PATH_NETWORK , FILENAME_NETWORK];
load(filename_network);

P = Network.P;
Neighbors = Network.Neighbors;
Dp = Network.Degrees;
% =========================================================================


% =========================================================================
% Create Factor Graph and Omega
%
% Omega is a P x 1 cell where the pth entry contains the neighbors of node
% p plus node p, all ordered.

components = cell(P,1);
Omega = cell(P,1);

ind_u = cell(P,1);               % Given a node p gives the indices of u_p
for p = 1 : P
    ind_u{p} = 1 + (p-1)*m_p*T : p*m_p*T;
end

for p = 1 : P
    neighbs_p = Neighbors{p};    % Neighbors of node p (only x's)
    Dp_p = Dp(p);
    
    Omega{p} = sort( [Neighbors{p} , p] );
    
    comp_p = zeros(m_p*T*(Dp_p+1) , 1);
    
    for j = 1 : Dp_p
        ind = 1 + (j-1)*m_p*T : j*m_p*T;
        comp_p(ind) = ind_u{neighbs_p(j)};
    end
            
    % Add component u of node p
    ind = 1 + Dp_p*m_p*T : (Dp_p+1)*m_p*T;
    comp_p(ind) = ind_u{p}; 
         
    % Needs sorting (read DADMMp documentation)
    components{p} = sort(comp_p);
end
% =========================================================================


% =========================================================================
% First step: generate the following matrices (all cells Px1):
%
%   A_pp:    dynamics; the pth entry contains the internal matrix of
%            subsystem p
%
%   B_pj:    influence of input u_j on state x_p (j is as in components)
%
%   Q_bar:   state cost; each entry contains a vector representing the 
%            diagonal of a (diagonal) matrix 
%
%   Q_f_bar: state cost final; same as Q
%
%   R_bar:   cost of inputs; same as Q
%
%   x0:      initial/measured state; each entry is a vector


A_pp    = cell(P,1);
B_pj    = cell(P,1);
Q_bar   = cell(P,1);
Q_f_bar = cell(P,1);
R_bar   = cell(P,1);
x0      = cell(P,1);

for p = 1 : P
    
    % Matrix A_pp
    if STABLE == 1                             % Generate stable subsystems
        aux_A_pp = randn(n_p,n_p);
        [aux_Q,aux_R] = qr(aux_A_pp);
        eigenvals = 2*rand(n_p,1) - ones(n_p,1);       % each eig in [-1,1]
        A_pp{p} = aux_Q*diag(eigenvals)*aux_Q';
    else
        A_pp{p} = randn(n_p,n_p);
    end
            
    
    % Matrices B_pj (ordered according to components)
    Dp_p = Dp(p);
    B_pj{p} = cell(Dp_p+1,1);
    
    for j = 1 : Dp_p + 1   
        B_pj{p}{j} = randn(n_p,m_p);
    end        
            
    % Matrices Q_bar, Q_f_bar, and R_bar
    Q_bar{p}   = 9*rand(n_p,1) + ones(n_p,1);
    Q_f_bar{p} = 9*rand(n_p,1) + ones(n_p,1);
    R_bar{p}   = 9*rand(m_p,1) + ones(m_p,1);
    
    % State initialization s0
    x0{p} = randn(n_p,1);
end
% =========================================================================



% =========================================================================
% Second step: generate the following matrices (all cells Px1); see the 
%               documentation 
%
%               R_p, Q_p, C_p, D_p0, E_p, w_p
%

R_p  = cell(P,1);
Q_p  = cell(P,1);
C_p  = cell(P,1);
D_p0 = cell(P,1); 
E_p  = cell(P,1);
w_p  = cell(P,1);

Id_T = ones(T,1);

for p = 1 : P
   
    % *********************************************************************
    % Matrices R_p and Q_p (diagonal matrices, so we just store the diag)    
          
    R_p{p} = kron(Id_T, R_bar{p});
    Q_p{p} = [kron(Id_T, Q_bar{p}) ; Q_f_bar{p}];
    % *********************************************************************
            
    Dp_p      = Dp(p);                     % Degree of node p
    neighbs_p = Neighbors{p};              % Neighbors of node p
    B_pj_cat  = zeros(n_p,m_p*(Dp_p+1));   % Concatenation of matrices B_pj
    
    % Perform concatenation of B_pj
    for j = 1 : Dp_p + 1        
        B_pj_cat(: , 1 + (j-1)*m_p : j*m_p) = B_pj{p}{j};
    end
        
    % *********************************************************************
    % Matrices C_p and D_p0
    C_p_aux  = zeros( n_p*(T+1) , m_p*T*(Dp_p+1) );
    D_p0_aux = zeros( n_p*(T+1) , 1);
    
    generating_vec_C = zeros(n_p,m_p*(Dp_p+1)*T);    
    
    powers_of_A_pp = eye(n_p);
    
    D_p0_aux(1:n_p) = powers_of_A_pp * x0{p};
    
    for t = 1 : T     % Note: first line of C_p and E_p contains only zeros
                                     
        generating_vec_C = circshift(generating_vec_C, [0 m_p*(Dp_p+1)]);
        generating_vec_C(: , 1:m_p*(Dp_p+1)) = powers_of_A_pp*B_pj_cat;
         
        powers_of_A_pp = powers_of_A_pp*A_pp{p};
        
        C_p_aux(1 + t*n_p : (t+1)*n_p , :) = generating_vec_C;
        
        D_p0_aux(1 + t*n_p : (t+1)*n_p) = powers_of_A_pp * x0{p};
    end
    
    C_p{p}  = C_p_aux;
    D_p0{p} = D_p0_aux;
    % *********************************************************************
    
    
    % *********************************************************************
    % Finally, matrix E_p and vector w_p
    
    Omega_p = Omega{p};
    
    position_of_p = find(Omega_p == p);
    
    indices_of_p = 1 + (position_of_p - 1)*m_p*T : position_of_p*m_p*T;
    
    C_p_T_Q_p = C_p{p}'*diag(Q_p{p});                  % Will be used twice
    
    E_p_aux = C_p_T_Q_p*C_p{p};
    
    E_p_aux(indices_of_p, indices_of_p) = ...
        E_p_aux(indices_of_p, indices_of_p) + diag(R_p{p});
    
    E_p{p} = E_p_aux;
    
    w_p{p} = 2*C_p_T_Q_p*D_p0{p};
    % *********************************************************************
end


% =========================================================================
% Compute the solution

% Construct matrices E_global and w_global

dim_variable = P*T*m_p;
E_global = zeros(dim_variable , dim_variable);
w_global = zeros(dim_variable , 1);

for p = 1 : P
    ind_p = components{p};
    E_global(ind_p, ind_p) = E_global(ind_p, ind_p) + E_p{p};
    w_global(ind_p)        = w_global(ind_p)        + w_p{p};
end

solution = -0.5*(E_global \ w_global);
% =========================================================================


% =========================================================================
% Compute Lipschitz constant
Lipschitz = 2*max(eigs(E_global));
% =========================================================================


% =========================================================================
% Save

% Filename
if STABLE == 1
    filename_ending = [FILENAME_OUTPUT_BASE, '_P_', num2str(P), '_STABLE'];
else
    filename_ending = [FILENAME_OUTPUT_BASE, '_P_', num2str(P), '_UNSTABLE'];
end

filename_data = [filename_ending, '.mat'];

save(filename_data, 'R_p', 'Q_p', 'C_p', 'D_p0', 'E_p', 'w_p', ...
    'components', 'Omega', 'T', 'n_p', 'm_p', 'solution', 'Lipschitz', ...
    'A_pp', 'B_pj', 'Q_bar', 'Q_f_bar', 'R_bar', 'x0');
% =========================================================================

% =========================================================================
% Save Statistics to text file

filename_stat = [filename_ending, '.txt'];

all_text = [...
'--------------------------------------------------------------\n' ...
'Statistics about the Data in: ', filename_data, '\n\n' ...
'Path to Network: ../', filename_network, '\n\n' ...
'Parameters defined by the user: \n' ...
'  STABLE:    ', num2str(STABLE), '\n' ...
'  T:         ', num2str(T), '\n' ...
'  n_p:       ', num2str(n_p), '\n' ...
'  m_p:       ', num2str(m_p), '\n', ...
'  Lipschitz: ', num2str(Lipschitz), '\n\n', ...
'--------------------------------------------------------------\n'];

fid = fopen(filename_stat, 'w');
fprintf(fid, all_text);
fclose(fid);
% =========================================================================






