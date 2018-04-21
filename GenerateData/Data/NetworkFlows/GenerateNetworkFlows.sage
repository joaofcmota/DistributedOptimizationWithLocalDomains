#!/usr/lib/sagemath/sage -python

# Creates a random instance of the following network flow problem:
#
#     minimize    \sum_{(i,j) in edges} x_ij/(c_ij - x_ij)
#     subject to  B*x == vector_of_flows
#
# where the problem is associated to a Digraph (obtained from an
# undirected network specified by the user), c_ij is the capacity
# of edge ij, and B is the node arc incidence matrix. 
# This script also solves a multicommodity flow problem, from which
# the data was generated and also solves the problem above.


import numpy as np
from scipy.io import savemat              # To save matrices to Matlab
from scipy.io import loadmat              # To save matrices to Matlab
from scipy import *
import cvxopt as cx                       # To solve our problem
import cvxopt.solvers as sol
from sage.all import *


# ******************************************************************
# Parameters: 

# Number of commodities to be generated
NUM_COMMODITIES = 100

# Network directory
DIRECTORY_NETWORKS = '../../Networks/MPC/ProcessedNetwork/'

# Network name 
NAME_NETWORK = 'Barabasi_P_2000.mat'

# Directory to save the data
DIRECTORY_DATA_OUTPUT = 'Results/'

# Will be used as seed for everything random
random_seed = 1234

# When assigning flows, if we select (randomly) a node with an empty
# reachable set, we select another node. The maximum number of 
# trials is MAX_TRIALS
MAX_TRIALS = 1000

# The flow intensity for each selected node origin and node destination
# is drawn randomly from the set of possible_capacities, divided by
# DIV_FACTOR (has to be real)
DIV_FACTOR = 100.0      
# ******************************************************************


# ******************************************************************
# Preprocessing

# Filename of networks file 
FILENAME_NETWORK = DIRECTORY_NETWORKS + NAME_NETWORK

# Filename for the output
FILENAME_DATA_OUTPUT = DIRECTORY_DATA_OUTPUT + 'NFData_' + NAME_NETWORK

# Load network
matdata = loadmat(FILENAME_NETWORK)
Network = matdata['Network']
Adj_mat = Network['Adj']

# Create a sage graph (undirected)
G_undirected = Graph(Adj_mat[0][0])

num_nodes = G_undirected.num_verts()
num_edges = G_undirected.num_edges()
# ******************************************************************


# ******************************************************************
# Create a digraph from the directed graph by assigning a random 
# direction to each edge. At the same time, assign random capacities

# Possible capacities (GBit/s) for each link and respective probability
possible_capacities = {'10' : 0.2,
                       '20' : 0.2,
                       '30' : 0.2,
                       '40' : 0.2,
                       '50' : 0.1,
                       '100': 0.1 }

# Declare probability distributions
Uniform = RealDistribution('uniform', [0,1], seed=random_seed)
Capacities_distribution = GeneralDiscreteDistribution(possible_capacities.values(),
                                                      seed=random_seed)


# Create list of edges for the digraph 
list_edges = list()

list_edges_undir = G_undirected.edges()
for edge in list_edges_undir:
    node1 = edge[0]
    node2 = edge[1]

    # Determine capacity randomly
    cap = Capacities_distribution.get_random_element()

    if Uniform.get_random_element() >= 0.5:
        list_edges.append((node1,node2, int(possible_capacities.keys()[cap])))
    else:
        list_edges.append((node2,node1, int(possible_capacities.keys()[cap])))
        
G = DiGraph(list_edges)        
# ******************************************************************


# ******************************************************************
# Compute the distances between all pairs of nodes (some will be
# infinite)

distances = G.distance_all_pairs() 

# Compute the reachable set of nodes for each node, i.e., the nodes for which
# their distance to node p is less than infinite
reachable_sets = list()

for i in range(num_nodes):
    dist_i = distances[i]
    reachable_set_i = [node for node in range(num_nodes) if dist_i[node] != Infinity and node != i]
    reachable_sets.append(reachable_set_i)
# ******************************************************************


# ******************************************************************
# Create flows randomly: for each commodity, we select a node p that
# has a non-empty reachable set and pick randomly a node j in its 
# reachable set; then that commodity has to go from p to j

vector_of_flows = np.zeros(num_nodes)

# Triplet (node_in, node_out, intensity)
node_out_node_in_intensity = list()

# List of selected nodes
list_selected_nodes = list()

# Probability distribution
vec_ones = [1.0/num_nodes for i in range(num_nodes)]
nodes_distr = GeneralDiscreteDistribution(vec_ones, seed=random_seed)

for comm in range(NUM_COMMODITIES):
    for trial in range(MAX_TRIALS):
        node_p = nodes_distr.get_random_element()
        reachable_set_node_p = reachable_sets[node_p]
        if len(reachable_set_node_p) > 0 and list_selected_nodes.count(node_p) == 0:

            # Select a node at random from the reachable set
            len_reach = len(reachable_set_node_p)
            vec_ones = [1.0/len_reach for i in range(len_reach)]
            reach_distr = GeneralDiscreteDistribution(vec_ones, seed=random_seed)
            node_j = reachable_set_node_p[reach_distr.get_random_element()]

            flow_aux = Capacities_distribution.get_random_element()
            flow = int(possible_capacities.keys()[flow_aux])/DIV_FACTOR            

            break

    if trial == MAX_TRIALS - 1:
        print 'Error: MAX_TRIALS reached'

    vector_of_flows[node_p] = vector_of_flows[node_p] - flow
    vector_of_flows[node_j] = vector_of_flows[node_j] + flow

    node_out_node_in_intensity.append((node_p,node_j, flow))

    list_selected_nodes.append(node_p)
# ******************************************************************



# ******************************************************************
# Compute solution of a multicommodity flow problem

solution_mcfp_aux = G.multicommodity_flow(node_out_node_in_intensity, 
        integer=False, use_edge_labels=True)

solution_mcfp = np.zeros((NUM_COMMODITIES,), dtype=np.object)
for i in range(NUM_COMMODITIES):
    solution_mcfp[i] = solution_mcfp_aux[i].edges()
# ******************************************************************



# ******************************************************************
# Compute solution to our problem with cvxopt

vec_capacities = np.array(G.edge_labels())
vec_capacities = cx.matrix(vec_capacities, tc='d')
node_arc_incidence = cx.matrix(np.array(G.incidence_matrix()), tc='d')

def F(x=None, z=None):
    if x is None: 
        return 0, cx.matrix(vec_capacities/2.0)
    caps_x = cx.matrix(vec_capacities - x)
    if min(caps_x) <= 0.0:
        return None
    f = sum(cx.div(x,caps_x))
    Df = cx.div(vec_capacities, caps_x**2).T
    if z is None:
        return f, Df
    H = cx.spdiag(2.0*z[0]*cx.div(vec_capacities, caps_x**3))
    return f, Df, H

G_mat = cx.matrix(np.concatenate(( -np.identity(num_edges), np.identity(num_edges) )), tc='d')
h_mat = cx.matrix(np.concatenate((np.zeros((num_edges,1)) , vec_capacities)))

# NOTE: the linear system 'node_arc_incidence*x == vector_of_flows' has one dependent
# equation; so, we remove the last one.
b_mat = cx.matrix(vector_of_flows[0:num_nodes-1])
A_mat = cx.matrix(node_arc_incidence[0:num_nodes-1,:], tc='d')

solution_prob = sol.cp(F, G=G_mat, h=h_mat, A=A_mat, b=b_mat)

solution = solution_prob['x']

if solution_prob['status'] != 'optimal':
    print '\nPlease repeat the experiment with a different seed/parameters: optimal solution not found\n'
# ******************************************************************



# ******************************************************************
# Save data to a .mat file

savemat(FILENAME_DATA_OUTPUT, {'capacities': vec_capacities, 
                               'incidence_matrix': node_arc_incidence, 
                               'flows': vector_of_flows,
                               'solution_mcfp': solution_mcfp,
                               'solution': solution,
                               'node_out_node_in_intensity': node_out_node_in_intensity
                               })
# ******************************************************************


# ******************************************************************
# Print to stdout characteristics of the data

print '\n\n'
print '--------------------------------------------------'
print 'Data parameters and statistics:\n'
print '\n'
print 'Network file: ' + DIRECTORY_NETWORKS + NAME_NETWORK
print '\n'
print 'Number of commodities:   ' + str(NUM_COMMODITIES)
print 'Random seed:             ' + str(random_seed)
print 'Division factor:         ' + str(DIV_FACTOR)
print 'Possible capacities: '
for i in range(len(possible_capacities)):
    print possible_capacities.keys()[i] + ' : ' + str(possible_capacities.values()[i])
print '\n'
print 'Lengths of paths in solution of multicommodity flow: '
for i in range(NUM_COMMODITIES):
    print str(len(solution_mcfp[i]))
print '\n'
# ******************************************************************


# Plots
# G_undirected.show(layout='spring', save_pos=True)
# pos_graph = G_undirected.get_pos()
# G.show(pos=pos_graph, edge_labels=True)


