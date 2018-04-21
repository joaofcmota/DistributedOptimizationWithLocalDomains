#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 23:56:33 2012

Creates random bipartite networks, using the Erdos-Renyi model (see 
the networkx documentation).

This script should have 4 inputs:
  C: number of nodes in one group
  P: number of nodes in the other group
  prob: number between 0 and 1
  Filename: where the data is saved
"""



# main function
def main():

  import sys 
  import networkx as nx
  from networkx.algorithms import bipartite
    
  if len(sys.argv) == 5:
    C_in = sys.argv[1]
    P_in = sys.argv[2]
    prob_in = sys.argv[3]
    Filename = sys.argv[4]
  else:
    print 'Error: there should be 3 input arguments in gen_bipartite_erdos'
    return
  
  # check for input errors
  #if 
  #print 'C = ' + str(C)
  #print 'P = ' + str(P)
  #print 'par = ' + str(par)
  
  C = int(float(C_in))
  P = int(float(P_in))
  prob = float(prob_in)

  G = nx.generators.bipartite_random_graph(C,P,prob)
  Adj = nx.adj_matrix(G)
  if nx.is_connected(G):
    conn = 1
    # get a coloring scheme  
    colors = bipartite.color(G)
    diameter = nx.diameter(G)
  else:
    conn = 0
    diameter = 0
    print 'Error: generated network was not connected. Try increasing prob'
  
  # get a coloring scheme  
  colors = bipartite.color(G)

  number_of_nodes = C + P
  
  # write to file
  FILE=open(Filename, "w")
  FILE.write('conn = '+str(conn)+'\n')
  FILE.write('C = '+str(C)+'\n')
  FILE.write('P = '+str(P)+'\n')
  FILE.write('diameter = '+str(diameter)+'\n')
  FILE.write('Adj =\n')
  
  for node1 in range(number_of_nodes):
    for node2 in range(number_of_nodes):
      FILE.write(str(int(Adj[node1,node2]))+' ')  
    FILE.write('\n')
  
  if conn == 1:
    FILE.write('colors = \n')
    for node1 in range(number_of_nodes):
      FILE.write(str(colors[node1])+' ')
    
  FILE.write('\n')  
  FILE.close()
  
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
