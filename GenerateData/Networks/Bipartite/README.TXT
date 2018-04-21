
We have the following files:

  generate_bipartite_erdos.py: Python script that creates random bipartite 
                               networks, using the Erdos-Renyi model 
                               (networkx)

  generate_series_bipartite.m: calls generate_bipartite_erdos.py for several
                               network sizes. It also creates automatically
                               the structure FactorGraphs.

  Bipartite_Nets.mat: contains bipartite networks in 'Networks' and the 
                      corresponding 'FactorGraphs': each central node
                      is considered the center of a cluster.
