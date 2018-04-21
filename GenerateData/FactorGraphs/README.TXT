
Contains functions for creating a connected factor graph. The networks used here
can be found in the folder '/GlobalVar/GenerateData/Networks/'.

The algorithm in the function 'CreateVarPartialConnected.m' works as follows.

Given a network with P nodes and a variable with n components,

  For l = 1 : n

      select node p randomly, and node p will depend on x_l

      fringe = {neighbors_of_p}
 
      For dep = 1 : depth

          select one node j randomly in the fringe

          node j will depend on x_l

          add neighbors_of_j to the fringe, and remove j from the fringe

      EndFor
  
  EndFor

  If at the end of this algorithm there is at least one node that does not
  depend on any component, the algorithm is repeated.


The files in this folder are

  CreateVarPartialConnected: implements the above algorithm

  CreateVarScript: calls CreateVarPartialConnected for the networks
                   in /GlobalVar/GenerateData/Networks/ and puts them
                   in the folder /GeneratedFGs.

  
