----------------------------------------------------------------------------------
Distributed Optimization With Local Domains: Applications in MPC and Network Flows
Copyright (C) 2018 João Mota
----------------------------------------------------------------------------------

Thank you for downloading this software. 

The purpose of this package is to make available a Matlab implementation of the 
algorithm and the reproducibility of the experiments in the paper


  [1]  J. Mota, J. Xavier, P. Aguiar, M. Püschel,
	   Distributed Optimization With Local Domains: Applications in MPC and 
       Network Flows, to appear in IEEE Transactions on Automatic Control, 2015


The implementation of the proposed algorithm is the the file DADMM_Partial.m, in 
the main directory.

All code is implemented in Matlab (GNU/Octave). Although all algorithms are 
distributed, the code is centralized and runs sequentially. The purpose is to 
simulate distributed algorithms.

-------------------------------------------------------------------------------
This software is distributed under the GNU General Public License; see gpl.txt

This package is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
-------------------------------------------------------------------------------


Main algorithm: DADMM_Partial.m

	It solves the optimization problem

    	minimize  f1(x_S1) + f2(x_S2) + ... + fP(x_SP)         (1)
	       x

	where the optimization variable has n components x = (x1, x2, ..., xn), and
	each set Sp indexes the components of x that node p depends on. See the 
	documentation by typing help DADMM_Partial. And see [1] for terminology.


Description of each folder:

	UsageExamples: Contains simple examples illustrating how to use 
				   DAMM_Partial. 

	Compare Algs: Code to reproduce Figures 4 and 5 in [1].

	
	GenerateData: Contains data used in the experiments and also contains code
				  to generate other types of data. Some of the data is 
				  generated in Python/Sage.

	KekatosADMM: Our implementation of the algorithm in

		     	 [2]  V. Kekatos, G. Giannakis, Distributed Robust Power System
					  State Estimation, IEEE Transactions on Power Systems, 
 					  Vol. 28, No. 2, 2013

			 	 Note that the algorithm in [2] was proposed to solve (1) with
				 a star-shaped variable, but we generalized it to operate with 
				 any connected variable.


	SpecialPurposeAlgs: Our implementation of the algorithms in [3] and [4], 
						which solve (1) only in the case of a star-shaped 
						variable.
	
	                    [3]  S. Boyd, N. Parikh, E. Chu, B. Peleato, J. 
                             Eckstein, Distributed Optimization and Statistical
                             Learning via the Alternating Direction Method of 
                             Multipliers, Foundations and Trends on Machine 
                             Learning, Vol. 3, No. 4, 2010 (section 7.2)

                    	[4]  Y. Nesterov, Introductory Lectures on Convex 
                             Optimization: A Basic Course, Kluwer Academic 
                             Publishers, 2003


	
	


