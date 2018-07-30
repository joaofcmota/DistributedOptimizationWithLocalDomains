# Distributed Optimization With Local Domains: Applications in MPC and Network Flows

Code to replicate the experiments in

**[Distributed Optimization With Local Domains: Applications in MPC and Network
Flows](http://dx.doi.org/10.1109/TAC.2014.2365686)**.  
  J. F. C. Mota, J. M. F. Xavier, P. M. Q. Aguiar, M. Püschel.  
  IEEE Transactions on Automatic Control, Vol. 60, No. 7, pp. 2004-2009, 2015.  
  [link](http://dx.doi.org/10.1109/TAC.2014.2365686),
  [arXiv](http://arxiv.org/abs/1305.1885) 

The main algorithm is `DADMM_Partial.m`.

All code is implemented in Matlab. Although all algorithms are 
distributed, the code is centralized and runs sequentially. The purpose is to 
simulate distributed algorithms.

## Organization

* `DADMM_Partial.m`:
	solves the optimization problem

	<p align="center"><img src="svgs/6ec3283e6e3265a8f069c5185acc3866.svg" align=middle width=488.24269395000005pt height=23.866193999999997pt/></p>

	where the optimization variable has <img src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg" align=middle width=9.86687624999999pt height=14.15524440000002pt/> components <img src="svgs/b7dec1314efeadd1815bf346b3472200.svg" align=middle width=131.5961955pt height=24.65753399999998pt/>, and
	each set <img src="svgs/2da32a437e1e022c68095bcb359d92a2.svg" align=middle width=16.85640329999999pt height=22.465723500000017pt/> indexes the components of <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg" align=middle width=9.39498779999999pt height=14.15524440000002pt/> that node <img src="svgs/2ec6e630f199f589a2402fdf3e0289d5.svg" align=middle width=8.270567249999992pt height=14.15524440000002pt/> depends on. See the 
	documentation in `DADMM_Partial`, and also the above paper [1] for terminology.

* UsageExamples: simple examples illustrating how to use `DAMM_Partial.m`. 

* CompareAlgs: 
  reproduces Figures 4 and 5 in the paper.

	
* GenerateData: contains data used in the experiments and also contains code
				  to generate other types of data. Some data is generated in Python/Sage.

* KekatosADMM: our implementation of the algorithm in

  **[Distributed Robust Power System State Estimation](
  https://doi.org/10.1109/TPWRS.2012.2219629)**.  
  V. Kekatos, G. B. Giannakis.  
  IEEE Transactions on Power Systems, Vol. 28, No. 2, 2013.  
  [link](https://doi.org/10.1109/TPWRS.2012.2219629),
  [arXiv](https://arxiv.org/abs/1204.0991) 

  Note that the algorithm proposed in this paper was designed to solve (1) with
  a star-shaped variable, but we generalized it to operate with any connected
  variable.


* SpecialPurposeAlgs: 
  our implementation of the algorithms in the references below, which solve (1)
  only in the case of a star-shaped variable.
	
  **[Distributed Optimization and Statistical Learning via the Alternating
  Direction Method of Multipliers](
  https://www.nowpublishers.com/article/Details/MAL-016)**.  
  S. Boyd, N. Parikh, E. Chu, B. Peleato, J. Eckstein.  
  Foundations and Trends on Machine Learning, Vol. 3, No. 4, 2010.  
  *(section 7.2)*.  
  [link](https://www.nowpublishers.com/article/Details/MAL-016),
  [authors' version](https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf)

  **[Introductory Lectures on Convex Optimization: A Basic Course](
  https://doi.org/10.1007/978-1-4419-8853-9)**.  
  Y. Nesterov.  
  Springer Science+Business Media New York, 2004.  
  [link](https://doi.org/10.1007/978-1-4419-8853-9), 
  [Google
  Books](https://books.google.co.uk/books?hl=en&lr=&id=2-ElBQAAQBAJ&oi=fnd&pg=PA1&dq=Introductory+Lectures+on+Convex+Optimization:+A+Basic+Course&ots=wlrO5qqckz&sig=2LOforFMisXArmF_2AxYg6LvXXA#v=onepage&q=Introductory%20Lectures%20on%20Convex%20Optimization%3A%20A%20Basic%20Course&f=false)

---

License: [ GPLv3 ]( https://www.gnu.org/licenses/gpl-3.0.en.html )
	
	


