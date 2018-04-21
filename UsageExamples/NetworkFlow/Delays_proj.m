function [proj, vars_BB] = Delays_proj(x, vars_BB)

const = vars_BB.const;
d = vars_BB.d;
capacities = vars_BB.capacities;
NetFlowSolver_proj = vars_BB.NetFlowSolver_proj;

proj = NetFlowSolver_proj(x, const, d, capacities);




