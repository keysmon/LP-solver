Hang Ruan
V00923058

Description:
LP solver which reads a text-based representation of a linear program from standard input and outputs the result of solving the LP (including whether the LP turned out to be infeasible or unbounded). 


Library used:
1. sys: to get standard input from user
2. numpy: to store list and basic array operation
	

Guide to run:
1. cat input.txt | python simplex_alg.py
or 
1. python simplex_alg.py 
2. copy and paste LP into stdin

Solver architecture:
It uses dual simplex method to turn any initial infeasible dictionary feasible. Next, the program uses steepest edge pivot rule and falling back on Bland's rule when three degenerate pivot occur. The pivot operation are performed in dictionary form.

Extra features:
1. steepest edge pivot rule, Bland's rule as pivot strategy.
2. Primal-Dual methods to solve initially infeasible LPs using a two phased-dual method.
