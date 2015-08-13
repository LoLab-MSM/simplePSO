Particle Swarm optimization to train models written in PySB.

Examples include example_earm.py. This model uses a basic PSO implemented with DEAP. 

Quickly, PSO can be set up as follows
from pso import PSO
pso = PSO()
# The function must accept a particle, which should be a parameter vector in log space (See example_earm.py for an example). It returns a single scalar to be added to a tuple.

pso.set_cost_function(cost_function)

# This is not required if the solver is defined as a global variable. 
pso.set_solver(solver)

# The initial starting position. This value will be used to initate the population of particles. 
pso.set_start_position(xnominal)

# bounds in log space, If given single value upper and lower will be calculated based on starting postion +/- range. Otherwise can set lower and upper as a vector of len(xnominal)
pso.set_bounds(range, lower,upper)

# set maximum and minimum speed of particles. 
pso.set_speed(-0.5,0.5)

# To run you need to set number of particles and number of iterations.
pso.run(n_particles,n_iteration)

# After it has completed you can access the best particle
pso.best


You can view PSO_UML.png to explore more attributes and functions.


You can email me at james.ch.pino@gmail.com for any questions or comments.




