Particle Swarm optimization to train models written in PySB.

Examples include train_pso_earm.py. This model uses a basic PSO implemented with DEAP. 

ParticleBased Iterative Nonlinear Optimization (PINO, clever I know).
This uses multiple swarms of particles in parameter space. Each swarm returns a vector towards its local minimum. 
The next step is to cluster these vectors and reinitiate placement of the new sets of swarms, now into this smaller, more guided parameter space.
This is done until the vector returned from particles within a cutoff distance are all zero. This would imply that the local minimum is within this population.

So far, I have the multiple swarm implemented. I have two examples contained in particle_swarm_multi_layer.py. The images shows how the particles all start out in equil distance regions of parameter space. They then move in the direction of the best region. 

TODO
1.) Make general functions that allow the number of swarms, number of particles per swarm, and limits on swarm movements all proportional to the sampling space.
2.) Make particle speed a function of each parameter.
3.) Implement PyPar ?


