# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:34:54 2014

@author: pinojc
"""


import numpy as np

 
#class creation
class Particle(object):
  def __init__(self, swarm,max_vel , min_vel , max_pos, min_pos,dim,obj_function,
               best = None, position = None, 
               velocity = None, fitness = None,
               ):
    #initilize a particle by using the following  
    self.swarm = swarm
    self.obj_function = obj_function
    self.dim =  dim
    self.max_vel =  max_vel
    self.min_vel =  min_vel
    self.max_pos =  max_pos
    self.min_pos =  min_pos
    self.position = self._return_start() if position is None else position
    self.best     = self.position      if best is None     else best
    self.fitness  = self.obj_function(self.position)     if fitness is None     else fitness
    self.velocity = self._random_velocity() if velocity is None else velocity
    
  def _return_start(self):
    start = np.random.uniform(self.min_pos,self.max_pos,self.dim)
    fitness = self.obj_function(start)    
    for _ in range(10):
      tmp = np.random.uniform(self.min_pos,self.max_pos,self.dim)
      tmp2 = self.obj_function(tmp)
      if tmp2 < fitness:
        start = tmp
        fitness = tmp2
    return start
  
  def _random_velocity(self):
    # Get a random velocity vector of values your MIN and MAX scaled by 1/10th 
    vel = np.random.uniform(self.min_vel,self.min_vel,self.dim)
    return vel
    
  def step(self):


    # Step the algorithm for one timestep
    self.velocity = self.velocity + \
                      np.random.uniform(0,1)*(self.best - self.position)
    +                 np.random.uniform(0,1)*(self.swarm.globalbest - self.position)
    
    self.velocity[self.velocity > self.max_vel ] = self.max_vel
    self.velocity[self.velocity < self.min_vel ] = self.min_vel
    #print np.sum(self.velocity)
    self.position = self.position + self.velocity
    self.position[self.position > self.max_pos] = self.max_pos
    self.position[self.position < self.min_pos] = self.min_pos


    # Check to see if this new position is fitter than the previous best
    
    tmp = self.obj_function(self.position)
    if tmp < self.fitness:
      self.best = self.position
      self.fitness = tmp
    else:
      self.best=self.best
    return 
  def return_position(self):
      return self.best


 
class Swarm(object):
  def __init__(self,swarm_size ):
    # Get the number of arguments our function takes
    self.swarm_size  = swarm_size
    self.particles = [None] * swarm_size
    self.globalbest = None
    self.time = 0
    self.globalFitness = None
  
  def add(self, particle,i):
    self.particles[i] = particle
    if self.particles[-1] != None:
        self._update_best()
  
  def _update_best(self):
    if self.globalFitness == None:
      tmp_best = None
    else:
      tmp_best = self.globalFitness
    for particle in self.particles:
      tmp_fit =  particle.fitness
      if tmp_best is None or tmp_fit < tmp_best:
        print 'improved global',tmp_fit
        tmp_best = tmp_fit
        current_best = particle.position
    current_best_fitness = tmp_best
    if self.globalFitness is None or current_best_fitness < self.globalFitness:
      self.globalbest = current_best
      self.globalFitness = current_best_fitness
  
  def evaluate_fitness(self):
    self.time += 1
    for particle in self.particles:
      particle.step()
    self._update_best()


# main
# Generate a swarm
def PSO(swarm_size,iterations,best,obj_function,max_vel,
        min_vel,max_pos,min_pos,dim):
  swarm = Swarm(swarm_size)
  swarm_size = swarm_size

  for j in xrange(swarm_size):
    swarm.add(Particle(swarm,position=None,best=best,obj_function=obj_function,
            max_vel = max_vel , min_vel = min_vel,
            max_pos = max_pos, min_pos = min_pos,  dim = dim ),j)
    # Run the algorithm a few times
    num_steps = iterations
  for _ in range(num_steps):
    swarm.evaluate_fitness()
    if _ % 10 == 0:
      print _

  print "After ",swarm.time,"iterations with "+str(swarm_size)+\
  " particles, PSO converged at with a value of ",swarm.globalFitness
  population = [None] * swarm_size
  count = 0
  for particle in swarm.particles:
      population[count] = particle.return_position()
      count+=1
  return swarm.globalbest,swarm.globalFitness, population