# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""

from deap import base
from deap import creator
from deap import tools
import numpy as np
import multiprocessing


class PSO():
    def __init__(self,cost_function=None,start=None,save_sampled=False,method = 'single_min'):

        
        self.cost_function = cost_function
        self.save_sampled = False if save_sampled is None else save_sampled
        self.start = None if start is None else start
        self.method = method
        self.best = None
        self.speedMax = None
        self.speedMin = None
        self.size = None
        self.lb = None
        self.ub = None
        self.bounds_set = False
        self.range = 2
        self.population = []
        self.w = 1
        self.update_w = True
        self.update_scheme = 'constriction'
        if self.update_scheme == 'constriction':
            fi = 2.05 + 2.05
            self.w = 2.0/np.abs(2.0 - fi - np.sqrt(np.power(fi,2) - 4 * fi))
        
    def generate(self):
        start_position = np.random.uniform(self.lb,self.ub,self.size)
        part = creator.Particle(start_position)
        part.speed = np.random.uniform(self.speedMin,self.speedMax,self.size)
        part.smin = self.speedMin
        part.smax = self.speedMax
        return part
        
    def set_w(self,option):
        self.update_w = option
    def get_best_value(self):
        return self.best.fitness.values
    
    def get_history(self):
        return self.best_history
    
    def updateParticle(self,part, phi1, phi2):
        v_u1 = np.random.uniform(0,1,self.size)* phi1 * (part.best - part)
        v_u2 = np.random.uniform(0,1,self.size)* phi2 * (self.best - part)
        part.speed = self.w* (part.speed + v_u1 + v_u2)
        np.place(part.speed, part.speed < part.smin, part.smin)
        np.place(part.speed, part.speed > part.smax, part.smax)
        part += part.speed
        for i, pos in enumerate(part):
            if pos < self.lb[i]:
                part[i] = self.lb[i]
            elif pos > self.ub[i]:
                part[i] =  self.ub[i]
                
    def set_cost_function(self,cost_function):
        self.cost_function = cost_function
              
    def setup_pso(self):
        if self.speedMax == None or self.speedMin == None:
            self.set_speed()

        if self.bounds_set == False:
            self.set_bounds()
            
        if self.size == None:
            try:
                self.size = len(self.best)
            except:
                print "Error: Must provide a starting position in order to set size of each particle"
                print "**** Provide PSO.set_start_position() your initial starting coordinates ****"
                print "Exiting due to failure"
                quit()
        if self.cost_function == None:
            print "Error: Must set a cost function. Use PSO.set_cost_function()."
            print "Exiting due to failure"
            quit()
        self.toolbox = base.Toolbox()
        if self.method == 'single_min':
            creator.create("FitnessMin", base.Fitness,weights=(-1.00,))
        creator.create("Particle", np.ndarray, fitness=creator.FitnessMin, \
            speed=list,smin=list, smax=list, best=None)
        self.toolbox.register("particle", self.generate)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.particle)
        self.toolbox.register("update", self.updateParticle, phi1=2.05, phi2=2.05)
        self.toolbox.register("evaluate", self.cost_function)
        self.stats = tools.Statistics(lambda ind: ind.fitness.values)
        self.stats.register("avg", np.mean, axis=0)
        self.stats.register("std", np.std, axis=0)
        self.stats.register("min", np.min, axis=0)
        self.stats.register("max", np.max, axis=0)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        pool = multiprocessing.Pool(1)
        self.toolbox.register("map", pool.map)
        self.toolbox.register("close", pool.close)
    
    def update_connected(self):
        for part in self.population:
            if part.best is None or part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if self.best is None or self.best.fitness < part.fitness:
                self.best = creator.Particle(part)
                self.best.fitness.values = part.fitness.values
    
    def update_neighbor(self):
        
        for i in range(len(self.population)):
            self.population[i].best
            if self.population[i].best is None or self.population[i].best.fitness < self.population[i].fitness:
                self.population[i].best = creator.Particle(self.population[i])
                self.population[i].best.fitness.values = self.population[i].fitness.values
            if self.best is None or self.best.fitness < part.fitness:
                self.best = creator.Particle(part)
                self.best.fitness.values = part.fitness.values
    
    
    def return_ranked_populations(self):
        positions = np.zeros(np.shape(self.population))
        fitnesses = np.zeros(len(self.population))
        for n,part in enumerate(self.population):
            fitnesses[n] = part.best.fitness.values[0]
            positions[n] = part.best
        idx = np.argsort(fitnesses)
        return fitnesses[idx],positions[idx]
    
    def set_start_position(self,position):
        self.start = position
        self.size = len(position)
        
    def set_speed(self,speedMin=-10000,speedMax=10000):
        self.speedMin = speedMin
        self.speedMax = speedMax
        
    def set_bounds(self,range=None,lower = None, upper = None):
        if range == None:
            range = self.range
        if len(self.start)==1:
            if self.start == None:
                print "Error: Please set starting position before setting bounds"
                quit()
        if lower==None:
            lower = self.start - range
        if upper == None:
            upper = self.start + range
        self.lb = lower
        self.ub = upper
        self.bounds_set = True
        
    def run(self,num_particles,num_iterations):
        self.setup_pso()
        try:
            print "Initial cost function value " , self.cost_function(self.start)
        except:
            print "Error: Unable to calculate cost function value of initial condition."
            print "Make sure that the cost function is callable and you have set_start_position()."
            if self.solver == None:
                print "Solver is not set. Is this the source of error?"
                print "Try setting solver that is used in cost_function."
                print "**** PSO.set_solver(solver goes here) ****"
            print "Exiting due to failure"
            quit()
        history = np.zeros((num_iterations,len(self.start)))
        if self.save_sampled == True:
            self.all_history = np.zeros((num_iterations,num_particles, len(self.start)))
            self.all_fitness = np.zeros((num_iterations, num_particles))
        values = np.zeros(num_iterations)
        self.population = self.toolbox.population(num_particles)
        for g in range(1,num_iterations+1):
            if self.update_w == True:
                self.w = (num_iterations -g + 1. )/ num_iterations
            fitnesses = self.toolbox.map(self.toolbox.evaluate, self.population)
            for ind, fit in zip(self.population, fitnesses):
                ind.fitness.values = fit
            self.update_connected()
            #self.update_neighbor()
            for part in self.population:
                self.toolbox.update(part)
            values[g-1]=self.best.fitness.values[0]
            history[g-1] = self.best
            if self.save_sampled == True:
                curr_fit, curr_pop = self.return_ranked_populations()
                self.all_history[g-1,:,:] = curr_pop
                self.all_fitness[g-1,:] = curr_fit
            self.logbook.record(gen=g, evals=len(self.population), **self.stats.compile(self.population))
            if self.logbook.select('std')[-1] < 1e-6:
                break
            print(self.logbook.stream),self.best.fitness.values
        self.toolbox.close()
        return values,history

        
        
