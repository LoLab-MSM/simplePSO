# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""

import operator
import random
import scipy.optimize
import numpy
import pylab
from deap import base
from deap import creator
from deap import tools
import pysb.integrate
import pysb.util
import numpy as np
import scipy.optimize
import scipy.interpolate
import matplotlib.pyplot as plt
import os
import inspect
import multiprocessing
import multiprocessing as mp
import sys
from pysb.util import load_params

class PSO():
    def __init__(self):
        self.solver = None
        self.speedMax = None
        self.speedMin = None
        self.cost_function = None
        self.best = None
        self.size = None
        self.lb = None
        self.ub = None
        self.best_history = None
        self.best_value_history = None
        self.range = 2
        self.population = []
    def generate(self,size):
        u1 = np.random.uniform(self.lb,self.ub,size)
        start_position = u1
        part = creator.Particle(start_position)
        part.speed = np.random.uniform(self.speedMin,self.speedMax,size)
        part.smin = self.speedMin
        part.smax = self.speedMax
        return part
    
    def get_history_value(self):
        return self.best_value_history
    
    def get_history(self):
        return self.best_history
    
    def updateParticle(self,part, best, phi1, phi2):
    
        u1 = numpy.random.uniform(0, 1, len(part))
        u2 = numpy.random.uniform(0, 1, len(part))
        v_u1 = u1 * phi1 * (part.best - part)
        v_u2 = u2 * phi2 * (best - part)
        part.speed += v_u1 + v_u2
        np.place(part.speed,part.speed<part.smin,part.smin)
        np.place(part.speed,part.speed>part.smax,part.smax)
        part += part.speed
        for i, pos in enumerate(part):
            if pos < self.lb[i]:
                part[i] = self.lb[i]
            elif pos > self.ub[i]:
                part[i] =  self.ub[i]
    def set_cost_function(self,cost_function):
        self.cost_function = cost_function
    
    def set_solver(self,solver_init):
        global solver
        solver = solver_init        
    def setup_pso(self):
        if self.speedMax == None or self.speedMin == None:
            self.set_speed()

        if self.lb == None or self.ub == None:
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
        toolbox = base.Toolbox()
        creator.create("FitnessMin", base.Fitness,weights=(-1.00,))
        creator.create("Particle", np.ndarray, fitness=creator.FitnessMin, \
            speed=list,smin=list, smax=list, best=self.best)
        
        toolbox.register("particle", self.generate, size=self.size)
        toolbox.register("population", tools.initRepeat, list, toolbox.particle)
        toolbox.register("update", self.updateParticle, phi1=2, phi2=2)
        toolbox.register("evaluate", self.cost_function)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", numpy.mean, axis=0)
        stats.register("std", numpy.std, axis=0)
        stats.register("min", numpy.min, axis=0)
        stats.register("max", numpy.max, axis=0)
        logbook = tools.Logbook()
        logbook.header = ["gen", "evals"] + stats.fields
        pool = multiprocessing.Pool(processes=4)
        toolbox.register("map", pool.map)
        globals().update(locals())
        return toolbox,stats,logbook,creator
    
    def set_start_position(self,position):
        self.best = position
        self.size = len(position)
        
    def set_speed(self,speedMin=-10000,speedMax=10000):
        self.speedMin = speedMin
        self.speedMax = speedMax
    
    def set_bounds(self,range=None):
        if range == None:
            range = self.range
        if len(self.best)==1:
            if self.best == None:
                print "Error: Please set starting position before setting bounds"
                quit()
        self.lb = self.best - range
        self.ub = self.best + range
        print self.lb,self.ub
        
    def run(self,num_particles,num_iterations):
        toolbox, stats, logbook,creator = self.setup_pso()
        try:
            print "Initial cost function value %s" % self.cost_function(self.best)
        except:
            print "Error: Unable to calculate cost function value of initial condition."
            print "Make sure that the cost function is callable and you have set_start_position()."
            if self.solver == None:
                print "Solver is not set. Is this the source of error?"
                print "Try setting solver that is used in cost_function."
                print "**** PSO.set_solver(solver goes here) ****"
            print "Exiting due to failure"
            quit()
        num_iterations = num_iterations
        num_particles = num_particles
        best = creator.Particle(self.best)
        best.fitness.values = self.cost_function(self.best)
        self.population = toolbox.population(num_particles)
        for g in range(1,num_iterations+1):
            fitnesses = toolbox.map(toolbox.evaluate, self.population)
            for ind, fit in zip(self.population, fitnesses):
                ind.fitness.values = fit
            for part in self.population:
                if not g == 1:
                    if part.fitness.values < part.best.fitness.values:
                        part.best = creator.Particle(part)
                        part.best.fitness.values = part.fitness.values
                else:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
    
                if part.fitness.values < best.fitness.values:
                    best = creator.Particle(part)
                    best.fitness.values = part.fitness.values
            for part in self.population:
                toolbox.update(part, best)
            logbook.record(gen=g, evals=len(self.population), **stats.compile(self.population))
            print(logbook.stream),best.fitness.values
        if self.best_value_history == None or self.best_value_history> best.fitness.values:
            self.best_value_history = best.fitness.values
            self.best_history = best
        
        