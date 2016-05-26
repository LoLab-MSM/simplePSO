# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""

import multiprocessing

import numpy as np
from deap import base
from deap import creator
from deap import tools


class PSO:
    def __init__(self, cost_function=None, start=None, save_sampled=False, verbose=False):
        """

        :param cost_function:
        :param start:
        :param save_sampled:
        """

        self.cost_function = cost_function
        self.save_sampled = save_sampled
        self.start = start
        self.size = None
        if self.start is not None:
            self.set_start_position(start)
        self.verbose = verbose
        self.method = 'single_min'
        self.best = None
        self.max_speed = None
        self.min_speed = None
        self.lb = None
        self.ub = None
        self.bounds_set = False
        self.range = 2
        self.population = []
        self.toolbox = base.Toolbox()
        self.stats = tools.Statistics(lambda ind: ind.fitness.values)
        self.logbook = tools.Logbook()
        self.all_history = None
        self.all_fitness = None
        self.values = []
        self.history = []
        self.w = 1
        self.update_w = True
        self._is_setup = False
        self.update_scheme = 'constriction'
        if self.update_scheme == 'constriction':
            fi = 2.05 + 2.05
            self.w = 2.0 / np.abs(2.0 - fi - np.sqrt(np.power(fi, 2) - 4 * fi))

    def generate(self):
        start_position = np.random.uniform(self.lb, self.ub, self.size)
        part = creator.Particle(start_position)
        part.speed = np.random.uniform(self.min_speed, self.max_speed, self.size)
        part.smin = self.min_speed
        part.smax = self.max_speed
        return part

    def set_w(self, option):
        self.update_w = option

    def get_best_value(self):
        return self.best.fitness.values

    def get_history(self):
        return self.best_history

    def update_particle_position(self, part, phi1, phi2):
        v_u1 = np.random.uniform(0, 1, self.size) * phi1 * (part.best - part)
        v_u2 = np.random.uniform(0, 1, self.size) * phi2 * (self.best - part)
        part.speed = self.w * (part.speed + v_u1 + v_u2)
        np.place(part.speed, part.speed < part.smin, part.smin)
        np.place(part.speed, part.speed > part.smax, part.smax)
        part += part.speed
        for i, pos in enumerate(part):
            if pos < self.lb[i]:
                part[i] = self.lb[i]
            elif pos > self.ub[i]:
                part[i] = self.ub[i]

    def set_cost_function(self, cost_function):
        self.cost_function = cost_function

    def setup_pso(self):
        if self.max_speed is None or self.min_speed is None:
            self.set_speed()

        if not self.bounds_set:
            self.set_bounds()

        assert self.start is not None, "Error: Must provide a starting position in order to set size of each particle\
                       **** Provide PSO.set_start_position() your initial starting coordinates ****\
                       Exiting due to failure"

        assert self.cost_function is not None, "Error: Must set a cost function. Use PSO.set_cost_function()."

        if self.method == 'single_min':
            creator.create("FitnessMin", base.Fitness, weights=(-1.00,))
        creator.create("Particle", np.ndarray, fitness=creator.FitnessMin, speed=list, smin=list, smax=list, best=None)
        self.toolbox.register("particle", self.generate)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.particle)
        self.toolbox.register("update", self.update_particle_position, phi1=2.05, phi2=2.05)
        self.toolbox.register("evaluate", self.cost_function)

        self.stats.register("avg", np.mean, axis=0)
        self.stats.register("std", np.std, axis=0)
        self.stats.register("min", np.min, axis=0)
        self.stats.register("max", np.max, axis=0)

        self.logbook.header = ["gen", "evals"] + self.stats.fields
        pool = multiprocessing.Pool(1)
        self.toolbox.register("map", pool.map)
        self.toolbox.register("close", pool.close)
        self._is_setup = True

    def update_connected(self):
        for part in self.population:
            if part.best is None or part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if self.best is None or self.best.fitness < part.fitness:
                self.best = creator.Particle(part)
                self.best.fitness.values = part.fitness.values

    def return_ranked_populations(self):
        positions = np.zeros(np.shape(self.population))
        fitnesses = np.zeros(len(self.population))
        for n, part in enumerate(self.population):
            fitnesses[n] = part.best.fitness.values[0]
            positions[n] = part.best
        idx = np.argsort(fitnesses)
        return fitnesses[idx], positions[idx]

    def set_start_position(self, position):
        self.start = position
        self.size = len(position)

    def set_speed(self, speed_min=-10000, speed_max=10000):
        self.min_speed = speed_min
        self.max_speed = speed_max

    def set_bounds(self, parameter_range=None, lower=None, upper=None):
        assert self.start is not None, "Must provide starting array: %r" % self.start
        all_set_to_none = False
        if parameter_range is None and upper is None and lower is None:
            all_set_to_none = True
        assert all_set_to_none == False,'Need to provide parameter range or upper and lower bounds'
        if parameter_range is None:
            assert self.range is not None
            parameter_range = self.range

        if lower is None:
            lower = self.start - parameter_range
        else:
            assert self.size == len(lower), "If providing array for bounds, must equal length of starting position"
        if upper is None:
            upper = self.start + parameter_range
        else:
            assert self.size == len(upper), "If providing array for bounds, must equal length of starting position"
        self.lb = lower
        self.ub = upper
        self.bounds_set = True

    def run(self, num_particles, num_iterations):
        if self._is_setup:
            pass
        else:
            self.setup_pso()
        assert type(self.cost_function(self.start)) == tuple, "Cost function must return a tuple"
        try:
            print "Initial cost function value ", self.cost_function(self.start)
        except:
            print "Error: Unable to calculate cost function value of initial condition."
            print "Make sure that the cost function is callable and you have set_start_position()."
            print "Exiting due to failure"
            quit()
        history = np.zeros((num_iterations, len(self.start)))
        if self.save_sampled:
            self.all_history = np.zeros((num_iterations, num_particles, len(self.start)))
            self.all_fitness = np.zeros((num_iterations, num_particles))
        values = np.zeros(num_iterations)
        self.population = self.toolbox.population(num_particles)
        for g in range(1, num_iterations + 1):
            if self.update_w:
                self.w = (num_iterations - g + 1.) / num_iterations
            fitnesses = self.toolbox.map(self.toolbox.evaluate, self.population)
            for ind, fit in zip(self.population, fitnesses):
                ind.fitness.values = fit
            self.update_connected()
            for part in self.population:
                self.toolbox.update(part)
            values[g - 1] = self.best.fitness.values[0]
            history[g - 1] = self.best
            if self.save_sampled:
                curr_fit, curr_pop = self.return_ranked_populations()
                self.all_history[g - 1, :, :] = curr_pop
                self.all_fitness[g - 1, :] = curr_fit
            self.logbook.record(gen=g, evals=len(self.population), **self.stats.compile(self.population))
            if self.logbook.select('std')[-1] < 1e-6:
                break
            if self.verbose:
                print(self.logbook.stream), self.best.fitness.values
        self.toolbox.close()
        self.values = values
        self.history = history

