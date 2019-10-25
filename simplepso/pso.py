# -*- coding: utf-8 -*-
import concurrent.futures
import os
from copy import deepcopy

import numpy as np

from simplepso.logging import setup_logger

os.environ['OMP_NUM_THREADS'] = "1"


class Particle(object):
    def __init__(self, pos):
        self.pos = pos
        self.fitness = None
        self.best = None
        self.smin = None
        self.smax = None


# Used for printing below
header = '{:<10}' + "\t".join(['{:>12}'] * 5)
stats_header = header.format('iteration', 'best', 'mean', 'min', 'max', 'std')
stats_output = '{:<10}' + '\t'.join(['{:>12.3f}'] * 5)


class PSO(object):
    """ Simple interface to run particle swarm optimization

    This class provides a simple interface to run particle swarm optimization.
    It builds off the deap package, but provides a simple interface.

    Examples:

        optimizer = simplepso.pso.PSO(cost_function,start_position)
        optimizer.set_bounds()
        optimizer.set_speed()
        optimizer.run(number_of_particles, number_of_iterations)


    Notes:
        Must set 1.) cost_function, 2.) starting position, 3.) bounds
        To run need to supply number of particles and number of iterations.
        20 particles is a good starting place.

    """

    def __init__(self, cost_function=None, start=None, num_proc=1,
                 save_sampled=False, verbose=False, shrink_steps=True,
                 simulator=None):
        """

        Parameters
        ----------
        cost_function : function
            Callable function that takes a parameter and returns a tuple
        start : list_list
            Starting position
        num_proc : int
            Number of processors to run on. If using scipy, note that you may
            need to set OMP_NUM_THREADS=1 to prevent each process from using
            more than one CPU.
        save_sampled : bool
            Save each position samples
        verbose : bool

        """
        self.cost_function = cost_function
        self.save_sampled = save_sampled
        if start is not None:
            self.set_start_position(start)
        else:
            self.start = None
            self.size = None
        if simulator is not None:
            self.simulator = simulator
        self.verbose = verbose
        self.num_proc = num_proc
        self.best = None
        self.max_speed = None
        self.min_speed = None
        self.lb = None
        self.ub = None
        self.bounds_set = False
        self.range = 2
        self.all_history = None
        self.all_fitness = None
        self.population = []
        self.values = []
        self.history = []
        self.update_w = shrink_steps
        self._is_setup = False
        self.log = setup_logger(verbose)
        fi = 2.05 + 2.05  # value from kennedy paper
        self.w = 2.0 / np.abs(2.0 - fi - np.sqrt(np.power(fi, 2) - 4 * fi))

    def get_best_value(self):
        """ Returns the best fitness value of the population """
        return self.best.fitness

    def get_history(self):
        """ Returns the history of the run"""
        return self.history

    def _generate(self):
        """ Creates Particles and sets their speed """
        part = Particle(np.random.uniform(self.lb, self.ub, self.size))
        part.speed = np.random.uniform(self.min_speed, self.max_speed,
                                       self.size)
        part.smin = self.min_speed
        part.smax = self.max_speed
        return part

    def set_cost_function(self, cost_function):
        """ Cost function must return a scalar followed by a , (tuple)

        Parameters
        ----------
        cost_function : callable function

        """
        self.cost_function = cost_function

    def _setup_pso(self):
        """ Sets up everything for PSO. Only does once """
        if self.max_speed is None or self.min_speed is None:
            self.set_speed()

        if not self.bounds_set:
            self.set_bounds()

        if self.start is None:
            raise Exception("Must provide a starting position\n"
                            "Use PSO.set_start_position() \n")

        if self.cost_function is None:
            raise Exception("No cost function. "
                            "Set with PSO.set_cost_function().")

        self._is_setup = True

    def _update_connected(self):
        """ Update the population of particles"""
        for part in self.population:
            if part.best is None or part.best.fitness > part.fitness:
                part.best = deepcopy(part)
                part.best.fitness = deepcopy(part.fitness)
            if self.best is None or self.best.fitness > part.fitness:
                self.best = deepcopy(part)
                self.best.fitness = deepcopy(part.fitness)

    def _update_particle_position(self, part):
        """ Updates an individual particles position """
        phi1 = 2.05
        phi2 = 2.05
        v_u1 = np.random.uniform(0, 1, self.size) * phi1 * (
                part.best.pos - part.pos)
        v_u2 = np.random.uniform(0, 1, self.size) * phi2 * (
                self.best.pos - part.pos)
        part.speed = self.w * (part.speed + v_u1 + v_u2)
        np.place(part.speed, part.speed < part.smin, part.smin)
        np.place(part.speed, part.speed > part.smax, part.smax)
        part.pos += part.speed
        for i, pos in enumerate(part.pos):
            if pos < self.lb[i]:
                part.pos[i] = self.lb[i]
            elif pos > self.ub[i]:
                part.pos[i] = self.ub[i]

    def return_ranked_populations(self):
        """ Returns a ranked population of the particles

        Returns
        -------
        np.array, np.array
        """
        positions = np.zeros(np.shape(self.population))
        fitnesses = np.zeros(len(self.population))
        for n, part in enumerate(self.population):
            fitnesses[n] = part.best.fitness
            positions[n] = part.best.pos
        idx = np.argsort(fitnesses)
        return fitnesses[idx], positions[idx]

    def set_start_position(self, position):
        """ Set the starting position for the population of particles.

        Parameters
        ----------
        position : array

        """
        self.start = np.array(position)
        self.size = len(position)

    def set_speed(self, speed_min=-10000, speed_max=10000):
        """ Sets the max and min speed of the particles.

        This is usually a fraction of the range of values that the particles
        can travel.

        Parameters
        ----------
        speed_min : float
        speed_max : float

        """
        self.min_speed = speed_min
        self.max_speed = speed_max

    def set_bounds(self, parameter_range=None, lower=None, upper=None):
        """ Set the search space bounds that the parameters may search.

        This can be a single float, in which the particles can search plus or
        minus the starting values around this. It can also be an array that
        must be the same length of the starting position.

        Parameters
        ----------
        parameter_range : float
        lower : array of lower bounds for parameters
        upper : array of upper bounds for parameters

        """
        if self.start is None:
            raise Exception("Must provide starting array: %r" % self.start)

        if parameter_range is None and upper is None and lower is None:
            raise Exception('Need to provide parameter ranges')

        if parameter_range is None:
            if self.range is None:
                raise Exception("Please provide parameter range")
            parameter_range = self.range

        if lower is None:
            lower = self.start - parameter_range
        else:
            if not self.size == len(lower):
                raise Exception("Size of bounds must equal size of particle")
        if upper is None:
            upper = self.start + parameter_range
        else:
            if not self.size == len(upper):
                raise Exception("Size of bounds must equal size of particle")
        self.lb = lower
        self.ub = upper
        self.bounds_set = True

    def calc_fitness(self, traj, num_sim):
        population_fitness = np.zeros((len(self.population)))
        for i, part in enumerate(self.population):
            start = i * num_sim
            end = (i + 1) * num_sim
            tmp_results = traj[list(range(start, end))]
            tmp_results = tmp_results.T.unstack('simulation')
            error = self.cost_function(tmp_results)
            part.fitness = error
            population_fitness[i] = error
        return population_fitness

    def get_parameters_from_population(self, num_sims, total_sims, n_params):
        rate_vals = np.zeros((total_sims, n_params))
        # create parameters for each particle
        for i, part in enumerate(self.population):
            start = i * num_sims
            end = (i + 1) * num_sims
            rate_vals[start:end, :] = part.pos
        return rate_vals

    def run(self, num_particles, num_iterations, num_processes=1,
            save_samples=False, stop_threshold=1e-5):
        """ Run optimization

        Parameters
        ----------
        num_particles : int
            Number of particles in population, ~20 is a good starting place
        num_iterations : int
        num_processes : int
            Number of processes (cpu cores) to use
        save_samples : bool
            Save positions of particles over time, can require large memory
            if num_particles, num_iterations, and len(parameters) is large.
        stop_threshold : float
            Threshold of standard deviation of all particles cost function.

        """
        if self._is_setup:
            pass
        else:
            self._setup_pso()

        history = np.zeros((num_iterations, len(self.start)))
        if self.save_sampled or save_samples:
            self.all_history = np.zeros(
                (num_iterations, num_particles, len(self.start)))
            self.all_fitness = np.zeros((num_iterations, num_particles))
        values = np.zeros(num_iterations)
        self.population = [self._generate() for _ in range(num_particles)]

        with concurrent.futures.ProcessPoolExecutor(num_processes) as executor:
            for g in range(num_iterations):
                if self.update_w:
                    self.w = (num_iterations - g + 1.) / num_iterations
                fitness = list(executor.map(
                    self.cost_function, [p.pos for p in self.population],
                    chunksize=int(num_particles / num_processes))
                )
                population_fitness = np.array(fitness)
                for ind, fit in zip(self.population, population_fitness):
                    ind.fitness = fit
                self._update_connected()
                for part in self.population:
                    self._update_particle_position(part)
                values[g] = self.best.fitness
                history[g] = self.best.pos
                if self.save_sampled or save_samples:
                    curr_fit, curr_pop = self.return_ranked_populations()
                    self.all_history[g, :, :] = curr_pop
                    self.all_fitness[g, :] = curr_fit

                if self.verbose:
                    self.print_stats(g + 1, fitness=population_fitness)

                if population_fitness.std() < stop_threshold:
                    self.log.info("Stopping criteria reached.")
                    break
        self.values = values[:g]
        self.history = history[:g, :]

    def run_ssa(self, model, num_particles, num_iterations, num_sim,
                save_samples=False, stop_threshold=0):
        """ Run PSO for a stochastic simulator

        Parameters
        ----------
        model : pysb.Model
        num_particles : int
        num_iterations : int
        num_sim : int
        save_samples : bool
        stop_threshold : float

        """
        if self._is_setup:
            pass
        else:
            self._setup_pso()
        if self.simulator is None:
            raise Exception("Must provide an SSA simulator to use this method")

        history = np.zeros((num_iterations, len(self.start)))
        values = np.zeros(num_iterations)
        # self.population = self.toolbox.population(num_particles)
        self.population = [self._generate() for _ in range(num_particles)]
        total_sims = num_particles * num_sim
        rate_p = model.parameters_rules()
        rate_mask = np.array([p in rate_p for p in model.parameters])
        nominal_values = np.array([p.value for p in model.parameters])
        log10_original_values = np.log10(nominal_values[rate_mask])
        all_param_vals = np.repeat([nominal_values], total_sims, axis=0)

        for g in range(0, num_iterations):
            if self.update_w:
                self.w = (num_iterations - g) / num_iterations
            rate_values = self.get_parameters_from_population(
                num_sim, total_sims, len(log10_original_values)
            )

            # convert back from log10 space
            all_param_vals[:, rate_mask] = 10 ** rate_values
            # reset initials and param_values so simulator won't try to
            # duplicate any
            self.simulator.initials = None
            self.simulator.param_values = None
            traj = self.simulator.run(param_values=all_param_vals).dataframe.T
            population_fitness = self.calc_fitness(traj, num_sim)
            self._update_connected()
            for part in self.population:
                self._update_particle_position(part)

            values[g] = self.best.fitness
            history[g] = self.best.pos

            if self.save_sampled or save_samples:
                curr_fit, curr_pop = self.return_ranked_populations()
                self.all_history[g, :, :] = curr_pop
                self.all_fitness[g, :] = curr_fit

            if self.verbose:
                self.print_stats(g + 1, fitness=population_fitness)
            if population_fitness.std() < stop_threshold:
                self.log.info("Stopping criteria reached.")
                break
        self.values = values[:g]
        self.history = history[:g, :]

    def print_stats(self, iteration, fitness):
        if iteration == 1:
            self.log.info(stats_header)
        self.log.info(
            stats_output.format(iteration, self.best.fitness, fitness.mean(),
                                fitness.min(), fitness.max(), fitness.std())
        )
