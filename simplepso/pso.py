# -*- coding: utf-8 -*-
import os
from concurrent.futures import ProcessPoolExecutor, Executor, Future
from copy import deepcopy

import numpy as np

from simplepso.logging import setup_logger

# set OMP_NUM_THREADS to 1 to prevent multi-processing to expand no each core
os.environ['OMP_NUM_THREADS'] = "1"


class Particle(object):
    """
    Particle to be used in the Swarm
    """

    def __init__(self, pos, fitness=None, smin=None, smax=None):
        self.pos = pos
        self.fitness = fitness
        self.smin = smin
        self.smax = smax
        self.__best = None

    @property
    def best(self):
        return self.__best

    @best.setter
    def best(self, new_best):
        self.__best = Particle(
            deepcopy(new_best.pos), deepcopy(new_best.fitness),
            deepcopy(new_best.smin), deepcopy(new_best.smax)
        )

# Used for printing below
header = '{:<10}' + "\t".join(['{:>12}'] * 5)
stats_header = header.format('iteration', 'best', 'mean', 'min', 'max', 'std')
stats_output = '{:<10}' + '\t'.join(['{:>12.3f}'] * 5)


class PSO(object):
    """
    Simple interface to run particle swarm optimization

    It can optimize parameters for a function using two run methods.
        run() : cost function gets a parameter vector and returns a scalar. Can
            be used with any callable function.
        run_ssa() : cost function gets multiple trajectories of a PySB model
            and returns a scalar. These trajectories are provided as
            a pandas.DataFrame. PySB dependent function.
    """

    def __init__(self, start=None, save_sampled=False, verbose=False,
                 shrink_steps=True):
        """

        Parameters
        ----------
        start : list
            Starting position
        save_sampled : bool
            Save each particles position and fitness over all iterations.
        verbose : bool
        shrink_steps : bool
            Shrink max distance traveled by each particle as the number of
            iterations increases
        """
        self.save_sampled = save_sampled
        if start is not None:
            self.set_start_position(start)
        else:
            self.start = None
            """Starting position"""
            self.size = None
        self.verbose = verbose
        self.best = None
        """Best Particle of population"""
        self.max_speed = None
        self.min_speed = None
        self.lb = None
        self.ub = None
        self.bounds_set = False
        self.range = 2
        self.all_history = None
        """History of all particles positions over all iterations. Only saved
        if save_sampled=True"""
        self.all_fitness = None
        """History of all particles finesses over all iterations. Only saved
        if save_sampled=True"""
        self.population = []
        """Population of particles"""
        self.population_fitness = []
        """Fitness values of population of particles"""
        self.values = []
        """Fitness values of the best particle for each iteration"""
        self.history = []
        """Best particle for each iteration"""
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

    def _setup_pso(self):
        """ Sets up everything for PSO. Only does once """
        if self.max_speed is None or self.min_speed is None:
            self.set_speed()

        if not self.bounds_set:
            self.set_bounds()

        if self.start is None:
            raise Exception("Must provide a starting position\n"
                            "Use PSO.set_start_position() \n")

        self._is_setup = True

    def _update_connected(self):
        """ Update the population of particles"""
        for part in self.population:
            if part.best is None or part.best.fitness > part.fitness:
                part.best = part
                # part.best = deepcopy(part)
                # part.best.fitness = deepcopy(part.fitness)
            if self.best is None or self.best.fitness > part.fitness:
                self.best = Particle(
                    deepcopy(part.pos),
                    deepcopy(part.fitness),
                    deepcopy(part.smin),
                    deepcopy(part.smax)
                )
                # self.best = deepcopy(part)
                # self.best.fitness = deepcopy(part.fitness)

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
        positions = []
        finesses = []
        for n, part in enumerate(self.population):
            finesses.append(part.best.fitness)
            positions.append(part.best.pos)

        positions = np.array(positions)
        finesses = np.array(finesses)
        idx = np.argsort(finesses)
        return finesses[idx], positions[idx]

    def set_start_position(self, position):
        """
        Set the starting position for the population of particles.

        Parameters
        ----------
        position : array

        """
        self.start = np.array(position)
        self.size = len(position)

    def set_speed(self, speed_min=-10000, speed_max=10000):
        """ Sets the max and min speed of the particles.

        This is usually a fraction of the bounds set with `set_bounds`. So if
        one sets the bound to be +/- 1 order of magnitude, you can set the
        speed to be -.1 and .1, allow the particles to update in 1/10 the
        parameter space. This keeps particles near their local position rather
        than jumping across parameter space quickly.

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
            If provided parameter_range, the range will be set by the starting
            position +/- this value. To set each parameter manually, use
            `lower` and `upper` args
        lower : array
            Lower bounds for parameters
        upper : array
            Upper bounds for parameters

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

    def run(self, num_particles, num_iterations, cost_function=None,
            num_processors=1, save_samples=False, stop_threshold=1e-5,
            max_iter_no_improv=None):
        """
        Run optimization

        Parameters
        ----------
        num_particles : int
            Number of particles in population, ~20 is a good starting place
        num_iterations : int
            Number of iterations to perform.
        cost_function : callable function
            Takes a parameter set and returns a scalar (particles fitness)
        num_processors : int
            Number of processors to run on. If using scipy, note that you may
            need to set OMP_NUM_THREADS=1 to prevent each process from using
            more than one CPU.
        save_samples : bool
            Save ALL positions of particles over time, can require large memory
            if num_particles, num_iterations, and len(parameters) is large.
        stop_threshold : float
             Standard deviation of the particles’ cost function at which the
             optimization is stopped
        max_iter_no_improv: int
            Maximum steps allowed without improvement before the optimization
            stops.
        """
        if self._is_setup:
            pass
        else:
            self._setup_pso()

        if not callable(cost_function):
            raise TypeError("Provide a callable cost function")

        history = np.zeros((num_iterations, len(self.start)))
        if self.save_sampled or save_samples:
            self.all_history = np.zeros(
                (num_iterations, num_particles, len(self.start)))
            self.all_fitness = np.zeros((num_iterations, num_particles))
        values = np.zeros(num_iterations)
        self.population = [self._generate() for _ in range(num_particles)]
        self.population_fitness = np.zeros(len(self.population))
        if max_iter_no_improv is None:
            max_iter_no_improv = np.inf
        iter_without_improvement = 0
        best_fitness = np.inf
        with SerialExecutor() if num_processors == 1 else \
                ProcessPoolExecutor(max_workers=num_processors) as executor:
            for g in range(num_iterations):
                if self.update_w:
                    self.w = (num_iterations - g + 1.) / num_iterations
                self.population_fitness = np.array(
                    list(executor.map(
                        cost_function,
                        [p.pos for p in self.population],
                    )
                    )
                )

                for ind, fit in zip(self.population, self.population_fitness):
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
                    self.print_stats(g + 1, fitness=self.population_fitness)

                if self.population_fitness.std() < stop_threshold:
                    self.log.info("Stopping. STD < stop_threshold.")
                    break
                if self.best.fitness < best_fitness:
                    best_fitness = self.best.fitness
                    iter_without_improvement = 0
                else:
                    iter_without_improvement += 1
                if iter_without_improvement > max_iter_no_improv:
                    self.log.info("Stopping. max_iter_no_improv reached.")
                    break
        self.values = np.array(values[:g])
        self.history = np.array(history[:g, :])

    def _calc_fitness_from_array(self, traj, num_sim, cost_function):
        index_names = traj.index.names
        traj.reset_index(inplace=True)
        for i, part in enumerate(self.population):
            start = i * num_sim
            end = (i + 1) * num_sim
            tmp_results = traj.loc[
                traj.simulation.isin(list(range(start, end)))]
            tmp_results.set_index(index_names, inplace=True)
            error = cost_function(tmp_results)
            part.fitness = error
            self.population_fitness[i] = error

    def _get_parameters_from_population(self, num_sims, total_sims, n_params):
        """ Create param_array for GPU based simulators """
        rate_vals = np.zeros((total_sims, n_params))
        # create parameters for each particle, creates blocks per num_sims
        for i, part in enumerate(self.population):
            start = i * num_sims
            end = (i + 1) * num_sims
            rate_vals[start:end, :] = part.pos
        return rate_vals

    def run_ssa(self, model, num_particles, num_iterations, num_sim,
                cost_function, simulator, save_samples=False,
                stop_threshold=0):
        """
        Run PSO for a stochastic simulator

        Parameters
        ----------
        model : pysb.Model
        num_particles : int
            Number of particles in the swarm.
        num_iterations : int
            Number of iterations to perform
        num_sim : int
            Number of SSA simulations to run for each particle.
        cost_function : function
            Takes a pandas dataframe of PySB trajectories created by running
            multiple SSA simulations. Function must return a scalar.
        simulator : pysb.Simulator
            PySB simulator (CudaSSASimulator or OpenclSSASimulator)
        save_samples : bool
        stop_threshold : float
        """
        if self._is_setup:
            pass
        else:
            self._setup_pso()
        if simulator is None:
            raise ValueError(
                "Must provide an SSA simulator to use this method")

        history = []
        values = []
        self.population = [self._generate() for _ in range(num_particles)]
        self.population_fitness = np.zeros((len(self.population)))
        total_sims = num_particles * num_sim
        rate_p = model.parameters_rules()
        rate_mask = np.array([p in rate_p for p in model.parameters])
        nominal_values = np.array([p.value for p in model.parameters])
        log10_original_values = np.log10(nominal_values[rate_mask])
        all_param_vals = np.repeat([nominal_values], total_sims, axis=0)

        for g in range(0, num_iterations):
            if self.update_w:
                self.w = (num_iterations - g) / num_iterations
            rate_values = self._get_parameters_from_population(
                num_sim, total_sims, len(log10_original_values)
            )

            # convert back from log10 space
            all_param_vals[:, rate_mask] = 10 ** rate_values
            # reset initials and param_values so simulator won't try to
            # duplicate any
            simulator.initials = None
            simulator.param_values = None
            traj = simulator.run(param_values=all_param_vals).dataframe
            self._calc_fitness_from_array(traj, num_sim, cost_function)
            self._update_connected()
            for part in self.population:
                self._update_particle_position(part)

            values.append(self.best.fitness)
            history.append(deepcopy(self.best.pos))

            if self.save_sampled or save_samples:
                curr_fit, curr_pop = self.return_ranked_populations()
                self.all_history[g, :, :] = curr_pop
                self.all_fitness[g, :] = curr_fit

            if self.verbose:
                self.print_stats(g + 1, fitness=self.population_fitness)
            if self.population_fitness.std() < stop_threshold:
                self.log.info("Stopping criteria reached.")
                break

        self.values = np.array(values)
        self.history = np.array(history)

    def print_stats(self, iteration, fitness):
        if iteration == 1:
            self.log.info(stats_header)
        self.log.info(
            stats_output.format(iteration, self.best.fitness, fitness.mean(),
                                fitness.min(), fitness.max(), fitness.std())
        )


class SerialExecutor(Executor):
    """ Execute tasks in serial (immediately on submission)
    Code originally from PySB.
    """

    def submit(self, fn, *args, **kwargs):
        f = Future()
        try:
            result = fn(*args, **kwargs)
        except BaseException as e:
            f.set_exception(e)
        else:
            f.set_result(result)

        return f
