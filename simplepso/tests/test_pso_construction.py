import numpy as np
from nose.tools import raises

from simplepso.pso import PSO


def h1(individual):
    """ Simple two-dimensional function containing several local maxima.
    Found in deap.benchmarks.h1
    Defined range of [-100, 100]
    minimum is at f(8.6998, 6.7665) = 2
    """
    num = (np.sin(individual[0] - individual[1] / 8)) ** 2 + (np.sin(individual[1] + individual[0] / 8)) ** 2
    denum = ((individual[0] - 8.6998) ** 2 + (individual[1] - 6.7665) ** 2) ** 0.5 + 1
    return -1 * num / denum,


def himmelblau(individual):
    """The Himmelblau's function is multimodal with 4 defined minimums in
    :math:`[-6, 6]^2.

        range [-6, 6]
        x_1 = (3.0, 2.0), = 0
        x_2 = (-2.805118, 3.131312), = 0
        x_3 = (-3.779310, -3.283186), = 0
        x_4 = (3.584428, -1.848126), = 0
    """
    return (individual[0] * individual[0] + individual[1] - 11) ** 2 + \
           (individual[0] + individual[1] * individual[1] - 7) ** 2,


def test_population_creation():
    pso = PSO(cost_function=h1, start=[10, 0], verbose=False)
    pso.set_bounds(lower=[-100, -100], upper=[100, 100])
    pso.run(num_iterations=100, num_particles=10)
    pso.return_ranked_populations()
    error = np.sum((pso.best - [8.6998, 6.7665]) ** 2)
    print('True value: [8.6998, 6.7665]. Found:{0}. Error^2 = {1}'.format(pso.best, error))
    assert (error < 0.1)


def test_himmelblau():
    """ test to see if PSO can find simple minimum
    """
    minimums = [[3.0, 2.0],
                [-2.805118, 3.131312],
                [-3.779310, -3.283186],
                [3.584428, -1.848126]]

    pso = PSO(cost_function=himmelblau, start=[10, 0], verbose=False)
    pso.set_bounds(lower=[-100, -100], upper=[100, 100])
    pso.run(num_iterations=100, num_particles=10)
    good_min = False
    for i in minimums:
        if np.sum((pso.best - i) ** 2) < .1:
            good_min = True
            error = np.sum((pso.best - i) ** 2)
            found_min = i
    if good_min:
        print('Found minimum')
        print('True value: {0}. Found:{1}. Error^2 = {2}'.format(found_min, pso.best, error))


@raises(AssertionError)
def test_missing_cost_function():
    pso = PSO(start=[10, 0], verbose=False)
    pso.set_bounds(lower=[-100, -100], upper=[100, 100])
    pso.run(num_iterations=100, num_particles=10)


@raises(AssertionError)
def test_mismatched_bounds():
    pso = PSO(start=[10, 0], cost_function=himmelblau, verbose=False)
    pso.set_bounds(lower=[-100, 0, -100], upper=[100, 100])
    pso.run(num_iterations=100, num_particles=10)

@raises(AssertionError)
def test_no_bounds():
    pso = PSO(start=[10, 0], cost_function=himmelblau, verbose=False)
    pso.run(num_iterations=100, num_particles=10)
