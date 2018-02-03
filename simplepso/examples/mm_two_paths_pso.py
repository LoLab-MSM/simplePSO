import numpy as np
import os
from mm_two_paths_model import model
from pysb.simulator.scipyode import ScipyOdeSimulator
from simplepso.pso import PSO
import matplotlib.pyplot as plt


directory = os.path.dirname(__file__)
data_path = os.path.join(directory, 'data', 'product_data.npy')
product_data = np.load(data_path)
t = np.linspace(0, 50, 51)

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
k_ids = [p.value for p in model.parameters_rules()]
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 2
solver = ScipyOdeSimulator(model, tspan=t)


def display(position):
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    sim = solver.run(param_values=param_values)
    plt.plot(t, product_data, color='r', marker='.', linestyle=':')
    plt.plot(t, sim.observables['Product'], color='b')

def likelihood(position):
    Y = np.copy(position)
    param_values[rate_mask] = 10 ** Y
    sim = solver.run(param_values=param_values)
    e = np.sum((sim.observables['Product'] - product_data) ** 2)
    return e,


def run_example():
    best_pars = np.zeros((10000, len(model.parameters)))

    counter = 0
    for i in range(10000):
         pso = PSO(save_sampled=False, verbose=False, num_proc=25)
         pso.set_cost_function(likelihood)
         pso.set_start_position(xnominal)
         pso.set_bounds(2)
         pso.set_speed(-.25, .25)
         pso.run(25, 100)
         Y=np.copy(pso.best)
         param_values[rate_mask] = 10 ** Y
         if pso.values.min() < 0.03:
            best_pars[counter] = param_values
            counter += 1
         print (i, counter)


    np.save('/home/oscar/main/tropical_earm_local_signatures/new_pars', best_pars)

if __name__ == '__main__':
    run_example()