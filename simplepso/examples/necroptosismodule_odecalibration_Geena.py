import datetime
now = datetime.datetime.now()

try:
    import matplotlib

    matplotlib.use('TkAgg')
    matplotlib.use('Agg') #this is needed only if you run it in a cluster
    import matplotlib.pyplot as plt
    plot = True #this is because in the run_example function you have the "if plot" expression
except ImportError:
    plot = False
    pass

import numpy as np
from necroptosismodule import model # I don't know why this is reporterd as error, because this works just fine with importing
from pysb.simulator import ScipyOdeSimulator

from simplepso.pso import PSO

model.enable_synth_deg()
"""Add components needed to support synthesis and degradation rules."""
obs_names = ['MLKLa_obs']

# Defining a few helper functions to use
def normalize(trajectories):
    """even though this is not really needed, if the data is already between 1 and 0!"""
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)

# t = np.array([0, 60, 120, 240, 360, 480]) # tspan of the data, in minutes
t = np.array([0., 30,  60,   120,  180, 270,  480,  960, 1440])
solver = ScipyOdeSimulator(model, tspan=t) #, rtol=1e-6, # rtol : float or sequence relative tolerance for solution
                            #atol=1e-6) #atol : float or sequence absolute tolerance for solution

rate_params = model.parameters_rules() # these are only the parameters involved in the rules
param_values = np.array([p.value for p in model.parameters]) # these are all the parameters
rate_mask = np.array([p in rate_params for p in model.parameters])  # this picks the element of intersection
# original_values = np.array([p.value for p in model.parameters]) #also not needed here

# Data
# data = np.array([0., 0., 0., 0.10, 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.,1.])
# data = np.array([0, 0.017, 0.09, 0.25, 0.488, 0.994, 1])
# data = np.array([0., 0., 0., 0., 0.03, 0.10, 0.5, 0.99, 1.])

# t = np.array([0., 30,  60,   120,  180, 270,  480,  960, 1440])
#
data = np.array([0., 0., 0., 0., 0.01, 0.05, 0.5, 0.99, 1.])

# We search in log10 space for the parameters - this relates to the parameters at the beginning!
log10_original_values = np.log10(param_values[rate_mask]) # this with the rate mask is needed because the solver wants the parameter values to have the same dimention as the whole parameter values for the model, same at *; it is a way of writing that only takes the parameters for which rate_mask is True

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# It must return a tuple
def obj_function(params):
    params_tmp = np.copy(params) #here you pass the parameter vector; the point of making a copy of it is in order not to modify it
    param_values[rate_mask] = 10 ** params_tmp # see comment above *
    result = solver.run(param_values = param_values)
    ysim_norm = normalize(result.observables['MLKLa_obs'])
    error = np.sum((data - ysim_norm) ** 2)
    return error, #the comma is for returning a touple

def run_example():
    for i in range(10000):
    # Here we initialize the class
    # We must proivde the cost function and a starting value
        optimizer = PSO(cost_function=obj_function, start=log10_original_values, verbose=True)
        # We also must set bounds. This can be a single scalar or an array of len(start_position)
        optimizer.set_bounds(parameter_range=2)
        optimizer.set_speed(speed_min=-.25, speed_max=.25)
        optimizer.run(num_particles=50, num_iterations=100)
        fitness, positions = optimizer.return_ranked_populations()  # at end of PSO for all # particles, rank by cost function value
        hist_all = optimizer.all_history
        fit_all = optimizer.all_fitness
        np.save('position_pso', positions)  # param vectors for 1000 particles
        np.save('values_cost_pso', fitness)  # cost function for each iteration of 1000 particles
        np.save('his_all_pso', hist_all)
        np.save('fit_all_pso', fit_all)
        # print(fitness)
        print(optimizer.best) # what it takes here is the best set of parameter
        np.save('optimizer_best_10000_%s' % i, optimizer.best)
        # np.save('./optimizer.best/Geena_complete_model/optimizer_best_50_all_new_mil_'+str(now.strftime('%Y-%m-%d_%H%M')), optimizer.best)
        # if plot:
        #     display(optimizer.best)

# def display(parameter):
#     Y = np.copy(parameter)
#     param_values[rate_mask] = 10 ** Y
#     result = solver.run(param_values=param_values)
#     ysim_norm_1 = normalize(result.observables['MLKLa_obs'])
#
#     param_values[rate_mask] = 10 ** log10_original_values
#     result = solver.run(param_values=param_values)
#     ysim_norm_2 = normalize(result.observables['MLKLa_obs'])
#
#     plt.figure()
#     plt.plot(t, ysim_norm_1, '-^', linewidth=5, label='ODE optimization')
#     plt.plot(t, ysim_norm_2, '->', label='Initial parameters')
#     plt.plot(t, data, '-*', label='Experimental data')
#     plt.legend(loc=0)
#     plt.ylabel('concentration')
#     plt.xlabel('time (min)')
#     plt.tight_layout()
#     plt.savefig('./Figures/PSO/Geena/results'+str(now.strftime('%Y-%m-%d_%H%M')))
#     plt.show(block=True)
#     plt.close()

if '__main__' == __name__:
    run_example()
