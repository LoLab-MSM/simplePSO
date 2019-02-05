try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass

import numpy as np
from necro_mn_marco import model
# from pysb.integrate import Solver
import scipy.interpolate
from pysb.integrate import *
from pysb.simulator.bng import BngSimulator
from simplepso.pso import PSO
import collections

model.enable_synth_deg()
obs_names = ['MLKLu', 'RIP3u', 'TRADDu', 'FADDu', 'C8u']

# Defining a few helper functions to use
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)

def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return np.vstack([recarray[name] for name in names]).T

t = np.linspace(0, 1440, 2)
solver1 = BngSimulator(model, tspan=t)

#make an array for each of the kd made up data for mlklp
#switching at 5 hours
wtx = np.array([60., 1440.])
wty = np.array([0.53, 1.])

#A20 data switching at 3 hours
rip3x = np.array([60., 1440.])
rip3y = np.array([1., 0.91])

#Tradd data switching at 5
tdx = np.array([60., 1440.])
tdy = np.array([1., 0.39])

#Fadd Data switching at 4 hours
fdx = np.array([60., 1440.])
fdy = np.array([1., 0.30])

#C8 Data switching at 4
c8x = np.array([60., 1440.])
c8y = np.array([1., 0.50])

# data = collections.OrderedDict([('wt', wty), ('rip3', rip3y), ('td', tdy),('fd', fdy), ('c8', c8y)])

ydata_norm = wty
rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
original_values = np.array([p.value for p in model.parameters])

# We search in log10 space for the parameters
log10_original_values = np.log10(original_values[rate_mask])

# We will use a best guess starting position for the model, up or down 1 order of magnitude
#start_position = log10_original_values #+ np.random.uniform(-3., 3., size=np.shape(log10_original_values)) #[-1.5, 1.3,
                                          #-.75]  # np.random.uniform(-1.5, 1.5, size=np.shape(log10_original_values))

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# We choose a chi square cost function, but you can provide any metric of you choosing
# It must return a tuple
def obj_function(params):
    params_tmp = np.copy(params)
    # rate_params = 10 ** params_tmp #don't need to change
    param_values[rate_mask] = 10 ** params_tmp  # don't need to change
    # print(param_values)
    result = solver1.run(method= 'ode', param_values=param_values)

    ysim_array1 = result.observables['MLKLu'][:]
    ysim_array2 = result.observables['RIP3u'][:]
    ysim_array3 = result.observables['TRADDu'][:]
    ysim_array4 = result.observables['FADDu'][:]
    ysim_array5 = result.observables['C8u'][:]

    ysim_norm1 = normalize(ysim_array1)
    ysim_norm2 = normalize(ysim_array2)
    ysim_norm3 = normalize(ysim_array3)
    ysim_norm4 = normalize(ysim_array4)
    ysim_norm5 = normalize(ysim_array5)

    e1 = np.sum((ydata_norm - ysim_norm1) ** 2)
    e2 = np.sum((rip3y - ysim_norm2) ** 2)
    e3 = np.sum((tdy - ysim_norm3) ** 2)
    e4 = np.sum((fdy - ysim_norm4) ** 2)
    e5 = np.sum((c8y - ysim_norm5) ** 2)

    error = e1 + e2 + e3 + e4 + e5
    return error,

def run_example():
    # print('run_example')
    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(cost_function=obj_function,start = log10_original_values, verbose=True)
    # print('optimizing')
    # We also must set bounds. This can be a single scalar or an array of len(start_position)
    optimizer.set_bounds(parameter_range=2)
    optimizer.set_speed(speed_min=-.25, speed_max=.25)
    optimizer.run(num_particles=25, num_iterations=100)
    print(optimizer.best)
    np.save('optimizer_best_50_necro',optimizer.best)

if '__main__' == __name__:
    run_example()
