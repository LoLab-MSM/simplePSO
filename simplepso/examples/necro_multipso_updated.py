try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass

import numpy as np
from necro_uncal_new import model
# from pysb.integrate import Solver
import scipy.interpolate
from pysb.integrate import *

from simplepso.pso import PSO
import collections

new_start =  np.load('optimizer_best_100_100_9_20_necromulti_pso3.npy')

model.enable_synth_deg()
obs_names = ['MLKLa_obs']
mlkl_obs = 'MLKLa_obs'


# Defining a few helper functions to use
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)


def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return np.vstack([recarray[name] for name in names]).T



x100 = np.array([30, 90, 270, 480, 600, 720, 840, 960])
y100 = np.array([0.00885691708746097,0.0161886154261265,0.0373005242261882,0.2798939020159581,0.510, .7797294067, 0.95,1])

# x10 = np.array([.5, 1.5, 4.5, 8, 10, 12, 14, 16])
y10 = np.array([0.0106013664572332,0.00519576571714913,0.02967443048221,0.050022163974868,0.108128107774737, 0.25,0.56055140114867, 0.77])

t = np.linspace(0, 960, 7)
solver1 = ScipyOdeSimulator(model, tspan=x100)
# x100 = np.array([0., .5, 1.5, 4.5, 8, 12, 16])
# y100 = np.array([0.,0.0088569170874609,0.0161886154261265,0.0373005242261882,0.07989390201595810, .639729406776,1])
#
# x10 = np.array([0., .5, 1.5, 4.5, 8, 12, 16])
# y10 = np.array([0., 0.0106013664572332,0.00519576571714913,0.02967443048221,0.050022163974868,0.198128107774737,0.56055140114867])

#
#
# #make an array for each of the kd made up data for mlklp
# x100 = np.array([0., 0.5, 1.5, 4.5, 8., 12., 16.])
# y100 = np.array([0., 0.008, 0.016, 0.037, 0.079, 0.639, 1.0])
#
# # x10 = np.array([0,2,4,6,8,10,12,14,16,18, 20])
# # y10 = np.array([0., 0.001, 0.02, 0.03, 0.04, 0.06, .09, .21, .40, .65, .81])
#
# x10 = np.array([0, .5, 1.5, 4.5, 8, 12, 16])
# # y10 = np.array([0, 0.001, 0.02, .09, .15, .41, 0.641])
#
# #
# y10 = np.array([0., 0.0106013664572332,
# 0.00519576571714913,
# 0.02967443048221,
# 0.050022163974868,
# 0.198128107774737,
# 0.640551401148672])
# y1 = np.array([0., 0.0060632926030143,
# 0.00942691670200054,
# 0.0113342231044983,
# 0.0312716821493274,
# 0.10956134602783,
# 0.36105885980804])
#
# #A20 data switching at 3 hours
# x1 = np.array([0, 0.5, 1.5, 4.5, 8., 12., 16.])
# y1 = np.array([0, 0.006, 0.009, 0.011, 0.0313, 0.110, 0.361])

data = collections.OrderedDict([('y100', y100), ('y10', y10)])
# data = collections.OrderedDict([('wt', wty), ('fd', fdy)])
# data = collections.OrderedDict()
# data = {'wt': wty, 'a20': a20y, 'td': tdy, 'fd': fdy, 'c8': c8y}
# data = collections.OrderedDict(sorted(data.items(), key = lambda t:t[1]))


# ydata_norm =

rate_params = model.parameters_rules()
# print(len(rate_params))
param_values = np.array([p.value for p in model.parameters])
# print(len(param_values))
rate_mask = np.array([p in rate_params for p in model.parameters])
# print(len(rate_mask))
# quit()

original_values = np.array([p.value for p in model.parameters])

# We search in log10 space for the parameters
log10_original_values = np.log10(original_values[rate_mask])

#Defining some functions to help plot the output of the parameters
def display(parameter_2):

    Y = np.copy(parameter_2)
    param_values[rate_mask] = 10 ** Y
    # rate_params = 10 ** Y

    x1_params = np.copy(param_values)
    x1_params[0] = 233
    ko_pars = [param_values, x1_params]

    result = solver1.run(param_values=ko_pars)

    ysim_array11 = result.observables[0]['MLKLa_obs']
    ysim_array22 = result.observables[1]['MLKLa_obs']

    # ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm11 = normalize(ysim_array11)
    ysim_norm22 = normalize(ysim_array22)

    mlkl_10 = np.array([0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
    mlkl_1 = np.array([0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])

    var_data = collections.OrderedDict([('var_00', mlkl_10), ('var_1', mlkl_1)])

    ysim = collections.OrderedDict([('sim_10', ysim_norm11), ('sim_1', ysim_norm22)])

    c = ['red', 'green', 'black','purple', 'cyan']
    d = ['red', 'green', 'black','purple', 'cyan']
    # colors = ['red', 'green', 'black', 'purple', 'orange']
    for i,j,k,l,m in zip(data, ysim,c, d, var_data):
        plt.figure()
        plt.plot(t / 60, wty,'--', color='black', label='Ctrl Mlklp')
        # plt.errorbar(t/60,data[i] , yerr=var_data[m],fmt='o', capsize=5, label='Noisy Mlklp{}'.format(i))
        plt.fill_between(t/60, data[i] - var_data[m], data[i] + var_data[m], facecolor='purple', label = '{}'.format(m), interpolate=True, alpha=.3)

        # plt.plot(t, norm_noisy_data_C, label='Noisy C')
        plt.plot(t/60, ysim[j],color = 'black',label='Best fit Mlklp {}'.format(j))
        # plt.plot(t, ysim_norm_2[:, 1], 'o', label='Best fit C')
        plt.legend(loc=0)
        plt.ylabel('molecules/cell')
        plt.xlabel('time (hrs)')
        plt.tight_layout()
        plt.savefig('necroptosis_pmlkl_25_100_start.png', format='png')
    plt.show()
    plt.close()

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# We choose a chi square cost function, but you can provide any metric of you choosing
# It must return a tuple
def obj_function(params):
    # print('obj function')
    # Y = np.copy(parameter_2)
    # param_values[rate_mask] = 10 ** Y
    params_tmp = np.copy(params)
    rate_params = 10 ** params_tmp #don't need to change
    param_values[rate_mask] = 10 ** params_tmp  # don't need to change
    #make a new parameter value set for each of the KD
    x1_params = np.copy(param_values)
    x1_params[0] = 23
    ko_pars = [x1_params, param_values]

    result = solver1.run(param_values=ko_pars)

    ysim_array11 = result.observables[0]['MLKLa_obs']
    ysim_array22 = result.observables[1]['MLKLa_obs']

    # ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm11 = normalize(ysim_array11)
    ysim_norm22 = normalize(ysim_array22)

    # mlkl_10 = np.array([0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
    # mlkl_1 = np.array([0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
    # mlkl_10 = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    # mlkl_1 = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

    # e1 = np.sum((y10 - ysim_norm11) ** 2 / (mlkl_10))
    # e2 = np.sum((y1 - ysim_norm22) ** 2 / (mlkl_1))
    #
    e1 = np.sum((y100 - ysim_norm11) ** 2)
    e2 = np.sum((y10 - ysim_norm22) ** 2)

    error = e1 + e2
    return error,


def run_example():
    # print('run_example')
    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(cost_function=obj_function,start = new_start, verbose=True)
    # We also must set bounds. This can be a single scalar or an array of len(start_position)
    optimizer.set_bounds(parameter_range=2)
    optimizer.set_speed(speed_min=-.25, speed_max=.25)
    optimizer.run(num_particles=75, num_iterations=100)
    print(optimizer.best)
    np.save('optimizer_best_100_100_9_20_necromulti_pso4',optimizer.best)

if '__main__' == __name__:
    run_example()
