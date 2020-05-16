try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass

import numpy as np
from necro_uncal_new_updated import model
# from pysb.integrate import Solver
import scipy.interpolate
from pysb.integrate import *

from simplepso.pso import PSO
import collections


# new_start =  np.load('values_cost_pso.npy')
# print(new_start)
# quit()

model.enable_synth_deg()
obs_names = ['MLKLa_obs']
mlkl_obs = 'MLKLa_obs'


# Defining a few helper functions to use
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)


# def extract_records(recarray, names):
#     """Convert a record-type array and list of names into a float array"""
#     return np.vstack([recarray[name] for name in names]).T


# t = np.linspace(0, 960, 7)
t = [0, 30, 90, 270, 480,600, 720, 840, 960]
solver1 = ScipyOdeSimulator(model, tspan=t)

#make an array for each of the kd made up data for mlklp
#switching at 5 hours
# wtx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
# wty = np.array([0., 0., 0., 0.10, 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.,1.])

# x10 = np.array([0,2,4,6,8,10,12,14,16,18, 20])
# y10 = np.array([0., 0.001, 0.02, 0.03, 0.04, 0.06, .09, .21, .40, .65, .81])
#
# x100 = np.array([0., 0.5, 1.5, 4.5, 8., 12., 16.])
# # y100 = np.array([0., 0.008, 0.016, 0.037, 0.079, 0.639, 1.0])
# y100 = np.array([0., 0.01, 0.0152, 0.03, 0.437, 0.693, 0.99])


y10 = np.array([0, 0.0106013664572332,0.00519576571714913,0.02967443048221,0.050022163974868,
                0.088128107774737, 0.17, 0.30055140114867, 0.47])

y100 = np.array([0, 0.00885691708746097,0.0161886154261265,0.0373005242261882,
                  0.2798939020159581, 0.510, .8097294067, 0.95,0.98])
ydata_norm = y100

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

# We will use a best guess starting position for the model, up or down 1 order of magnitude
#start_position = log10_original_values #+ np.random.uniform(-3., 3., size=np.shape(log10_original_values)) #[-1.5, 1.3,
                                          #-.75]  # np.random.uniform(-1.5, 1.5, size=np.shape(log10_original_values))

# Defining some functions to help plot the output of the parameters
# def display(parameter_2):
#     # Y = np.copy(parameter_1)
#     # param_values[rate_mask] = 10 ** Y
#     # solver.run(param_values)
#     # ysim_array_1 = extract_records(solver.yobs, obs_names)
#     # ysim_norm_1 = normalize(ysim_array_1)
#     Y = np.copy(parameter_2)
#     param_values[rate_mask] = 10 ** Y
#     # rate_params = 10 ** Y
#     result = solver1.run(param_values=param_values)
#
#     ysim_array11 = result.observables['MLKLa_obs'][:]
#
#     # ysim_array = extract_records(solver.yobs, obs_names)
#     ysim_norm11 = normalize(ysim_array11)
#
#     mlkl_wt = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
#
#     var_data = collections.OrderedDict([('wt_var', mlkl_wt), ('a20_var', mlkl_a20), ('td_var', mlkl_td), ('fd_var', mlkl_fd),('c8_var', mlkl_c8)])
#
#     ysim = collections.OrderedDict([('wt_sim', ysim_norm11), ('a20_sim', ysim_norm22), ('td_sim', ysim_norm33), ('fd_sim', ysim_norm44),('c8_sim', ysim_norm55)])
#     # ysim = collections.OrderedDict([('wt_sim', ysim_norm11), ('fd_sim', ysim_norm44)])
#     # ysim = collections.OrderedDict(sorted(ysim.items(), key=lambda t: t[1]))
#
#
#     # solver1.run(param_values)
#     # ysim_array_2 = solver.yobs['MLKLa_obs']
#     # ysim_norm_2 = normalize(ysim_array_2)
#     # ['red', 'green', 'black', 'purple', 'orange']
#     # param_values[rate_mask] = 10 ** log10_original_values
#     # solver.run(param_values)
#     # ysim_array_3 = solver.yobs['MLKLa_obs']
#     # ysim_norm_3 = normalize(ysim_array_3)
#
#     # colors = [cmap(i) for i in np.linspace(0, 1, 5)]
#
#     c = ['red', 'green', 'black','purple', 'cyan']
#     d = ['red', 'green', 'black','purple', 'cyan']
#     # colors = ['red', 'green', 'black', 'purple', 'orange']
#     for i,j,k,l,m in zip(data, ysim,c, d, var_data):
#         plt.figure()
#         # plt.plot(t, ysim_norm_3[:, 0], '-^', linewidth=5, label='Ideal P')
#         # plt.plot(t, ysim_norm_3[:, 1], '-^', linewidth=5, label='Ideal C')
#         # plt.plot(t, ysim_norm_1[:, 0], '->', label='Starting P')
#         # plt.plot(t, ysim_norm_1[:, 1], '->', label='Starting C')
#         # plt.plot(t/60, data[i], color = 'purple',label='Noisy Mlklp{}'.format(i))
#         plt.plot(t / 60, wty,'--', color='black', label='Ctrl Mlklp')
#         # plt.errorbar(t/60,data[i] , yerr=var_data[m],fmt='o', capsize=5, label='Noisy Mlklp{}'.format(i))
#         plt.fill_between(t/60, data[i] - var_data[m], data[i] + var_data[m], facecolor='purple', label = '{}'.format(m), interpolate=True, alpha=.3)
#
#         # plt.plot(t, norm_noisy_data_C, label='Noisy C')
#         plt.plot(t/60, ysim[j],color = 'black',label='Best fit Mlklp {}'.format(j))
#         # plt.plot(t, ysim_norm_2[:, 1], 'o', label='Best fit C')
#         plt.legend(loc=0)
#         plt.ylabel('molecules/cell')
#         plt.xlabel('time (hrs)')
#         plt.tight_layout()
#         plt.savefig('necroptosis_kds_all_75_5000_start2l.png', format='png')
#     plt.show()
#     plt.close()

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# We choose a chi square cost function, but you can provide any metric of you choosing
# It must return a tuple
def obj_function(params):
    # print('obj function')
    # Y = np.copy(parameter_2)
    # param_values[rate_mask] = 10 ** Y
    params_tmp = np.copy(params)
    # rate_params = 10 ** params_tmp #don't need to change
    param_values[rate_mask] = 10 ** params_tmp  # don't need to change

    result = solver1.run(param_values=param_values)

    ysim_array1 = result.observables['MLKLa_obs'][:]

    ysim_norm1 = normalize(ysim_array1)

    #mlkl_wt = np.array([0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
    #mlkl_10 = np.array([0.001, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
    #e1 = np.sum(((ydata_norm - ysim_norm1) ** 2)/(mlkl_10))
    mlkl_wt = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05])

    e1 = np.sum(((ydata_norm - ysim_norm1) ** 2)/(mlkl_wt))
	
    return e1,

# def run_example():
#     best_pars = np.zeros((10000, len(model.parameters)))
#

def run_example2():
    # print('run_example')
    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(cost_function=obj_function, start=log10_original_values, verbose=True)
    # We also must set bounds. This can be a single scalar or an array of len(start_position)
    optimizer.set_bounds(parameter_range=2)
    optimizer.set_speed(speed_min=-.25, speed_max=.25)
    optimizer.run(num_particles=100, num_iterations=100)
    print(optimizer.best)
    np.save('optimizer_best_100_100_5_6_100tnf', optimizer.best)
    # print('whatever')
    # if plot:
    # display(optimizer.best)

if '__main__' == __name__:
    run_example2()
    # counter = 0
    # for i in range(10000):
    #      pso = PSO(save_sampled=False, verbose=False, num_proc=8)
    #      pso.set_cost_function(obj_function)
    #      pso.set_start_position(log10_original_values)
    #      pso.set_bounds(2)
    #      pso.set_speed(-.25, .25)
    #      pso.run(1000, 500)
    #      # np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples/pso_vals', pso.values) #cost function for PSO runs
    #      # ranked_values = PSO.return_ranked_populations() #at end of PSO for all # particles, rank by cost function value
    #      #
    #      # np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples', )
    #
    #      Y=np.copy(pso.best)
    #      param_values[rate_mask] = 10 ** Y
    #      print(param_values)
    #      if pso.values.min() < 10.0:
    #         best_pars[counter,:] = param_values
    #         print(best_pars[0:10,:])
    #         counter += 1
    #      print (i, counter)

#     np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples', best_pars)
# #
# def run_example():
#     pso = PSO(verbose=True, save_sampled=True)
#     pso.set_cost_function(obj_function)
#     pso.set_start_position(log10_original_values)
#     pso.set_bounds(lower=2, upper=2)
#     pso.set_speed(-.25, .25)
#     pso.run(25, 5)
#     fitness,positions = pso.return_ranked_populations()  # at end of PSO for all # particles, rank by cost function value
#     hist_all = pso.all_history
#     fit_all = pso.all_fitness
#     # print('fitness, particles')
#     # print(fitness)
#     # print('position rank')
#     # print(positions)
#     # print('all history')
#     # print(hist_all)
#     # print('fit all iterations')
#     # print(fit_all)
#     if plot:
#         display(pso.best)
#     display(pso.best)
#     np.save('position_pso', positions) # param vectors for 1000 particles
#     np.save('values_cost_pso', fitness) #cost function for each iteration of 1000 particles
#     np.save('his_all_pso', hist_all)
#     np.save('fit_all_pso', fit_all)
#
#
# if __name__ == '__main__':
#     run_example()


# def run_example():
#     print('run_example')
#     best_pars = np.zeros((100, len(model.parameters)))
#
#     counter = 0
#     # Here we initial the class
#     # We must proivde the cost function and a starting value
#     for i in range(10):
#         optimizer = PSO(cost_function=obj_function,start = log10_original_values, verbose=True)
#         # We also must set bounds. This can be a single scalar or an array of len(start_position)
#         optimizer.set_bounds(parameter_range=2)
#         optimizer.set_speed(speed_min=-.25, speed_max=.25)
#         optimizer.run(num_particles=25, num_iterations=25)
#         best_pars[i] = optimizer.best
#         print(optimizer.best)
#
#     np.save('optimizer_best_75_50_100TNF',best_pars)
#
#
# # def run_example2():
# #     best_pars = np.zeros((10, len(model.parameters)))
# #
# #     counter = 0
# #     for i in range(10):
# #          pso = PSO(save_sampled=False, verbose=False, num_proc=8)
# #          pso.set_cost_function(obj_function)
# #          pso.set_start_position(log10_original_values)
# #          pso.set_bounds(2)
# #          pso.set_speed(-.25, .25)
# #          pso.run(num_particles=10, num_iterations=50)
# #          # np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples/pso_vals', pso.values) #cost function for PSO runs
# #          # ranked_values = PSO.return_ranked_populations() #at end of PSO for all # particles, rank by cost function value
# #          #
# #          # np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples', )
# #
# #          # Y=np.copy(pso.best)
# #          # param_values[rate_mask] = 10 ** Y
# #          best_pars[i] = pso.best
# #          # print(param_values)
# #          # if pso.values.min() < 10.0:
# #          #    best_pars[counter,:] = param_values
# #          #    print(best_pars[0:10,:])
# #          #    counter += 1
# #          print (i, counter)
# #
# #     np.save('10_pso_nfkb_marco', best_pars)
#
# if '__main__' == __name__:
#     run_example()
