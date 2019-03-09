try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass

import numpy as np
from necroptosismodule import model
# from pysb.integrate import Solver
import scipy.interpolate
from pysb.integrate import *

from simplepso.pso import PSO
import collections

new_start =  np.load('optimizer_best_5000_all_new_2.npy')

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

t = np.linspace(0, 720, 13)
solver1 = ScipyOdeSimulator(model, tspan=t)

#make an array for each of the kd made up data for mlklp
#switching at 5 hours
wtx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
wty = np.array([0., 0., 0., 0., 0.10, 0.25, 0.5, 0.75, 1., 1., 1., 1.,1.])

#A20 data switching at 3 hours
a20x = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
a20y = np.array([0.,0., 0.2, 0.5, 0.75, 1., 1., 1., 1., 1., 1., 1., 1.])

#Tradd data switching at 5
tdx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
tdy = np.array([0., 0., 0., 0.10, 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.,1.])

#Fadd Data switching at 4 hours
fdx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
fdy = np.array([0., 0., 0.10,0.25, 0.5, 0.75, 1., 1., 1., 1., 1., 1., 1.])

#C8 Data switching at 4
c8x = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
c8y = np.array([0., 0., 0.10,0.25, 0.5, 0.75, 1., 1., 1., 1., 1., 1.,1.])

data = collections.OrderedDict([('wt', wty), ('a20', a20y), ('td', tdy),('fd', fdy), ('c8', c8y)])
# data = collections.OrderedDict([('wt', wty), ('fd', fdy)])
# data = collections.OrderedDict()
# data = {'wt': wty, 'a20': a20y, 'td': tdy, 'fd': fdy, 'c8': c8y}
# data = collections.OrderedDict(sorted(data.items(), key = lambda t:t[1]))


ydata_norm = wty

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
def display(parameter_2):
    # Y = np.copy(parameter_1)
    # param_values[rate_mask] = 10 ** Y
    # solver.run(param_values)
    # ysim_array_1 = extract_records(solver.yobs, obs_names)
    # ysim_norm_1 = normalize(ysim_array_1)
    Y = np.copy(parameter_2)
    param_values[rate_mask] = 10 ** Y
    # rate_params = 10 ** Y

    a20_params = np.copy(param_values)
    a20_params[6] = 2700
    tradd_params = np.copy(param_values)
    tradd_params[2] = 2700
    fadd_params = np.copy(param_values)
    fadd_params[8] = 2424
    c8_params = np.copy(param_values)
    c8_params[11] = 2700
    ko_pars = [param_values, a20_params, tradd_params,fadd_params, c8_params]

    result = solver1.run(param_values=ko_pars)

    ysim_array11 = result.observables[0]['MLKLa_obs']
    ysim_array22 = result.observables[1]['MLKLa_obs']
    ysim_array33 = result.observables[2]['MLKLa_obs']
    ysim_array44 = result.observables[3]['MLKLa_obs']
    ysim_array55 = result.observables[4]['MLKLa_obs']

    # ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm11 = normalize(ysim_array11)
    ysim_norm22 = normalize(ysim_array22)
    ysim_norm33 = normalize(ysim_array33)
    ysim_norm44 = normalize(ysim_array44)
    ysim_norm55 = normalize(ysim_array55)


    mlkl_wt = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_a20= np.array([0.05, 0.05,0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_td = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_fd = np.array([0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_c8 = np.array([0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

    var_data = collections.OrderedDict([('wt_var', mlkl_wt), ('a20_var', mlkl_a20), ('td_var', mlkl_td), ('fd_var', mlkl_fd),('c8_var', mlkl_c8)])

    ysim = collections.OrderedDict([('wt_sim', ysim_norm11), ('a20_sim', ysim_norm22), ('td_sim', ysim_norm33), ('fd_sim', ysim_norm44),('c8_sim', ysim_norm55)])
    # ysim = collections.OrderedDict([('wt_sim', ysim_norm11), ('fd_sim', ysim_norm44)])
    # ysim = collections.OrderedDict(sorted(ysim.items(), key=lambda t: t[1]))


    # solver1.run(param_values)
    # ysim_array_2 = solver.yobs['MLKLa_obs']
    # ysim_norm_2 = normalize(ysim_array_2)
    # ['red', 'green', 'black', 'purple', 'orange']
    # param_values[rate_mask] = 10 ** log10_original_values
    # solver.run(param_values)
    # ysim_array_3 = solver.yobs['MLKLa_obs']
    # ysim_norm_3 = normalize(ysim_array_3)

    # colors = [cmap(i) for i in np.linspace(0, 1, 5)]

    c = ['red', 'green', 'black','purple', 'cyan']
    d = ['red', 'green', 'black','purple', 'cyan']
    # colors = ['red', 'green', 'black', 'purple', 'orange']
    for i,j,k,l,m in zip(data, ysim,c, d, var_data):
        plt.figure()
        # plt.plot(t, ysim_norm_3[:, 0], '-^', linewidth=5, label='Ideal P')
        # plt.plot(t, ysim_norm_3[:, 1], '-^', linewidth=5, label='Ideal C')
        # plt.plot(t, ysim_norm_1[:, 0], '->', label='Starting P')
        # plt.plot(t, ysim_norm_1[:, 1], '->', label='Starting C')
        # plt.plot(t/60, data[i], color = 'purple',label='Noisy Mlklp{}'.format(i))
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
        plt.savefig('necroptosis_kds_all_75_5000_start2l.png', format='png')
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
    # rate_params = 10 ** params_tmp #don't need to change
    param_values[rate_mask] = 10 ** params_tmp  # don't need to change
    # print(len(param_values[rate_mask]))
    # quit()
    #make a new parameter value set for each of the KD
    a20_params = np.copy(param_values)
    a20_params[6] = 2700
    tradd_params = np.copy(param_values)
    tradd_params[2] = 2700
    fadd_params = np.copy(param_values)
    fadd_params[8] = 2424
    c8_params = np.copy(param_values)
    c8_params[11] = 2700
    ko_pars = [param_values, a20_params, tradd_params,fadd_params, c8_params]

    result = solver1.run(param_values=ko_pars)
    # solver2.run(param_values=a20_params)
    # solver3.run(param_values=tradd_params)
    # solver4.run(param_values=fadd_params)
    # solver5.run(param_values=c8_params)

    # list = [y1, y2, y3, y4, y5]
    # for i in list:
    ysim_array1 = result.observables[0]['MLKLa_obs']
    ysim_array2 = result.observables[1]['MLKLa_obs']
    ysim_array3 = result.observables[2]['MLKLa_obs']
    ysim_array4 = result.observables[3]['MLKLa_obs']
    ysim_array5 = result.observables[4]['MLKLa_obs']

    # ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm1 = normalize(ysim_array1)
    ysim_norm2 = normalize(ysim_array2)
    ysim_norm3 = normalize(ysim_array3)
    ysim_norm4 = normalize(ysim_array4)
    ysim_norm5 = normalize(ysim_array5)

    # mlkl_var = np.var(y)
    mlkl_wt = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_a20= np.array([0.05, 0.05,0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_td = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_fd = np.array([0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    mlkl_c8 = np.array([0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

    e1 = np.sum((ydata_norm - ysim_norm1) ** 2 / (mlkl_wt))
    e2 = np.sum((a20y - ysim_norm2) ** 2 / (mlkl_a20))
    e3 = np.sum((tdy - ysim_norm3) ** 2 / (mlkl_td))
    e4 = np.sum((fdy - ysim_norm4) ** 2 / (mlkl_fd))
    e5 = np.sum((c8y - ysim_norm5) ** 2 / (mlkl_c8))

    error = e1 + e2 + e3 + e4 + e5
    return error,

# def run_example():
#     best_pars = np.zeros((10000, len(model.parameters)))
#
#     counter = 0
#     for i in range(10000):
#          pso = PSO(save_sampled=False, verbose=False, num_proc=8)
#          pso.set_cost_function(obj_function)
#          pso.set_start_position(log10_original_values)
#          pso.set_bounds(2)
#          pso.set_speed(-.25, .25)
#          pso.run(1000, 500)
#          # np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples/pso_vals', pso.values) #cost function for PSO runs
#          # ranked_values = PSO.return_ranked_populations() #at end of PSO for all # particles, rank by cost function value
#          #
#          # np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples', )
#
#          Y=np.copy(pso.best)
#          param_values[rate_mask] = 10 ** Y
#          print(param_values)
#          if pso.values.min() < 10.0:
#             best_pars[counter,:] = param_values
#             print(best_pars[0:10,:])
#             counter += 1
#          print (i, counter)
#
#     np.save('/home/ildefog/ParticleSwarmOptimization/simplepso/examples', best_pars)
#
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


def run_example():
    # print('run_example')
    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(cost_function=obj_function,start = new_start, verbose=True)
    # We also must set bounds. This can be a single scalar or an array of len(start_position)
    optimizer.set_bounds(parameter_range=2)
    optimizer.set_speed(speed_min=-.25, speed_max=.25)
    optimizer.run(num_particles=50, num_iterations=500)
    print(optimizer.best)
    np.save('optimizer_best_50_500_mar8',optimizer.best)
#     # print('whatever')
#     if plot:
# 	 display(optimizer.best)
#     #
#     #     print("Original values {0}".format(log10_original_values ** 10))
#     #     print("Starting values {0}".format(start_position ** 10))
#     #     print("Best PSO values {0}".format(optimizer.best ** 10))
#     #     fig = plt.figure()
#     #     fig.add_subplot(221)
#     #     plt.scatter(log10_original_values[0], log10_original_values[1], marker='>', color='b', label='ideal')
#     #     plt.scatter(start_position[0], start_position[1], marker='^', color='r', label='start')
#     #     plt.scatter(optimizer.history[:, 0], optimizer.history[:, 1], c=optimizer.values, cmap=plt.cm.coolwarm)
#     #
#     #     fig.add_subplot(223)
#     #     plt.scatter(log10_original_values[0], log10_original_values[2], marker='>', color='b', label='ideal')
#     #     plt.scatter(start_position[0], start_position[2], marker='^', color='r', label='start')
#     #     plt.scatter(optimizer.history[:, 0], optimizer.history[:, 2], c=optimizer.values, cmap=plt.cm.coolwarm)
#     #
#     #     fig.add_subplot(222)
#     #     plt.scatter(log10_original_values[1], log10_original_values[2], marker='>', color='b', label='ideal')
#     #     plt.scatter(start_position[1], start_position[2], marker='^', color='r', label='start')
#     #     plt.scatter(optimizer.history[:, 1], optimizer.history[:, 2], c=optimizer.values, cmap=plt.cm.coolwarm)
#     #
#     #     fig.add_subplot(224)
#     #     plt.legend(loc=0)
#     #     plt.colorbar()
#     #     plt.tight_layout()
#     #     plt.savefig('population_necro.png')
#     #     plt.show()
#
#
if '__main__' == __name__:
    run_example()
