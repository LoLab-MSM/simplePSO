try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass

import numpy as np
from necro import model
from pysb.integrate import Solver
import scipy.interpolate

from simplepso.pso import PSO

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
solver = Solver(model, t, integrator='lsoda', rtol=1e-6, atol=1e-6)
solver.run()

# Creating ideal data
ysim_array = extract_records(solver.yobs, obs_names)
norm_data = normalize(ysim_array)

x = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
y = np.array([0., 0., 0., 0., 0., 0.5, 1., 1., 1., 1., 1., 1., 1.])

mlkl_obs_total = model.parameters['MLKLa_0'].value
mlkl_data = np.array([9810.0, 180.0, mlkl_obs_total])

# noisy_data_A = ysim_array[:, 0]
# norm_noisy_data_A = normalize(noisy_data_A) + np.random.uniform(-0.1, 0.1, np.shape(ysim_array[:, 0]))
ydata_norm = y

# np.random.seed(0)

# Creating noisy data
# prod = np.array([ 0.,   7.,  13.,  24.,  30.,  39.,  43.,  46.,  49.,  57.,  59.,  63.,  65.,  68.,  73.,
#           77.,  79.,  82.,  83.,  85.,  88.,  88.,  88.,  88.,  90.,  90.,  91.,  91.,  91.,  91.,
#           92.,  93.,  93.,  93.,  94.,  95.,  96.,  97.,  98.,  98.,  99.,  99.,  99.,  99.,  99.,
#           99.,  99.,  99.,  99.,  99.,  99.])
# prod_norm = normalize(prod)
# ydata_norm = prod_norm

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])

original_values = np.array([p.value for p in model.parameters])

# We search in log10 space for the parameters
log10_original_values = np.log10(original_values[rate_mask])

# We will use a best guess starting position for the model, up or down 1 order of magnitude
start_position = log10_original_values + np.random.uniform(-3., 3., size=np.shape(log10_original_values)) #[-1.5, 1.3,
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
    solver.run(param_values)
    ysim_array_2 = solver.yobs['MLKLa_obs']
    ysim_norm_2 = normalize(ysim_array_2)

    # param_values[rate_mask] = 10 ** log10_original_values
    # solver.run(param_values)
    # ysim_array_3 = solver.yobs['MLKLa_obs']
    # ysim_norm_3 = normalize(ysim_array_3)

    plt.figure()
    plt.subplot(111)
    # plt.plot(t, ysim_norm_3[:, 0], '-^', linewidth=5, label='Ideal P')
    # plt.plot(t, ysim_norm_3[:, 1], '-^', linewidth=5, label='Ideal C')
    # plt.plot(t, ysim_norm_1[:, 0], '->', label='Starting P')
    # plt.plot(t, ysim_norm_1[:, 1], '->', label='Starting C')
    plt.plot(t, ydata_norm, label='Noisy Mlklp')
    # plt.plot(t, norm_noisy_data_C, label='Noisy C')
    plt.plot(t, ysim_norm_2, 'o', label='Best fit Mlklp')
    # plt.plot(t, ysim_norm_2[:, 1], 'o', label='Best fit C')
    plt.legend(loc=0)
    plt.ylabel('molecules/cell')
    plt.xlabel('time (min)')
    plt.tight_layout()
    plt.savefig('necroptosis.png', format = 'png')
    plt.show()
    plt.close()

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# We choose a chi square cost function, but you can provide any metric of you choosing
# It must return a tuple
def obj_function(params):
    params_tmp = np.copy(params)
    param_values[rate_mask] = 10 ** params_tmp
    solver.run(param_values)
    ysim_array = solver.yobs['MLKLa_obs']
    # ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    mlkl_var = np.var(y)
    mlkl_v = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])	
    e1 = np.sum((ydata_norm - ysim_norm) ** 2 / (mlkl_v))

#    st, sc, sk = scipy.interpolate.splrep(t, ysim_norm)
#    t10 = scipy.interpolate.sproot((st, sc - 0.10, sk))[0]
#    t90 = scipy.interpolate.sproot((st, sc - 0.90, sk))[0]
#    td = (t10 + t90) / 2
#    ts = t90 - t10
#    yfinal = ysim_array[-1]
#    mlkl_sim = [td, ts, yfinal]
#    e2 = np.sum((mlkl_data - mlkl_sim) ** 2)

#    error = e1 + e2
    return e1,


def run_example():
    # print('run_example')
    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(cost_function=obj_function, start= log10_original_values, verbose=True)
    # We also must set bounds. This can be a single scalar or an array of len(start_position)
    optimizer.set_bounds(parameter_range=3)
    optimizer.set_speed(speed_min=-.25, speed_max=.25)
    optimizer.run(num_particles=50, num_iterations=100)
    # print('whatever')
    if plot:
	 display(optimizer.best)
    #
    #     print("Original values {0}".format(log10_original_values ** 10))
    #     print("Starting values {0}".format(start_position ** 10))
    #     print("Best PSO values {0}".format(optimizer.best ** 10))
    #     fig = plt.figure()
    #     fig.add_subplot(221)
    #     plt.scatter(log10_original_values[0], log10_original_values[1], marker='>', color='b', label='ideal')
    #     plt.scatter(start_position[0], start_position[1], marker='^', color='r', label='start')
    #     plt.scatter(optimizer.history[:, 0], optimizer.history[:, 1], c=optimizer.values, cmap=plt.cm.coolwarm)
    #
    #     fig.add_subplot(223)
    #     plt.scatter(log10_original_values[0], log10_original_values[2], marker='>', color='b', label='ideal')
    #     plt.scatter(start_position[0], start_position[2], marker='^', color='r', label='start')
    #     plt.scatter(optimizer.history[:, 0], optimizer.history[:, 2], c=optimizer.values, cmap=plt.cm.coolwarm)
    #
    #     fig.add_subplot(222)
    #     plt.scatter(log10_original_values[1], log10_original_values[2], marker='>', color='b', label='ideal')
    #     plt.scatter(start_position[1], start_position[2], marker='^', color='r', label='start')
    #     plt.scatter(optimizer.history[:, 1], optimizer.history[:, 2], c=optimizer.values, cmap=plt.cm.coolwarm)
    #
    #     fig.add_subplot(224)
    #     plt.legend(loc=0)
    #     plt.colorbar()
    #     plt.tight_layout()
    #     plt.savefig('population_necro.png')
    #     plt.show()


if '__main__' == __name__:
    run_example()
