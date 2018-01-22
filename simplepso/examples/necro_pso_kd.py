try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False
    pass

import numpy as np
from correct_necro_molecules import model
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
solver1 = Solver(model, t, integrator='lsoda', rtol=1e-6, atol=1e-6)
solver2 = Solver(model, t, integrator='lsoda', rtol=1e-6, atol=1e-6)
solver3 = Solver(model, t, integrator='lsoda', rtol=1e-6, atol=1e-6)
solver4 = Solver(model, t, integrator='lsoda', rtol=1e-6, atol=1e-6)
solver5 = Solver(model, t, integrator='lsoda', rtol=1e-6, atol=1e-6)

solver1.run()
solver2.run()
solver3.run()
solver4.run()
solver5.run()

# Creating ideal data
ysim_array = extract_records(solver1.yobs, obs_names)
norm_data = normalize(ysim_array)


#make an array for each of the kd made up data for mlklp
wtx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
wty = np.array([0., 0., 0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.])

#A20 data
a20x = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
a20y = np.array([0., 0., 0.25, 0.5, 0.75, 1., 1., 1., 1., 1., 1., 1., 1.])

#Tradd data
tdx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
tdy = np.array([0., 0., 0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.])

#Fadd Data
fdx = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
fdy = np.array([0., 0., 0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.])

#C8 Data
c8x = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
c8y = np.array([0., 0., 0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1., 1., 1.])



mlkl_obs_total = model.parameters['MLKLa_0'].value
mlkl_data = np.array([2000.0, 180.0, mlkl_obs_total])

# noisy_data_A = ysim_array[:, 0]
# norm_noisy_data_A = normalize(noisy_data_A) + np.random.uniform(-0.1, 0.1, np.shape(ysim_array[:, 0]))
ydata_norm = wty

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
    plt.plot(t, ysim_norm_2, label='Best fit Mlklp')
    # plt.plot(t, ysim_norm_2[:, 1], 'o', label='Best fit C')
    plt.legend(loc=0)
    plt.ylabel('molecules/cell')
    plt.xlabel('time (min)')
    plt.tight_layout()
    plt.savefig('necroptosis_30_new.png', format = 'png')
    plt.show()
    plt.close()

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# We choose a chi square cost function, but you can provide any metric of you choosing
# It must return a tuple
def obj_function(params):
    params_tmp = np.copy(params)
    param_values[rate_mask] = 10 ** params_tmp #don't need to change
    #make a new parameter value set for each of the KD
    a20_params = np.copy(param_values)
    a20_params[6] = 2700
    tradd_params = np.copy(param_values)
    tradd_params[2] = 2700
    fadd_params = np.copy(param_values)
    fadd_params[8] = 2424
    c8_params = np.copy(param_values)
    c8_params[11] = 2700

    solver1.run(param_values)
    solver2.run(a20_params)
    solver3.run(tradd_params)
    solver4.run(fadd_params)
    solver5.run(c8_params)

    # list = [y1, y2, y3, y4, y5]
    # for i in list:
    ysim_array1 = solver1.yobs['MLKLa_obs']
    ysim_array2 = solver2.yobs['MLKLa_obs']
    ysim_array3 = solver3.yobs['MLKLa_obs']
    ysim_array4 = solver4.yobs['MLKLa_obs']
    ysim_array5 = solver5.yobs['MLKLa_obs']

    # ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm1 = normalize(ysim_array1)
    ysim_norm2 = normalize(ysim_array2)
    ysim_norm3 = normalize(ysim_array3)
    ysim_norm4 = normalize(ysim_array4)
    ysim_norm5 = normalize(ysim_array5)

    # mlkl_var = np.var(y)
    mlkl_v = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    # mlkl_v = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    # mlkl_v = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    # mlkl_v = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    # mlkl_v = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

    e1 = np.sum((ydata_norm - ysim_norm1) ** 2 / (mlkl_v))
    e2 = np.sum((a20y - ysim_norm2) ** 2 )
    e3 = np.sum((tdy - ysim_norm3) ** 2 )
    e4 = np.sum((fdy - ysim_norm4) ** 2 )
    e5 = np.sum((c8y - ysim_norm5) ** 2 )

    error = e1 + e2 + e3 + e4 + e5
    return error,


def run_example():
    # print('run_example')
    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(cost_function=obj_function, start= log10_original_values, verbose=True)
    # We also must set bounds. This can be a single scalar or an array of len(start_position)
    optimizer.set_bounds(parameter_range=3)
    optimizer.set_speed(speed_min=-.25, speed_max=.25)
    optimizer.run(num_particles=25, num_iterations=100)
    print(optimizer.best)
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
