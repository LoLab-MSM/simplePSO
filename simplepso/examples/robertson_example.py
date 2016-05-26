import matplotlib.pyplot as plt
import numpy as np
from pysb.examples.robertson import model
from pysb.integrate import Solver

from simplepso.pso import PSO


obs_names = ['A_total', 'C_total']


t = np.linspace(0, 40, 100)
solver = Solver(model, t, integrator='vode', rtol=1e-8, atol=1e-8)
solver.run()

# Defining a few helper functions to use
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)


def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return np.vstack([recarray[name] for name in names]).T


ysim_array = extract_records(solver.yobs, obs_names)
norm_data = normalize(ysim_array)
plt.plot(t, norm_data)
plt.legend(['A_Total', 'C_Total'], loc=0)
plt.show()


noisy_data_A = ysim_array[:, 0] + np.random.uniform(-0.05, 0.05, np.shape(ysim_array[:, 0]))
norm_noisy_data_A = normalize(noisy_data_A)
noisy_data_C = ysim_array[:, 1] + np.random.uniform(-.02, .02, np.shape(ysim_array[:, 1]))
norm_noisy_data_C = normalize(noisy_data_C)
ydata_norm = np.column_stack((norm_noisy_data_A, norm_noisy_data_C))

plt.plot(t, norm_noisy_data_A)
plt.plot(t, norm_noisy_data_C)
plt.plot(t, norm_data)
plt.legend(['A_total_noisy', 'C_total_noisy', 'A_total', 'C_total'], loc=0)
plt.show()

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])

original_values = np.array([p.value for p in model.parameters])

# We search in log10 space for the parameters
log10_original_values = np.log10(original_values[rate_mask])

# We will use a best guess starting position for the model, up or down 1 order of magnitude
start_position = log10_original_values + [1, -1,.75]# np.random.uniform(-1.0, 1.0, size=np.shape(log10_original_values))

# Defining some functions to help plot the output of the parameters
def display(parameter_1, parameter_2):
    Y = np.copy(parameter_1)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    ysim_array_1 = extract_records(solver.yobs, obs_names)
    ysim_norm_1 = normalize(ysim_array_1)
    Y = np.copy(parameter_2)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    ysim_array_2 = extract_records(solver.yobs, obs_names)
    ysim_norm_2 = normalize(ysim_array_2)

    plt.figure()
    plt.subplot(211)
    plt.title('Starting position')
    plt.plot(t, ysim_norm_1[:, 0], label='A')
    plt.plot(t, ysim_norm_1[:, 1], label='C')
    plt.plot(t, norm_noisy_data_A, label='Noisy A')
    plt.plot(t, norm_noisy_data_C, label='Noisy C')
    plt.subplot(212)
    plt.title('PSO best position')
    plt.plot(t, ysim_norm_2[:, 0], label='A')
    plt.plot(t, ysim_norm_2[:, 1], label='C')
    plt.plot(t, norm_noisy_data_A, label='Noisy A')
    plt.plot(t, norm_noisy_data_C, label='Noisy C')
    plt.legend(loc=0)
    plt.ylabel('concentration')
    plt.xlabel('time (s)')
    plt.tight_layout()
    plt.show()
    plt.close()

# Here we define the cost function. We pass it the parameter vector, unlog it, and pass it to the solver.
# We choose a chi square cost function, but you can provide any metric of you choosing
# It must return a tuple
def obj_function(params):
    params_tmp = np.copy(params)
    param_values[rate_mask] = 10 ** params_tmp
    solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    err = np.sum((ydata_norm - ysim_norm) ** 2)
    return err,

# Here we initial the class
# We must proivde the cost function and a starting value
optimizer = PSO(cost_function=obj_function, start=start_position, verbose=True)
# We also must set bounds. This can be a single scalar or an array of len(start_position)
optimizer.set_bounds(parameter_range=1)

optimizer.run(num_particles=20, num_iterations=50)

display(start_position, optimizer.best)
print(log10_original_values ** 10)
print(start_position ** 10)
print(optimizer.best ** 10)
