"""
Example usage of the simplePSO package
"""
import matplotlib.pyplot as plt
import numpy as np

from pysb.examples.robertson import model
from pysb.simulator import ScipyOdeSimulator
from simplepso.logging import get_logger
from simplepso.pso import PSO

logger = get_logger()


# Defining a few helper functions to use
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)


def display(params, title=None):
    local_params = np.copy(params)
    param_values[rate_mask] = 10 ** local_params
    traj = solver.run(param_values=param_values)
    ysim_array = traj.dataframe[obs_names].values
    ysim_norm = normalize(ysim_array)

    plt.figure(figsize=(6, 4))
    if title is not None:
        plt.title(title)
    plt.plot(t, ysim_norm[:, 0], '^r', linestyle='--', label='A')
    plt.plot(t, ysim_norm[:, 1], 'ok', linestyle='--', label='C')
    plt.plot(t, norm_noisy_data_A, 'r-', label='Noisy A')
    plt.plot(t, norm_noisy_data_C, 'k-', label='Noisy C')
    plt.legend(loc=0)
    plt.ylabel('Normalized concentration')
    plt.xlabel('Time (s)')


'''
Here we define the cost function. We pass it the parameter vector, unlog it,
and pass it to the solver. We choose a chi square cost function, but you can 
provide any metric of you choosing
'''


def obj_function(params):
    # create copy of parameters
    params_tmp = np.copy(params)
    # convert back into regular base
    param_values[rate_mask] = 10 ** params_tmp
    traj = solver.run(param_values=param_values)
    ysim_array = traj.dataframe[obs_names].values
    ysim_norm = normalize(ysim_array)
    # chi^2 error
    err = np.sum((ydata_norm - ysim_norm) ** 2)
    if np.isnan(err):
        return 1000
    return err


# setup model and simulator
t = np.linspace(0, 50, 51)
obs_names = ['A_total', 'C_total']

# create pysb simulator instance
solver = ScipyOdeSimulator(
    model,
    t,
    integrator='lsoda',
    integrator_options={'rtol': 1e-8, 'atol': 1e-8}
)
traj = solver.run()

# Creating ideal data
ysim_array = traj.dataframe[obs_names].values
norm_data = normalize(ysim_array)

noise = 0.05
noisy_data_A = norm_data[:, 0] + np.random.uniform(-1 * noise, noise, len(t))
norm_noisy_data_A = normalize(noisy_data_A)

noisy_data_C = norm_data[:, 1] + np.random.uniform(-noise, noise, len(t))
norm_noisy_data_C = normalize(noisy_data_C)

ydata_norm = np.column_stack((norm_noisy_data_A, norm_noisy_data_C))

'''
In this example we are going to optimize the rate parameters of the model. 
We do this in log10 so we can span across large parameter spaces easily.
'''
rate_params = model.parameters_rules()
rate_mask = np.array([p in rate_params for p in model.parameters])
param_values = np.array([p.value for p in model.parameters])
log_original_values = np.log10(param_values[rate_mask])

if '__main__' == __name__:
    # We will use a best guess starting position for the model,
    start_position = log_original_values + \
                     np.random.uniform(-1, 1,
                                       size=len(log_original_values))

    display(start_position, "Before optimization")
    plt.tight_layout()
    plt.savefig("fit_before_pso.png", bbox_inches='tight')
    logger.info("Saving pre-fit figure as fit_before_pso.png")

    # Here we initial the class
    # We must proivde the cost function and a starting value
    optimizer = PSO(start=start_position, verbose=True, shrink_steps=False)

    # We also must set bounds of the parameter space, and the speed PSO will
    # travel (max speed in either direction)
    optimizer.set_bounds(parameter_range=4)
    optimizer.set_speed(speed_min=-.05, speed_max=.05)

    # Now we run the pso algorithm
    optimizer.run(num_particles=50, num_iterations=500, num_processors=12,
                  cost_function=obj_function, max_iter_no_improv=25)

    best_params = optimizer.best.pos
    display(best_params, "After optimization")
    plt.tight_layout()
    plt.savefig("fit_after_pso.png", bbox_inches='tight')
    logger.info("Saving post-fit figure as fit_after_pso.png")
    plt.show()
