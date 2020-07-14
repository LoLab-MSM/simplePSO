import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import wasserstein_distance

from pysb.integrate import odesolve
from pysb.simulator import OpenCLSSASimulator
from simplepso import PSO


def run_params_for_plot(params):
    simulator.param_values = None
    simulator.initials = None
    traj = simulator.run(tspan, param_values=params, number_sim=num_sim)
    return traj.dataframe[name].unstack(0).values


def obj_function(traj_dist):
    tmp = traj_dist[name].T.unstack('simulation').values[-1, :]
    return wasserstein_distance(tmp, actual)


def run_pso():
    optimizer = PSO(start=log10_original_values, verbose=True)
    optimizer.set_bounds(parameter_range=.2)
    optimizer.set_speed(speed_min=-.05, speed_max=.05)
    optimizer.run_ssa(model, num_particles=num_particles, num_iterations=200,
                      cost_function=obj_function, num_sim=num_sim,
                      simulator=simulator)

    print("Ending parameters")
    print(10 ** optimizer.best.pos)
    return optimizer.best.pos


def add_subplot(traj, title, params):
    plt.title(title)
    plt.plot(tspan, traj, '0.5', lw=2, alpha=0.25)
    plt.plot(tspan, traj.mean(1), 'k-*', lw=3, label="Mean")
    plt.plot(tspan, traj.min(1), 'b--', lw=3, label="Minimum")
    plt.plot(tspan, traj.max(1), 'r--', lw=3, label="Maximum")
    y = odesolve(model, tspan, params)
    plt.plot(tspan, y[name], 'g--', lw=3, label="ODE")
    plt.ylim(0, 800)
    plt.xlabel('Time')
    plt.ylabel('X molecules')


def run():
    fit_params = run_pso()

    plt.figure(figsize=(9, 5))
    plt.subplot(231)
    add_subplot(actual_traj, "Training data", orig_values)

    plt.subplot(232)
    orig_values[rate_mask] = 10 ** noisy_start
    start_sim = run_params_for_plot(orig_values)

    add_subplot(start_sim, "Starting state", orig_values)

    orig_values[rate_mask] = 10 ** fit_params
    best_fit_sim = run_params_for_plot(orig_values)

    plt.subplot(233)
    add_subplot(best_fit_sim, "After training", orig_values)

    plt.subplot(212)
    sns.kdeplot(actual, label='actual')
    sns.kdeplot(start_sim[-1, :], label='before')
    sns.kdeplot(best_fit_sim[-1, :], label='after')

    plt.legend()
    plt.tight_layout()
    savename = 'trained_model'
    plt.savefig("{}.png".format(savename), bbox_inches='tight', dpi=300)
    plt.show()
    plt.close()


if __name__ == '__main__':
    from pysb.examples.schloegl import model

    tspan = np.linspace(0, 100, 101)
    model.parameters['X_0'].value = 500
    name = 'X_total'
    num_sim = 100
    num_particles = 8
    # uncomment to use CUDA device
    # simulator = CudaSSASimulator(model, tspan=tspan, verbose=False)

    simulator = OpenCLSSASimulator(model, tspan=tspan, verbose=False)

    actual_traj = run_params_for_plot(None)
    actual = actual_traj[-1, :]
    rate_mask = np.array(
        [p in model.parameters_rules() for p in model.parameters])

    orig_values = np.array([p.value for p in model.parameters])
    log10_original_values = np.log10(orig_values[rate_mask])
    n_params = len(log10_original_values)
    print("Ideal parameters")
    print(log10_original_values)

    noisy_start = log10_original_values + np.random.uniform(-.1, .1, n_params)

    run()
