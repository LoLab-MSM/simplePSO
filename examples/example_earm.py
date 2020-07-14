# -*- coding: utf-8 -*-
"""
simplePSO example using the EARM 1.0 model
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate

from pysb.examples.earm_1_0 import model
from pysb.simulator import ScipyOdeSimulator
from simplepso.pso import PSO

# load experimental data
data_path = os.path.join(os.path.dirname(__file__), 'data',
                         'EC-RP_IMS-RP_IC-RP_data_for_models.csv')

exp_data = pd.read_csv(data_path, index_col=False)

# Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# yfinal (the last value of the IMS-RP trajectory)
momp_data = np.array([9810.0, 180.0, model.parameters['mSmac_0'].value])
momp_var = np.array([7245000.0, 3600.0, 1e4])

# timepoints for simulation. These must match the experimental data.
tspan = exp_data['Time'].values.copy()

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
starting_position = np.log10(param_values[rate_mask])

solver = ScipyOdeSimulator(
    model,
    tspan,
    integrator='lsoda',
    use_analytic_jacobian=True,
    compiler='python',
    integrator_options={"rtol": 1e-6, "atol": 1e-6}
)


def likelihood(position):
    param_values[rate_mask] = 10 ** position.copy()
    traj = solver.run(param_values=param_values)

    # normalize trajectories
    bid_traj = traj.observables['tBid_total'] / model.parameters['Bid_0'].value
    cparp_traj = traj.observables['CPARP_total'] / model.parameters[
        'PARP_0'].value
    momp_traj = traj.observables['cSmac_total']

    # calculate chi^2 distance for each time course
    e1 = np.sum((exp_data['norm_IC-RP'] - bid_traj) ** 2 /
                (2 * exp_data['nrm_var_IC-RP'])) / len(bid_traj)

    e2 = np.sum((exp_data['norm_EC-RP'] - cparp_traj) ** 2 /
                (2 * exp_data['nrm_var_EC-RP'])) / len(cparp_traj)

    # Here we fit a spline to find where we get 50% release of MOMP reporter
    if np.nanmax(momp_traj) == 0:
        print('No aSmac!')
        e3 = 10000000
    else:
        ysim_momp_norm = momp_traj / np.nanmax(momp_traj)
        st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
        try:
            t10 = scipy.interpolate.sproot((st, sc - 0.10, sk))[0]
            t90 = scipy.interpolate.sproot((st, sc - 0.90, sk))[0]
        except IndexError:
            t10 = 0
            t90 = 0

        # time of death  = halfway point between 10 and 90%
        td = (t10 + t90) / 2

        # time of switch is time between 90 and 10 %
        ts = t90 - t10

        # final fraction of aSMAC (last value)
        momp_sim = [td, ts, momp_traj[-1]]

        e3 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3

    return e1 + e2 + e3


def display(position, save_name):
    param_values[rate_mask] = 10 ** position
    traj = solver.run(param_values=param_values)

    # normalize trajectories
    bid_traj = traj.observables['tBid_total'] / model.parameters['Bid_0'].value
    cparp_traj = traj.observables['CPARP_total'] / model.parameters[
        'PARP_0'].value
    aSmac_traj = traj.observables['cSmac_total'] / model.parameters[
        'mSmac_0'].value

    # create all plots for each observable
    plt.figure(figsize=(3, 9))

    # plot cleaved parp
    plt.subplot(311)
    plt.plot(tspan, bid_traj, color='r', marker='^', label='tBID sim')
    plt.errorbar(exp_data['Time'], exp_data['norm_IC-RP'],
                 yerr=exp_data['nrm_var_IC-RP'] ** .5,
                 ecolor='black', color='black', elinewidth=0.5, capsize=0)
    plt.legend(loc=0)

    # plot cleaved parp
    plt.subplot(312)
    plt.plot(tspan, cparp_traj, color='blue', marker='*', label='cPARP sim')
    plt.errorbar(exp_data['Time'], exp_data['norm_EC-RP'],
                 yerr=exp_data['nrm_var_EC-RP'] ** .5,
                 ecolor='black', color='black', elinewidth=0.5, capsize=0)
    plt.legend(loc=0)

    # plot activated SMAC
    plt.subplot(313)
    plt.plot(tspan, aSmac_traj, color='g', label='aSMAC sim')
    plt.axvline(momp_data[0], -0.05, 1.05, color='black', linestyle=':',
                label='exp aSMAC')
    plt.legend(loc=0)
    plt.savefig('{}.png'.format(save_name), bbox_inches='tight')
    plt.close()


def create_gif_of_model_training(pso_instance):
    """
    Create a gif showing the improvements of parameter fits from a PSO instance
    """

    def create(position, save_name, values):
        param_values[rate_mask] = 10 ** position
        traj = solver.run(param_values=param_values)

        # normalize trajectories
        cparp_traj = traj.observables['CPARP_total'] / \
                     model.parameters['PARP_0'].value

        # create all plots for each observable
        plt.figure(figsize=(4, 6))
        time = tspan / 3600
        # plot cleaved parp
        plt.subplot(211)
        plt.plot(time, cparp_traj, color='r', label='CPARP sim')
        plt.errorbar(time, exp_data['norm_EC-RP'],
                     yerr=exp_data['nrm_var_EC-RP'] ** .5, label='CPARP exp',
                     ecolor='black', color='black', elinewidth=0.5, capsize=0)
        plt.legend(loc=2)
        plt.xlabel("Time(hr)")
        plt.subplot(212)
        plt.plot(range(len(values)), np.log10(values), '-k')
        plt.plot(range(len(values))[save_name], np.log10(values[save_name]),
                 'or')

        plt.xlabel("Iteration number")
        plt.ylabel("Objective function value")
        plt.subplots_adjust(hspace=.3)
        plt.savefig('images/{}.png'.format(save_name),
                    bbox_inches='tight', dpi=300)
        plt.close()

    try:
        import imageio
    except ImportError:
        raise ImportError("Please install imageio to create a gif")

    if not os.path.exists('images'):
        os.mkdir('images')
    for n, i in enumerate(pso_instance.history):
        create(i, n, pso_instance.values)

    images = ['images/{}.png'.format(i) for i in
              range(len(pso_instance.history) - 1)]

    imgs = [imageio.imread(i) for i in images]
    print("Creating gif")
    imageio.mimsave('training_earm.gif', imgs, duration=10 / len(images))
    print("Creating mp4")
    imageio.mimsave('training_earm.mp4', imgs, fps=len(images) / 10)


def run_example():
    # create PSO object
    pso = PSO(save_sampled=False, verbose=True, shrink_steps=False)
    pso.set_start_position(starting_position)

    # allows particles to move +/- 2 orders of magnitude
    pso.set_bounds(2)
    # sets maximum speed that a particle can travel
    pso.set_speed(-.1, .1)

    pso.run(
        num_particles=24,
        num_iterations=100,
        stop_threshold=1e-5,
        num_processors=18,
        max_iter_no_improv=20,
        cost_function=likelihood
    )

    display(pso.best.pos, save_name='best_fit')
    np.savetxt("pso_fit_for_model.csv", pso.best.pos)
    create_gif_of_model_training(pso)


if __name__ == '__main__':
    # Runs PSO
    run_example()
