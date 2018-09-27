# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate

from earm_lopez_embedded_flat import model
from pysb.simulator import ScipyOdeSimulator
from simplepso.pso import PSO

# load experimental data
data_path = os.path.join(os.path.dirname(__file__), 'data',
                         'EC-RP_IMS-RP_IC-RP_data_for_models.csv')

exp_data = pd.read_csv(data_path, index_col=False)

# timepoints for simulation. These must match the experimental data.
tspan = exp_data['Time'].values.copy()

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
starting_position = np.log10(param_values[rate_mask])
args = {'integrator': 'lsoda', 'use_analytic_jacobian': True,
        "rtol": 1e-6, "atol": 1e-6}

solver = ScipyOdeSimulator(model, tspan, **args)

# Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# yfinal (the last value of the IMS-RP trajectory)
momp_data = np.array([9810.0, 180.0, model.parameters['Smac_0'].value])
momp_var = np.array([7245000.0, 3600.0, 1e4])


def display(position, save_name):
    param_values[rate_mask] = 10 ** position
    traj = solver.run(param_values=param_values)

    # normalize trajectories
    bid_traj = traj.observables['mBid'] / model.parameters['Bid_0'].value
    cparp_traj = traj.observables['cPARP'] / model.parameters['PARP_0'].value
    aSmac_traj = traj.observables['aSmac'] / model.parameters['Smac_0'].value

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
    plt.savefig('{}.png'.format(save_name))
    plt.close()


def likelihood(position):
    param_values[rate_mask] = 10 ** position
    traj = solver.run(param_values=param_values)

    # normalize trajectories
    bid_traj = traj.observables['mBid'] / model.parameters['Bid_0'].value
    cparp_traj = traj.observables['cPARP'] / model.parameters['PARP_0'].value
    momp_traj = traj.observables['aSmac']

    # calculate chi^2 distance for each time course
    e1 = np.sum((exp_data['norm_IC-RP'] - bid_traj) ** 2 /
                (2 * exp_data['nrm_var_IC-RP'])) / len(bid_traj)

    e2 = np.sum((exp_data['norm_EC-RP'] - cparp_traj) ** 2 /
                (2 * exp_data['nrm_var_EC-RP'])) / len(cparp_traj)

    # Here we fit a spline to find where we get 50% release of MOMP reporter
    if np.nanmax(momp_traj) == 0:
        print('No aSmac!')
        t10 = 0
        t90 = 0
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
    yfinal = momp_traj[-1]
    momp_sim = [td, ts, yfinal]

    e3 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3
    # return sum of errors ( the ',' is required)
    return e1 + e2 + e3,


def run_example():
    # Runs the cost function to calculate error between model and data
    print("Error at start = {}".format(likelihood(starting_position)[0]))
    # Displays the model with defaul positions
    display(starting_position, save_name='starting_position')

    # create PSO object
    pso = PSO(save_sampled=False, verbose=True, num_proc=4)
    pso.set_cost_function(likelihood)
    pso.set_start_position(starting_position)
    # allows particles to move +/- 2 orders of magnitude
    pso.set_bounds(2)
    # sets maximum speed that a particle can travel
    pso.set_speed(-.25, .25)

    pso.run(num_particles=25, num_iterations=50, stop_threshold=1e-5)
    display(pso.best, save_name='best_fit')
    np.savetxt("pso_fit_for_model.csv", pso.best)


if __name__ == '__main__':
    # Runs PSO
    run_example()
