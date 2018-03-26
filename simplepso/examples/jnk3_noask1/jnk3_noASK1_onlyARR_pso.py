from jnk3_no_ask1 import model
from simplepso.pso import PSO
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

exp_data = pd.read_csv('../data/exp_data_arrestin_noarrestin.csv')

idx_pars_calibrate = [21, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 37]
rates_of_interest_mask = [i in idx_pars_calibrate for i, par in enumerate(model.parameters)]

# Initial conditions of mkk4, mkk7, uJNK3 respectively
arrestin_idx = [38]

param_values = np.array([p.value for p in model.parameters])
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rates_of_interest_mask])
# ntimes = len(exp_data['Time (secs)'])
# tmul = 10
# tspan = np.linspace(exp_data['Time (secs)'].values[0], exp_data['Time (secs)'].values[-1],
#                     (ntimes-1) * tmul + 1)
tspan = exp_data['Time (secs)'].values
solver = ScipyOdeSimulator(model, tspan=tspan)


def display(position):
    Y = np.copy(position)
    param_values[rates_of_interest_mask] = 10 ** Y

    sim = solver.run(param_values=param_values).all

    plt.semilogx(exp_data['Time (secs)'].values, sim['pTyr_jnk3'], color='red', label='pJNK3 by MKK4 sim')
    plt.errorbar(exp_data['Time (secs)'].values, exp_data['pTyr_arrestin_avg'].values,
                 exp_data['pTyr_arrestin_std'].values,
                 linestyle='None', marker='o', capsize=5, color='red', label='pJNK3 by MKK4 exp')
    plt.semilogx(exp_data['Time (secs)'].values, sim['pThr_jnk3'], color='blue', label='pJNK3 by MKK7 sim')
    plt.errorbar(exp_data['Time (secs)'].values, exp_data['pThr_arrestin_avg'].values,
                 exp_data['pThr_arrestin_std'].values,
                 linestyle='None', marker='o', capsize=5, color='blue', label='pJNK3 by MKK7 exp')

    plt.xlabel('Arrestin (microM)')
    plt.ylabel('pJNK3 (microM)')
    plt.legend()
    plt.savefig('jnk3_noASK1_onlyARR_trained.png')
    plt.show()


def likelihood(position):
    Y = np.copy(position)
    param_values[rates_of_interest_mask] = 10 ** Y

    pars = np.copy(param_values)

    sim = solver.run(param_values=pars).all

    e_mkk4 = np.sum((exp_data['pTyr_arrestin_avg'].values - sim['pTyr_jnk3']) ** 2 /
                    (2 * exp_data['pTyr_arrestin_std'].values)) / len(exp_data['pTyr_arrestin_std'].values)
    e_mkk7 = np.sum((exp_data['pThr_arrestin_avg'].values - sim['pThr_jnk3']) ** 2 /
                    (2 * exp_data['pThr_arrestin_std'].values)) / len(exp_data['pThr_arrestin_std'].values)
    error1 = e_mkk4 + e_mkk7

    return error1,

# new_nominal = np.load('jnk3_noASK1_calibrated_pars_arr_noarr.npy')

def run_example():
    pso = PSO(save_sampled=False, verbose=True, num_proc=4)
    pso.set_cost_function(likelihood)
    pso.set_start_position(xnominal)
    pso.set_bounds(5)
    pso.set_speed(-.25, .25)
    pso.run(50, 300)
    display(pso.best)
    np.save('jnk3_noASK1_onlyARR_calibrated_pars', pso.best)

if __name__ == '__main__':
    run_example()