from jnk3_ASK1_released_v2 import model
from simplepso.pso import PSO
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

mkk4_data = np.loadtxt('data/mkk4_data.csv', delimiter=',')
sd_mkk4_data = np.loadtxt('data/sd_mkk4_data.csv', delimiter=',')
mkk7_data = np.loadtxt('data/mkk7_data.csv', delimiter=',')
sd_mkk7_data = np.loadtxt('data/sd_mkk7_data.csv', delimiter=',')

rates_of_interest_mask = [('kr' in par.name) or ('kcat' in par.name) or ('keq' in par.name) for par in model.parameters]

# rates_of_interest_mask = [False, False, False, False, False, False, False, False,
#                           True, True, True, True, True, True, True,
#                           False, False, False, False, False]

# initials_mask = [False, False, False, False, False, False, False, False, False, False, False,
#                  False, False, False, False, False, False, True, True, True, True, True]

# Initial conditions of mkk4, mkk7, uJNK3 respectively
initials_experiment_mkk4 = [0.05, 0, 0.5]
initials_experiment_mkk7 = [0, 0.05, 0.5]
initials_exps_idxs = [33, 34, 35]
arrestin_idx = [32]

param_values = np.array([p.value for p in model.parameters])
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rates_of_interest_mask])
tspan = np.linspace(0, 10, 50)
solver = ScipyOdeSimulator(model, tspan=tspan)


def display(position):
    Y = np.copy(position)
    param_values[rates_of_interest_mask] = 10 ** Y

    jnk3_mkk4_sim = np.zeros(13, dtype='float64')
    jnk3_mkk7_sim = np.zeros(13, dtype='float64')

    param_values[initials_exps_idxs] = initials_experiment_mkk4
    for i, arr_conc in enumerate(mkk4_data[:, 0]):
        param_values[arrestin_idx] = arr_conc
        sim = solver.run(param_values=param_values).all
        jnk3_mkk4_sim[i] = sim['mkk4_pjnk3'][-1]

    param_values[initials_exps_idxs] = initials_experiment_mkk7
    for i, arr_conc in enumerate(mkk7_data[:, 0]):
        param_values[arrestin_idx] = arr_conc
        sim = solver.run(param_values=param_values).all
        jnk3_mkk7_sim[i] = sim['mkk7_pjnk3'][-1]

    plt.plot(mkk4_data[:, 0], jnk3_mkk4_sim, 'x', color='red', label='pJNK3 by MKK4 sim')
    plt.errorbar(mkk4_data[:, 0], mkk4_data[:, 1], sd_mkk4_data[:, 1], linestyle='None', marker='o',
                 capsize=5, color='red', label='pJNK3 by MKK4 exp')
    plt.plot(mkk7_data[:, 0], jnk3_mkk7_sim, 'x', color='blue', label='pJNK3 by MKK7 sim')
    plt.errorbar(mkk7_data[:, 0], mkk7_data[:, 1], sd_mkk7_data[:, 1], linestyle='None', marker='o',
                 capsize=5, color='blue', label='pJNK3 by MKK7 exp')
    plt.xlabel('Arrestin (microM)')
    plt.ylabel('pJNK3 (microM)')
    plt.legend()
    plt.savefig('jnk3_ASK1_released_trained.png')
    plt.show()


def likelihood(position):
    Y = np.copy(position)
    param_values[rates_of_interest_mask] = 10 ** Y

    jnk3_mkk4_pars = [0] * 13
    jnk3_mkk7_pars = [0] * 13

    param_values[initials_exps_idxs] = initials_experiment_mkk4
    for i, arr_conc in enumerate(mkk4_data[:, 0]):
        param_values[arrestin_idx] = arr_conc
        jnk3_mkk4_pars[i] = np.copy(param_values)

    param_values[initials_exps_idxs] = initials_experiment_mkk7
    for i, arr_conc in enumerate(mkk7_data[:, 0]):
        param_values[arrestin_idx] = arr_conc
        jnk3_mkk7_pars[i] = np.copy(param_values)

    all_pars = jnk3_mkk4_pars + jnk3_mkk7_pars

    sims = solver.run(param_values=all_pars).all
    jnk3_mkk4_sim = [sim['mkk4_pjnk3'][-1] for sim in sims[:13]]
    jnk3_mkk7_sim = [sim['mkk7_pjnk3'][-1] for sim in sims[13:26]]

    e_mkk4 = np.sum((mkk4_data[:, 1] - jnk3_mkk4_sim) ** 2 / (2 * sd_mkk4_data[:, 1])) / len(sd_mkk4_data[:, 1])
    e_mkk7 = np.sum((mkk7_data[:, 1] - jnk3_mkk7_sim) ** 2 / (2 * sd_mkk7_data[:, 1])) / len(sd_mkk7_data[:, 1])
    error = e_mkk4 + e_mkk7
    return error,

new_nominal = np.load('jnk3_ASK1_released_calibrated_pars3.npy')

def run_example():
    pso = PSO(save_sampled=False, verbose=True, num_proc=4)
    pso.set_cost_function(likelihood)
    pso.set_start_position(new_nominal)
    pso.set_bounds(2)
    pso.set_speed(-.25, .25)
    pso.run(25, 2)
    display(pso.best)
    np.save('jnk3_ASK1_released_calibrated_pars', pso.best)

if __name__ == '__main__':
    run_example()