import matplotlib
print matplotlib.rcParams['backend']
matplotlib.use('Qt4Agg')

import matplotlib.pyplot as plt
plt.interactive(False)
import numpy as np

from pysb.examples.robertson import model
from pysb.integrate import odesolve
import pysb
from simplepso.pso import PSO

# Next we will import the model and the odesolver from pysb
t = np.linspace(0, 40,100)
obs_names = ['A_total', 'C_total']

solver = pysb.integrate.Solver(model, t, integrator='vode',rtol=1e-8, atol=1e-8)
solver.run()

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

plt.plot(t,norm_data)
plt.legend(['A_Total','C_Total'], loc = 0)

noisy_data_A = ysim_array[:,0] + np.random.uniform(-0.02,0.02,np.shape(ysim_array[:,0]))
norm_noisy_data_A = normalize(noisy_data_A)
noisy_data_C = ysim_array[:,1] + np.random.uniform(-.02,.01,np.shape(ysim_array[:,1]))
norm_noisy_data_C = normalize(noisy_data_C)
ydata_norm = np.column_stack((norm_noisy_data_A,norm_noisy_data_C))

plt.plot(t,norm_noisy_data_A)
plt.plot(t,norm_noisy_data_C)
plt.plot(t,norm_data)
plt.legend(['A_total_noisy','C_total_noisy','A_total', 'B_total', 'C_total'], loc=0)


rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])

nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])


def display(x=None):
    if x is None:
        solver.run(param_values)
    else:
        Y=np.copy(x)
        param_values[rate_mask] = 10 ** Y
        solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    plt.figure(figsize=(8,6),dpi=200)
    plt.plot(t,ysim_norm[:,0],label='A')
    plt.plot(t,ysim_norm[:,1],label='C')
    plt.plot(t,norm_noisy_data_A,label='Noisy A')
    plt.plot(t,norm_noisy_data_C,label='Noisy C')
    plt.legend(loc=0)
    plt.ylabel('concentration')
    plt.xlabel('time (s)')
    plt.show()

def obj_function(params):
    params_tmp = np.copy(params)
    param_values[rate_mask] = 10 ** params_tmp
    solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    err = np.sum((ydata_norm - ysim_norm) ** 2 )
    return err,

print xnominal**10
start_position = xnominal +np.random.uniform(-1.0,1.0,size = np.shape(xnominal))
print start_position**10

optimizer = PSO(cost_function=obj_function,start=start_position,verbose=True)
optimizer.set_bounds(parameter_range=2)
optimizer.run(100,100)

display(start_position)
display(optimizer.best)