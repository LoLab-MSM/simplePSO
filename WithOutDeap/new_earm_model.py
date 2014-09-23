# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 12:10:47 2014

@author: pinojc
"""

#!/usr/bin/python

import numpy as np
import pylab as plt
import scipy.interpolate
import pysb.integrate
import pysb.util
import os
from earm.lopez_embedded import model
import PSO_CLASS_multiprocc


# List of model observables and corresponding data file columns for
# point-by-point fitting
obs_names = ['mBid', 'cPARP']
data_names = ['norm_ICRP', 'norm_ECRP']
var_names = ['nrm_var_ICRP', 'nrm_var_ECRP']
# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['Bid_0'].value, model.parameters['PARP_0'].value]

# Load experimental data file
#earm_path = os.path.dirname(__file__)
#data_path = os.path.join(earm_path, 'xpdata', 'forfits',
#                         'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
exp_data = np.genfromtxt('/home/pinojc/Projects/earm/xpdata/forfits/\
EC-RP_IMS-RP_IC-RP_data_for_models.csv', delimiter=',', names=True)


# Build time points for the integrator, using the same time scale as the
# experimental data but with greater resolution to help the integrator converge.
ntimes = len(exp_data['Time'])
# Factor by which to increase time resolution
tmul = 10
# Do the sampling such that the original experimental timepoints can be
# extracted with a slice expression instead of requiring interpolation.
tspan = np.linspace(exp_data['Time'][0], exp_data['Time'][-1],(ntimes-1) * tmul + 1)
# Initialize solver object
solver = pysb.integrate.Solver(model, tspan,rtol=1e-5, atol=1e-5,\
              integrator='lsoda')
rate_params = model.parameters_rules()

momp_obs = 'aSmac'
# Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# yfinal (the last value of the IMS-RP trajectory)
momp_obs_total = model.parameters['Smac_0'].value
momp_data = np.array([9810.0, 180.0, momp_obs_total])
momp_var = np.array([7245000.0, 3600.0, 1e4])

# Build time points for the integrator, using the same time scale as the
# experimental data but with greater resolution to help the integrator converge.
ntimes = len(exp_data['Time'])
# Factor by which to increase time resolution
tmul = 10
# Do the sampling such that the original experimental timepoints can be
# extracted with a slice expression instead of requiring interpolation.
tspan = np.linspace(exp_data['Time'][0], exp_data['Time'][-1],
                    (ntimes-1) * tmul + 1)
# Initialize solver object
solver = pysb.integrate.Solver(model, tspan,rtol=1e-5, atol=1e-5,\
              integrator='lsoda')

# Get parameters for rates only
rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
# Build a boolean mask for those params against the entire param list


def objective_func(x):
    # Simulate model with rates taken from x
    #create a copy to allow to be log transformed
    Y=np.copy(x)
    #log transform X
    #print param_values
    param_values[rate_mask] = 10 ** Y
    #print param_values
    solver.run(param_values)
    # Calculate error for point-by-point trajectory comparisons
    e1 = 0
    
    for obs_name, data_name, var_name, obs_total in \
            zip(obs_names, data_names, var_names, obs_totals):
        # Get model observable trajectory (this is the slice expression
        # mentioned above in the comment for tspan)
        ysim = solver.yobs[obs_name][::tmul]
        # Normalize it to 0-1
        ysim_norm = ysim / np.nanmax(ysim)
        # Get experimental measurement and variance
        ydata = exp_data[data_name]
        yvar = exp_data[var_name]
        # Compute error between simulation and experiment (chi-squared)
        e1 += np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
        #e1 += np.sum((ydata - ysim_norm) ** 2)
    ysim_momp = solver.yobs[momp_obs]
    ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
    # Build a spline to interpolate it
    st, sc, sk = scipy.interpolate.splrep(solver.tspan, ysim_momp_norm)
    if np.any(np.isnan(st)) or np.any(np.isnan(sc)) or np.any(np.isnan(sk)):
      e2=1000000.
    else:
      #print st,sc,sk
      # Use root-finding to find the point where trajectory reaches 10% and 90%
      t10 = scipy.interpolate.sproot((st, sc-0.10, sk))[0]
      t90 = scipy.interpolate.sproot((st, sc-0.90, sk))[0]
      # Calculate Td as the mean of these times
      td = (t10 + t90) / 2
      # Calculate Ts as their difference
      ts = t90 - t10
      # Get yfinal, the last element from the trajectory
      yfinal = ysim_momp[-1]
      # Build a vector of the 3 variables to fit
      momp_sim = [td, ts, yfinal]
      #print momp_sim
      # Perform chi-squared calculation against mean and variance vectors
      e2 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3 
    if np.isnan(e1):
      e1 = 1000000.
    #print e1,e2
    error = e1 + e2
    #print error
    return error


        
def display(x):
    Y=np.copy(x)
    param_values = np.array([p.value for p in model.parameters])
    rate_mask = np.array([p in rate_params for p in model.parameters])
    #log transform X
    #print param_values
    param_values[rate_mask] = 10 ** Y
    # Construct matrix of experimental data and variance columns of interest


    # Simulate model with new parameters and construct a matrix of the
    # trajectories of the observables of interest, normalized to 0-1.
    solver.run(param_values)
    obs_names_disp = obs_names + ['aSmac']
    sim_obs = solver.yobs[obs_names_disp].view(float).reshape(len(solver.yobs), -1)
    totals = obs_totals + [momp_obs_total]
    sim_obs_norm = (sim_obs / totals).T
    count=1
    for obs_name, data_name, var_name, obs_total in \
            zip(obs_names, data_names, var_names, obs_totals):
    # Get model observable trajectory (this is the slice expression
    # mentioned above in the comment for tspan)
      ysim = solver.yobs[obs_name][::tmul]
      # Normalize it to 0-1
      ysim_norm = ysim / np.nanmax(ysim)
      # Get experimental measurement and variance
      ydata = exp_data[data_name]
      yvar = exp_data[var_name]
      plt.subplot(3,1,count)
      count+=1
      plt.plot(solver.tspan[::tmul],ysim_norm)
      plt.errorbar(solver.tspan[::tmul],ydata,xerr=yvar)  
      plt.ylabel('concentration')
      plt.title(str(obs_name))
    plt.subplot(3,1,count)      
    plt.plot(solver.tspan,sim_obs_norm[2],color='g')
    plt.vlines(momp_data[0], -0.05, 1.05, color='g', linestyle=':',label='aSmac')
    plt.xlabel('time (s)')
    plt.show()
    plt.clf()

    return param_values
        
"""
def __init__(self, swarm, best = None, position = None, 
               velocity = None, fitness = None,obj_function = None,
               max_vel = None, min_vel = None,max_pos = None, min_pos = None,
              dim = None):      
"""

DIM = 106
nominal_values = np.array([p.value for p in model.parameters])
x0 = np.log10(nominal_values[rate_mask])
f,MIN,MAX,SAVE = objective_func,-10,3,'EARM'
best,fitness,pop=PSO_CLASS.PSO(swarm_size=50,iterations=500,pos=x0,best=x0,max_pos=MAX,
        min_pos=MIN,max_vel=2.,min_vel=-2.,obj_function=f,dim=DIM)



end=display(best)
#display(x0)
pop = np.asarray(pop)
for i in xrange(1,10):
    plt.plot(best[i],best[i+1],'o')
    plt.plot(x0[i],x0[i+1],'x')
    plt.plot(pop[:,i],pop[:,i+1],'.')
    plt.xlim(-10,3)
    plt.ylim(-10,3)
    plt.show()
