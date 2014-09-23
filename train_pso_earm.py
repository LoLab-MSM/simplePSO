# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""

import operator
import random
import scipy.optimize
import numpy
import pylab
from deap import base

from deap import creator
from deap import tools
import pysb.integrate
import pysb.util
import numpy as np
import scipy.optimize
import scipy.interpolate
import matplotlib.pyplot as plt
import os
import inspect
import multiprocessing
import multiprocessing as mp
from earm.lopez_embedded import model


obs_names = ['mBid', 'cPARP']
data_names = ['norm_ICRP', 'norm_ECRP']
var_names = ['nrm_var_ICRP', 'nrm_var_ECRP']
# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['Bid_0'].value, model.parameters['PARP_0'].value]

# Load experimental data file
earm_path = os.path.dirname(__file__)
data_path = os.path.join(earm_path, 'xpdata', 'forfits',
                         'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
exp_data = np.genfromtxt('/home/pinojc/Projects/earm/xpdata/forfits/\
EC-RP_IMS-RP_IC-RP_data_for_models.csv', delimiter=',', names=True)


rate_params = model.parameters_rules()

momp_obs = 'aSmac'
# Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# yfinal (the last value of the IMS-RP trajectory)
momp_obs_total = model.parameters['Smac_0'].value
momp_data = np.array([9810.0, 180.0, momp_obs_total])
momp_var = np.array([7245000.0, 3600.0, 1e4])


tspan = exp_data['Time']


# Initialize solver object
#solver = pysb.integrate.Solver(model, tspan, integrator='lsoda', rtol=1e-6, atol=1e-6, nsteps=20000)
solver = pysb.integrate.Solver(model, tspan, integrator='vode',  with_jacobian=True,rtol=1e-5, atol=1e-5, nsteps=20000)
# Get parameters for rates only
rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
# Build a boolean mask for those params against the entire param list
k_ids = [p.value for p in model.parameters_rules()]
solver.verbose = False     

    

        
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
      ysim = solver.yobs[obs_name]
      # Normalize it to 0-1
      ysim_norm = ysim / np.nanmax(ysim)
      # Get experimental measurement and variance
      ydata = exp_data[data_name]
      yvar = exp_data[var_name]
      plt.subplot(3,1,count)
      count+=1
      plt.plot(solver.tspan,ysim_norm)
      plt.errorbar(solver.tspan,ydata,xerr=yvar)  
      plt.ylabel('concentration')
      plt.title(str(obs_name))
    plt.subplot(3,1,count)      
    plt.plot(solver.tspan,sim_obs_norm[2],color='g')
    plt.vlines(momp_data[0], -0.05, 1.05, color='g', linestyle=':',label='aSmac')
    plt.xlabel('time (s)')
    plt.show()
    plt.clf()

def plot_all(pop):
    for part in pop:
          Y=np.copy(part)
          param_values = np.array([p.value for p in model.parameters])
          rate_mask = np.array([p in rate_params for p in model.parameters])
          #log transform X
          #print param_values
          param_values[rate_mask] = 10 ** Y
          # Construct matrix of experimental data and variance columns of interest
          # Simulate model with new parameters and construct a matrix of the
          # trajectories of the observables of interest, normalized to 0-1.
          solver.run(param_values)
          ysim = solver.yobs['mBid']
          plt.subplot(3,1,1)
          ysim_norm = ysim / np.nanmax(ysim) 
          plt.plot(solver.tspan,ysim_norm)
          plt.subplot(3,1,2)
          ysim1 = solver.yobs['cPARP']
          ysim1_norm = ysim1 / np.nanmax(ysim1)
          plt.plot(solver.tspan,ysim1_norm)
          plt.subplot(3,1,3)
          ysim2 = solver.yobs['aSmac']
          ysim2_norm = ysim2 / np.nanmax(ysim2) 
          plt.plot(solver.tspan,ysim2_norm)
    plt.show()
def objective_func(x):
  
    # Simulate model with rates taken from x
    #create a copy to allow to be log transformed
    Y=np.copy(x)
    #if np.any((x < lb) | (x > ub)):
        #print "bounds-check failed"
    #    return 10000,
    #log transform X
    #print param_values
    param_values[rate_mask] = 10 ** Y
    #print param_values
    solver.run(param_values)
    # Calculate error for point-by-point trajectory comparisons
    e1 = 0
    
    ysim = solver.yobs['mBid']
    ysim_norm = ysim / np.nanmax(ysim)
    ydata = exp_data['norm_ICRP']
    yvar = exp_data['nrm_var_ICRP']
    e1 += np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
    ysim2 = solver.yobs['cPARP']
    ysim_norm2 = ysim2 / np.nanmax(ysim2)
    ydata2 = exp_data['norm_ECRP']
    yvar2 = exp_data['nrm_var_ECRP']
    e1 += np.sum((ydata2 - ysim_norm2) ** 2 / (2 * yvar2)) / len(ydata2)
    ysim_momp = solver.yobs['aSmac']
    ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
    # Build a spline to interpolate it
    st, sc, sk = scipy.interpolate.splrep(solver.tspan, ysim_momp_norm)
    if np.any(np.isnan(st)) or np.any(np.isnan(sc)) or np.any(np.isnan(sk)):
      e2=100000.
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
      #print e2
    if np.isnan(e1):
      e1 = 100000.
    error = e1 + e2
    return error,


def generate(size, pmin, pmax, speedmin, speedmax):
    u1 = (random.uniform(-1.*bounds_radius, 1.*bounds_radius) for _ in range(size)) 
    #tmp = (random.uniform(pmin,pmax) for _ in range(size))
    tmp = map(operator.add,np.log10(k_ids),u1)
    part = creator.Particle(tmp)    
    part.speed = [random.uniform(10*speedmin, 10*speedmax) for _ in range(size)]
    part.smin = [random.uniform(0, speedmax) for _ in range(size)]
    part.smax = [random.uniform(speedmax, 0) for _ in range(size)]
    for i, pos in enumerate(part):    
      part.smin[i] += np.abs(xnominal[i])*speedmin
      part.smax[i] += np.abs(xnominal[i])*speedmax
    return part

def updateParticle(part, best, phi1, phi2):

    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
    v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
    part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))
    for i, speed in enumerate(part.speed):
        if speed < part.smin[i]:
            part.speed[i] = part.smin[i]
        elif speed > part.smax[i]:
            part.speed[i] =  part.smax[i]
    part[:] = list(map(operator.add, part, part.speed))
    for i, pos in enumerate(part):
        if pos < lb[i]:
            part[i] = lb[i]
        elif pos > ub[i]:
            part[i] =  ub[i]



toolbox = base.Toolbox()
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 2
lb = xnominal - bounds_radius
ub = xnominal + bounds_radius


creator.create("FitnessMin", base.Fitness,weights=(-1.00,))
creator.create("Particle", list, fitness=creator.FitnessMin, \
    speed=list,smin=list, smax=list, best=None)


toolbox.register("particle", generate, size=np.shape(rate_params)[0], \
      pmin=-9, pmax=5,speedmin=-.5,speedmax=.5)
toolbox.register("population", tools.initRepeat, list, toolbox.particle)
toolbox.register("update", updateParticle, phi1=1, phi2=1)
toolbox.register("evaluate", objective_func)

stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", numpy.mean)
stats.register("std", numpy.std)
stats.register("min", numpy.min)
stats.register("max", numpy.max)
logbook = tools.Logbook()
logbook.header = ["gen", "evals"] + stats.fields

def init(sample,dictionary):
    global Sample
    global Dictionary
    Sample,Dictionary = sample,dictionary
	
def OBJ(block):
    #print block
    obj_values[block]=objective_func(sample[block])
    
if __name__ == '__main__':
    GEN = 100
    pop = toolbox.population(n=100)
    best = None
    for g in range(1,GEN):
        m = mp.Manager()
        obj_values = m.dict()
        sample = []
        for p in pop:
            sample.append(p)
        p = mp.Pool(8,initializer = init, initargs=(sample,obj_values))
        allblocks =range(len(pop))    
        p.imap_unordered(OBJ,allblocks)
        p.close()
        p.join()
        count=0
        worst = 0
        for part in pop:
            part.fitness.values = obj_values[count]
            count+=1
            if not part.best  or part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if not best or best.fitness < part.fitness:
                best = creator.Particle(part)
                best.fitness.values = part.fitness.values
        for part in pop:
            toolbox.update(part, best)


          
        

        logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
        print(logbook.stream),best.fitness.values
    
    #print logbook
    print best.fitness.values
    avg =np.average(pop,axis=0)
    std = np.std(pop,axis=0)
    for i in range(0,np.shape(np.average(pop,axis=0))[0]):
        #print avg[i]," ",std[i]
        print avg[i],' ',std[i],'  ',xnominal[i]
    plt.errorbar(np.arange(0,len(avg)),avg,yerr=std)
    plt.show()
    
    plt.hold(True)
    plt.clf()
    empty=[]
    for i in pop:
      empty.append(i.fitness.values)
    empty=np.asarray(empty)
    plt.hist(empty,25)
    plt.show()
    plt.clf()
    #plot_all(pop)
    #POP = np.array(pop) - avg
    #[pcas,pcab] = numpy.linalg.eig(np.cov(POP));
    #si = np.argsort(-pcas.ravel()); print si;
    #pcas = numpy.diag(pcas);
    #pcab = pcab[:,si];

    #fig = plt.figure();
    #ax = fig.add_subplot(111);
    #plt.plot(pcas[0,:],pcas[1,:])
    scores = []
    positions  = []
    for part in pop:
      scores.append(part.best.fitness.values)
      positions.append(part.best)
    #plot_all(positions)
    Data = np.column_stack((scores,positions))
    LLHOOD = np.array(sorted(Data, key=lambda score_entry: score_entry[0]))
    plt.plot(LLHOOD[:10,0])    
    plt.show()
    #plot_all(LLHOOD[0:2,1:])
    for i in range(0,10):
        display(LLHOOD[i,1:])
    plt.plot(LLHOOD[:10,1:].T)
    for i in range(1,len(LLHOOD.T)):
        print np.average(LLHOOD[:10,i]),np.std(LLHOOD[:10,i])
    plt.plot(np.std(LLHOOD[:10,:].T,axis=1))
    #scores = np.array(np.log10(scores))
    #pcacoffs = numpy.dot(pcab.conj().T,POP);
    #pcacoffs =numpy.real(pcacoffs)
    #p= ax.scatter(pcacoffs[:,0],pcacoffs[:,1],c=scores)
    #fig.colorbar(p)
    #plt.show()
