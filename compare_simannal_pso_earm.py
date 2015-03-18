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
#from earm.lopez_direct import model
#from earm.lopez_indirect import model
import sys
from pysb.util import load_params

best=np.loadtxt("/home/pinojc/Copy/HHMI/embedded.txt")
#for param in model.parameters:
#    if param.name in param_dict:
#        param.value = param_dict[param.name]

data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data.npy')

ydata_norm = numpy.load(data_filename)

exp_var = 0.2

rate_params = model.parameters_rules()

tspan = np.linspace(0,5.5 * 3600,len(ydata_norm)*10)  # 5.5 hours, in seconds
obs_names = ['mBid', 'aSmac', 'cPARP']
# Initialize solver object
#solver = pysb.integrate.Solver(model, tspan, integrator='lsoda', rtol=1e-6, atol=1e-6, nsteps=20000)
solver = pysb.integrate.Solver(model, tspan, integrator='vode',  with_jacobian=True,rtol=1e-5, atol=1e-5,)
# Get parameters for rates only
rate_params = model.parameters_rules()

param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
# Build a boolean mask for those params against the entire param list
k_ids = [p.value for p in model.parameters_rules()]


def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)

def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return numpy.vstack([recarray[name] for name in names]).T

def likelihood(x):

    Y=np.copy(x)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    err = numpy.sum((ydata_norm - ysim_norm[::10]) ** 2 / (2 * exp_var ** 2))
    #err = (ydata_norm - ysim_norm) ** 2 / (2 * exp_var ** 2)
    #err= np.sum(err,axis=0)
    #print err
    return err,

def objective_func(x):
    if np.any((x < lb) | (x > ub)):
        print "bounds-check failed"
        return np.inf
    Y=np.copy(x)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    err = numpy.sum((ydata_norm - ysim_norm[::10]) ** 2 / (2 * exp_var ** 2))
    
    return err



def display(x):
    Y=np.copy(x)
    param_values = np.array([p.value for p in model.parameters])
    rate_mask = np.array([p in rate_params for p in model.parameters])
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    count=1
    for j,obs_name in enumerate(obs_names):
      plt.subplot(3,1,count)
      plt.plot(solver.tspan,ysim_norm[:,j])
      plt.plot(solver.tspan[::10],ydata_norm[:,j],'-x')
      plt.title(str(obs_name))
      count+=1
    plt.ylabel('concentration')
    plt.xlabel('time (s)')
    plt.show()

def generate(size, speedmin, speedmax):
    u1 = (random.uniform(-1., 1.) for _ in range(size))
    tmp = map(operator.add,np.log10(k_ids),u1)
    part = creator.Particle(tmp)
    part.speed = [random.uniform(speedmin, speedmax) for _ in range(size)]
    part.smin = speedmin
    part.smax = speedmax
    return part

def updateParticle(part, best, phi1, phi2):

    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
    v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
    part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))
    for i, speed in enumerate(part.speed):
        if speed < part.smin:
            part.speed[i] = part.smin
        elif speed > part.smax:
            part.speed[i] =  part.smax
    part[:] = list(map(operator.add, part, part.speed))
    for i, pos in enumerate(part):
        if pos < lb[i]:
            part[i] = lb[i]
        elif pos > ub[i]:
            part[i] =  ub[i]

toolbox = base.Toolbox()
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 10
lb = xnominal - bounds_radius
ub = xnominal + bounds_radius


creator.create("FitnessMin", base.Fitness,weights=(-1.00,))
creator.create("Particle", list, fitness=creator.FitnessMin, \
    speed=list,smin=list, smax=list, best=None)
toolbox.register("particle", generate, size=np.shape(rate_params)[0],\
                 speedmin=-.2,speedmax=.2)
toolbox.register("population", tools.initRepeat, list, toolbox.particle)
toolbox.register("update", updateParticle, phi1=2, phi2=2)
toolbox.register("evaluate", likelihood)

stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", numpy.mean)
stats.register("std", numpy.std)
stats.register("min", numpy.min)
stats.register("max", numpy.max)
logbook = tools.Logbook()
logbook.header = ["gen", "evals"] + stats.fields

def estimate(start_values=None):

        # Set starting position to nominal parameter values if not specified
    if start_values is None:
        start_values = nominal_values
    else:
        assert start_values.shape == nominal_values.shape
    # Log-transform the starting position
    x0 = np.log10(start_values[rate_mask])
    # Displacement size for annealing moves
    dx = .02
    # The default 'fast' annealing schedule uses the 'lower' and 'upper'
    # arguments in a somewhat counterintuitive way. See
    # http://projects.scipy.org/scipy/ticket/1126 for more information. This is
    # how to get the search to start at x0 and use a displacement on the order
    # of dx (note that this will affect the T0 estimation which *does* expect
    # lower and upper to be the absolute expected bounds on x).
    lower = x0 - dx / 2
    upper = x0 + dx / 2
    # Log-transform the rate parameter values
    xnominal = np.log10(nominal_values[rate_mask])
    # Hard lower and upper bounds on x
    lb = xnominal - bounds_radius
    ub = xnominal + bounds_radius

    # Perform the annealing
    args = [rate_mask, lb, ub]
    (fopt, xmin, direc, iter, funcalls, warnflag) = \
        scipy.optimize.fmin_powell(objective_func, x0, full_output=True)
    # Construct vector with resulting parameter values (un-log-transformed)
    params_estimated = start_values.copy()
    params_estimated[rate_mask] = 10 ** xmin
    #display(fopt)
    # Display annealing results
    print 'fmin =%f' % xmin
#     for v in (fopt, xmin, direc, iter, funcalls, warnflag):
#         print v

    return params_estimated




def init(sample,dictionary):
    global Sample
    global Dictionary
    Sample,Dictionary = sample,dictionary

def OBJ(block):
    #print block
    obj_values[block]=likelihood(sample[block])

if __name__ == '__main__':
    GEN = 500
    num_particles = 50
    pop = toolbox.population(n=num_particles)
    best = creator.Particle(xnominal)
    best.fitness.values = likelihood(xnominal)
    best_values =[]
    evals = []
    for g in range(1,GEN+1):
        m = mp.Manager()
        obj_values = m.dict()
        sample = []
        for p in pop:
            sample.append(p)
        p = mp.Pool(4,initializer = init, initargs=(sample,obj_values))
        allblocks =range(len(pop))
        p.imap_unordered(OBJ,allblocks)
        p.close()
        p.join()
        count=0
        
        for part in pop:
            part.fitness.values = obj_values[count]
            count+=1
            if g == 1:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            elif part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if best.fitness < part.fitness:
                best = creator.Particle(part)
                best.fitness.values = part.fitness.values
        for part in pop:
            toolbox.update(part, best)
 
        logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
        print(logbook.stream),best.fitness.values
        best_values.append(best.fitness.values)
        evals.append(g*num_particles)
    
    #plt.plot(evals,best_values)
    #plt.show()
    #quit()
    #display(best)
    Y=np.copy(best)
    param_values = np.array([p.value for p in model.parameters])
    rate_mask = np.array([p in rate_params for p in model.parameters])
    param_values[rate_mask] = 10 ** Y
   
    params_estimated = estimate()
    params_estimated = estimate(param_values)





  
#    scores = []
#    positions  = []
#    for part in pop:
#      scores.append(part.best.fitness.values)
#      positions.append(part.best)
#
#    positions = np.array(positions)
#    #for i in range(len(positions)):
#        #plt.title(model.parameters_rules()[i])
#        #plt.hist(positions[:,i],bins=25,normed=1)
#        #plt.vlines(np.log10(model.parameters_rules()[i].value), 0, 1)
#        #plt.show()
#     avg =np.average(positions,axis=0)
#     std = np.std(positions,axis=0)
    #for i in range(0,np.shape(np.average(pop,axis=0))[0]):
    #    print model.parameters_rules()[i],avg[i],' ',std[i],'  ',xnominal[i],best[i]
    #np.savetxt('indirect_'+str(sys.argv[1]),np.asarray(best))  
    #np.savetxt('indirect.txt',np.asarray(best))
      

#    [pcas,pcab] = numpy.linalg.eig(np.cov(np.array(positions)));
#     si = np.argsort(-pcas.ravel()); print si;
#     pcas = numpy.diag(pcas);
#     pcab = pcab[:,si];
#     cm = plt.cm.get_cmap('RdYlBu')
#     pcacoffs = numpy.dot(pcab.conj().T, positions-avg);
#     pcacoffs =numpy.real(pcacoffs)
#     scores = np.array(scores)
#     sc = plt.scatter(pcacoffs[:,0],pcacoffs[:,1],c=scores,s=35,cmap=cm,)
#     plt.colorbar(sc)
#     plt.show()
#main()
