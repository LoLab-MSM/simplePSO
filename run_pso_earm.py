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
#from earm.lopez_embedded import model
#from earm.lopez_embedded import model
import sys


data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data.npy')

ydata_norm = numpy.load(data_filename)

exp_var = 0.2

rate_params = model.parameters_rules()

tspan = np.linspace(0,5.5 * 3600,len(ydata_norm))  # 5.5 hours, in seconds
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
    err = numpy.sum((ydata_norm - ysim_norm) ** 2 / (2 * exp_var ** 2))
    #err = (ydata_norm - ysim_norm) ** 2 / (2 * exp_var ** 2)
    #err= np.sum(err,axis=0)
    #print err
    return err,

def display(x):
    Y=np.copy(x)
    param_values = np.array([p.value for p in model.parameters])
    rate_mask = np.array([p in rate_params for p in model.parameters])
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    ysim_array = extract_records(solver.yobs, obs_names)
    ysim_norm = normalize(ysim_array)
    count=1
    #for j,obs_name in enumerate(obs_names):
      #plt.subplot(3,1,count)
      #plt.plot(solver.tspan,ysim_norm[:,j])
      #plt.plot(solver.tspan,ydata_norm[:,j],'-x')
      #plt.title(str(obs_name))
      #count+=1
    #plt.ylabel('concentration')
    #plt.xlabel('time (s)')
    #plt.show()

def generate(size, speedmin, speedmax):
    u1 = (random.uniform(-1.*bounds_radius, 1.*bounds_radius) for _ in range(size))
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
bounds_radius = 1
lb = xnominal - bounds_radius
ub = xnominal + bounds_radius


creator.create("FitnessMin", base.Fitness,weights=(-1.00,))
creator.create("Particle", list, fitness=creator.FitnessMin, \
    speed=list,smin=list, smax=list, best=None)
toolbox.register("particle", generate, size=np.shape(rate_params)[0],\
                 speedmin=-.5,speedmax=.5)
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

def init(sample,dictionary):
    global Sample
    global Dictionary
    Sample,Dictionary = sample,dictionary

def OBJ(block):
    #print block
    obj_values[block]=likelihood(sample[block])

if __name__ == '__main__':
    GEN = 500
    pop = toolbox.population(n=25)
    best = creator.Particle(xnominal)
    best.fitness.values = likelihood(xnominal)

    for g in range(0,GEN):
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
        worst = 0
        #for part in pop:
            #part.fitness.values = obj_values[count]
            #count+=1
            #if g == 0:
                #part.best = creator.Particle(part)
                #part.best.fitness.values = part.fitness.values

            #elif part.fitness.dominates2(part.best.fitness):
                #part.best = creator.Particle(part)
                #part.best.fitness.values = part.fitness.values

            #if part.fitness.dominates2(best.fitness):
                #%best = creator.Particle(part)
                #best.fitness.values = part.fitness.values
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
        for part in pop:
            toolbox.update(part, best)

        logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
        print(logbook.stream),best.fitness.values

    #quit()
    display(best)

    scores = []
    positions  = []
    for part in pop:
      scores.append(part.best.fitness.values)
      positions.append(part.best)

    positions = np.array(positions)
    #for i in range(len(positions)):
        #plt.title(model.parameters_rules()[i])
        #plt.hist(positions[:,i],bins=25,normed=1)
        #plt.vlines(np.log10(model.parameters_rules()[i].value), 0, 1)
        #plt.show()
    avg =np.average(positions,axis=0)
    std = np.std(positions,axis=0)
    #for i in range(0,np.shape(np.average(pop,axis=0))[0]):
    #    print model.parameters_rules()[i],avg[i],' ',std[i],'  ',xnominal[i],best[i]
    np.savetxt(str(sys.argv[1]),np.asarray(best))
#main()
