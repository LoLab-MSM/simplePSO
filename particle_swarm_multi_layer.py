# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 15:44:28 2014

@author: pinojc
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""
from varsens import *
import operator
import random
import scipy.optimize
import numpy
import pylab
from deap import base
import ghalton
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
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl

import multiprocessing as mp

    
from matplotlib.colors import LogNorm
        
cmap = plt.cm.jet

cNorm  = colors.Normalize(vmin=np.min(0), vmax=np.max(100))
scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)

def h12(individual):

    num = (sin(individual[0] - individual[1] / 8))**2 + (sin(individual[1] + individual[0] / 8))**2
    denum = ((individual[0] - 8.6998)**2 + (individual[1] - 6.7665)**2)**0.5 + 1
    return num / denum
def himmelblau2(individual):

    return (individual[0] * individual[0] + individual[1] - 11)**2 + \
        (individual[0] + individual[1] * individual[1] - 7)**2
def himmelblau(individual):
    """The Himmelblau's function is multimodal with 4 defined minimums in 
    :math:`[-6, 6]^2`.

    .. list-table:: 
       :widths: 10 50
       :stub-columns: 1

       * - Type
         - minimization
       * - Range
         - :math:`x_i \in [-6, 6]`
       * - Global optima
         - :math:`\mathbf{x}_1 = (3.0, 2.0)`, :math:`f(\mathbf{x}_1) = 0`\n
           :math:`\mathbf{x}_2 = (-2.805118, 3.131312)`, :math:`f(\mathbf{x}_2) = 0`\n
           :math:`\mathbf{x}_3 = (-3.779310, -3.283186)`, :math:`f(\mathbf{x}_3) = 0`\n
           :math:`\mathbf{x}_4 = (3.584428, -1.848126)`, :math:`f(\mathbf{x}_4) = 0`\n
       * - Function
         - :math:`f(x_1, x_2) = (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 -7)^2`

    .. plot:: code/benchmarks/himmelblau.py
        :width: 67 %
    """
    return (individual[0] * individual[0] + individual[1] - 11)**2 + \
        (individual[0] + individual[1] * individual[1] - 7)**2,

def h1(individual):
    #Y= np.copy(individual)
    #Y = np.log10(Y)
    """ Simple two-dimensional function containing several local maxima.
    From: The Merits of a Parallel Genetic Algorithm in Solving Hard 
    Optimization Problems, A. J. Knoek van Soest and L. J. R. Richard 
    Casius, J. Biomech. Eng. 125, 141 (2003)

    .. list-table:: 
       :widths: 10 50
       :stub-columns: 1

       * - Type
         - maximization
       * - Range
         - :math:`x_i \in [-100, 100]`
       * - Global optima
         - :math:`\mathbf{x} = (8.6998, 6.7665)`, :math:`f(\mathbf{x}) = 2`\n
       * - Function
         - :math:`f(\mathbf{x}) = \\frac{\sin(x_1 - \\frac{x_2}{8})^2 + \
            \\sin(x_2 + \\frac{x_1}{8})^2}{\\sqrt{(x_1 - 8.6998)^2 + \
            (x_2 - 6.7665)^2} + 1}`

    .. plot:: code/benchmarks/h1.py
       :width: 67 %
    """
    num = (sin(individual[0] - individual[1] / 8))**2 + (sin(individual[1] + individual[0] / 8))**2
    denum = ((individual[0] - 8.6998)**2 + (individual[1] - 6.7665)**2)**0.5 + 1
    #print num / denum
    return num / denum,
    
def generateSample(size,samples,discard,lb,ub):
    seq = ghalton.Halton(size)
    scaling =  lambda x:scale.linear(x, lower_bound=lb, upper_bound=ub)
    seq.get(20*size + int(discard)) # Remove initial linear correlated points plus any additional specified by the user.
    x = numpy.array(seq.get(samples))
    M_1 = scaling(x[0:samples,...])
    return M_1
    


def generate(size,  sv,range_percent, smin, smax):
    part = creator.Particle(sv) 
    tmp = [random.uniform(-1*range_percent, range_percent) for _ in range(size)]
    part[:] = list(map(operator.add, part, list(map(operator.mul, part, tmp))))
    #print part
    part.speed = [random.uniform(-1, 1) for _ in range(size)]
    part.smin = smin
    part.smax = smax

    return part

def updateParticle(part, best, lb, ub,phi1, phi2,):

    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
    v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
    part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))
    #print 'speed',part.speed
    for i, speed in enumerate(part.speed):
        #print i,'spI',part.speed
        #print 'max',part.smin[i],part.smax[i]
        if speed < part.smin[i]:
            #print 's'
            part.speed[i] = part.smin[i]
        elif speed > part.smax[i]:
            #print 'x'
            part.speed[i] = part.smax[i]
        #print i,'after',part.speed
    part[:] = list(map(operator.add, part, part.speed))
    for i, pos in enumerate(part):
        #if i ==1:
        #    print 'sp',part.speed[i],
        #    print 'cur',pos,
        #    print 'lb',lb[i],
        #    print 'ub',ub[i],
        if pos < lb[i]:
            #print '1'
            part[i] = lb[i]
        elif pos > ub[i]:
            #print '2',
            part[i] =  ub[i]
        #if i ==1:
        #    print 'new',part[i]
            
def setPos(part, newPos):
    part[:] = newPos


toolbox = base.Toolbox()

creator.create("FitnessMax", base.Fitness,weights=(-1.00,))
creator.create("Particle", list, fitness=creator.FitnessMax, \
    speed=list,smin=list, smax=list, best=None)

#toolbox.register("particle", generate, size=2, lb = lb,ub=lb, smin=-3, smax=3)
#toolbox.register("population", tools.initRepeat, list, toolbox.particle)
toolbox.register("update", updateParticle, phi1=1, phi2=1)
toolbox.register("setPos", setPos)
toolbox.register("evaluate", himmelblau)

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
    obj_values[block]=h1(sample[block])
def sampleRegion(start_values,realrange,range_percent):
    lb = start_values + ( realrange[0])
    ub = start_values + (realrange[1])
    minspeed = ub-lb
    toolbox.register("particle", generate, size=2, sv=start_values,range_percent=range_percent, smin=-.1*minspeed, smax=.1*minspeed)
    toolbox.register("population", tools.initRepeat, list, toolbox.particle)
    GEN = 10
    nParticles = 100
    pop = toolbox.population(n=nParticles)
    best = None
    for g in range(1,GEN):
        for part in pop:
            part.fitness.values = toolbox.evaluate(part)

            if not part.best  or part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if not best or best.fitness < part.fitness:
                best = creator.Particle(part)
                best.fitness.values = part.fitness.values
        for part in pop:
            toolbox.update(part, best,lb,ub)
    return best ,best.fitness.values, pop

if __name__ == '__main__':
    nominal_values = np.array([-10,10])
    nParticles = 50
    lower,upper=-6,6
    Y=generateSample(2,nParticles,20,lb=lower,ub=upper)
    data=[]
    score = []
    POP = []
    plt.figure(figsize=(8,8)) 
    print int((upper-lower)*-.1) ,int((upper-lower)*.1)
    tmpbest=Y[0]
    for i in xrange(0,len(Y)):
        #tmpbest=Y[0]
        #Y[i,1]*=-1
        tmpbest,tmpscore ,tmppop= sampleRegion(Y[i].reshape(2,1),[int((upper-lower)*-.1) ,int((upper-lower)*.1)],.1)
        data.append(tmpbest)
        score.append(tmpscore)
        POP.append(tmppop)
        #tmpbest =np.array(tmpbest)
        #Y[i]=tmpbest.reshape(2)
        #plt.figure()
        #print Y[i][0], Y[i][1], np.float(data[-1][0]), np.float(data[-1][1])
        dx,dy = np.float(data[-1][0])- Y[i][0]  , np.float(data[-1][1])-Y[i][1]
        colorVal = scalarMap.to_rgba(tmpscore)
        plt.arrow(Y[i,0],Y[i,1],dx,dy,)#color=colorVal)
        #plt.quiver( Y[i][0], Y[i][1], dx,dy,color=np.float(score[-1][0]) )
        
        #plt.plot(Y[i][0],Y[i][1],'o')
        #plt.plot(data[-1][0],data[-1][1],'.')
        #print data[-1]
    #plt.show()
    
    #plt.figure()
    #score = np.asarray(score)
    #print data
    data=np.asarray(data)
    dx = data[:,0] - Y[:,0]
    dy = data[:,1] - Y[:,1]
    #plt.arrow(Y[:,0],Y[:,1],dx,dy,color=np.array(score))
    plt.scatter(data[:,0],data[:,1],c=np.array(score),s=100)
    plt.ylim(lower,upper)
    plt.xlim(lower,upper)
    #plt.colorbar()
    #plt.show()
    #plt.gray()
    #plt.figure(figsize=(8,8))    
    blank = np.zeros((100,100))
    
    #for countx,i in enumerate(np.linspace(lower,upper,100)):
    #    for county,j in enumerate(np.linspace(lower,upper,100)):
    #        blank[countx,county]= himmelblau2([i,j])

    plt.colorbar()
    #plt.savefig('multiSwarm_pso_x6-0.png',dpi=200)
    #plt.figure(figsize=(8,8)) 
    #plt.imshow(blank,extent=[lower,upper,lower,upper],interpolation='none',origin='lower',)#norm=LogNorm(vmin=0.01,vmax=1))#vmin=np.min(blank), vmax=np.max(blank))
    #plt.colorbar()    
    #plt.savefig('landscapex6-0.png',dpi=200)
    