# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 22:38:07 2015

@author: pinojc
"""
import pylab as plt
import numpy as np

import scipy.optimize
import numpy
import pysb.integrate
import pysb.util
import os
import inspect
import multiprocessing
import multiprocessing as mp
from earm.lopez_indirect import model
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
    return np.float(err)




data = np.loadtxt("indirect_1.txt")
for i in range(2,12):
    tmp = np.loadtxt("indirect_"+str(i)+'.txt')
    data = np.column_stack((data,tmp))
#for i in range(27,250):
#    tmp = np.loadtxt(str(i)+'.txt')
#    data = np.column_stack((data,tmp))
#print np.shape(data)


scores = np.zeros((np.shape(data)[1],1))
for i in range(np.shape(data)[1]):
    scores[i]=likelihood(data[:,i])
    
avg=np.mean(data,axis=1)
[pcas,pcab] = np.linalg.eig(np.cov(data.T));
si = np.argsort(-pcas.ravel());# print si;
pcas = np.diag(pcas);
pcab = pcab[si,:];
cm = plt.cm.get_cmap('RdYlBu')
pcacoffs = np.dot(pcab.conj().T, data.T-avg);
pcacoffs =np.real(pcacoffs)
#
from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#cNorm  = plt.matplotlib.colors.Normalize(vmin=0, vmax=np.max(scores))
#sc = ax.scatter(pcacoffs[:,0],pcacoffs[:,1],pcacoffs[:,2],c=scores,cmap=cm,norm=cNorm)
#plt.colorbar(sc)
#plt.show()

#for i in range(105):
#    tmp= np.std(data[i,:])/np.abs(np.average(data[i,:]))*100
#    if tmp > 100.:
#        print np.std(data[i,:]),np.abs(np.average(data[i,:]))
#        plt.hist(data[i,:],50)
#        plt.show()

from matplotlib.mlab import PCA as mlabPCA
pca = mlabPCA(data.T)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cNorm  = plt.matplotlib.colors.Normalize(vmin=0, vmax=np.max(scores))
sc = ax.scatter(pca.Y[:,0],pca.Y[:,1],pca.Y[:,2],c=scores,cmap=cm,norm=cNorm)
plt.colorbar(sc)
plt.show()
