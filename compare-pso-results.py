# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 22:38:07 2015

@author: pinojc
"""
import pylab as plt
import numpy as np
import scipy.optimize
import scipy.interpolate
import numpy
import pysb.integrate
import pysb.util
import os
import inspect
import sys

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import PCA as mlabPCA
cm = plt.cm.get_cmap('RdYlBu')

#opt = sys.argv[1]
opt = 'Direct'

if opt == 'Direct':
    print opt
    directory = 'Direct'
    from earm.lopez_direct import model
    name = model.name
if opt == 'Indirect':
    print opt
    directory = 'Indirect'
    from earm.lopez_indirect import model
    name = model.name
if opt =='Embedded':
    print opt
    directory = 'Embedded'
    from earm.lopez_embedded import model
    name = model.name
directory = '/home/pinojc/Projects/ParticleSwarmOptimization/Single_objc/Direct'
data = np.loadtxt("%s/1_%s_overall_best.txt" % (directory,name))
for i in range(2,500):
    if os.path.exists("%s/%s_%s_overall_best.txt"% (directory,str(i),name)):
        tmp = np.loadtxt("%s/%s_%s_overall_best.txt"% (directory,str(i),name))
        data = np.column_stack((data,tmp))

#print np.shape(data)
#for i in range(0,105):
#    plt.hist(data[i,:])
#    plt.show()

obs_names = ['mBid', 'cPARP']
data_names = ['norm_ICRP', 'norm_ECRP']
var_names = ['nrm_var_ICRP', 'nrm_var_ECRP']
# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['Bid_0'].value,
              model.parameters['PARP_0'].value]

earm_path = '/home/pinojc/Projects/earm'
data_path = os.path.join(earm_path, 'xpdata', 'forfits',
                         'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
exp_data = np.genfromtxt(data_path, delimiter=',', names=True)

# Model observable corresponding to the IMS-RP reporter (MOMP timing)
momp_obs = 'aSmac'
# Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# yfinal (the last value of the IMS-RP trajectory)
momp_obs_total = model.parameters['Smac_0'].value
momp_data = np.array([9810.0, 180.0, momp_obs_total])
momp_var = np.array([7245000.0, 3600.0, 1e4])



ntimes = len(exp_data['Time'])
tmul = 10
tspan = np.linspace(exp_data['Time'][0], exp_data['Time'][-1],
                    (ntimes-1) * tmul + 1)
from pysb.bng import generate_network
#generate_network(model,verbose=True,cleanup=False)
#for rule in model.rules:
#    print rule.name
solver = pysb.integrate.Solver(model, tspan, integrator='lsoda', rtol=1e-6, atol=1e-6, nsteps=20000)
#solver = pysb.integrate.Solver(model, tspan, integrator='vode',  rtol=1e-5, atol=1e-5,)


rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
k_ids = [p.value for p in model.parameters_rules()]



def likelihood(position):
    Y=np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    for obs_name, data_name, var_name, obs_total in \
            zip(obs_names, data_names, var_names, obs_totals):
        ysim = solver.yobs[obs_name][::tmul]
        ysim_norm = ysim / obs_total
        ydata = exp_data[data_name]
        yvar = exp_data[var_name]
        if obs_name == 'mBid':
            e1 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
        else:
            e2 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)

    ysim_momp = solver.yobs[momp_obs]
    if np.nanmax(ysim_momp) == 0:
        print 'No aSmac!'
        ysim_momp_norm = ysim_momp
        t10 = 0
        t90 = 0
    else:
        ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
        st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
        try:
            t10 = scipy.interpolate.sproot((st, sc-0.10, sk))[0]
            t90 = scipy.interpolate.sproot((st, sc-0.90, sk))[0]
        except IndexError:
            t10 = 0
            t90 = 0

    td = (t10 + t90) / 2
    ts = t90 - t10
    yfinal = ysim_momp[-1]
    momp_sim = [td, ts, yfinal]
    e3 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3

    error = e1 + e2 +e3
    return error



scores = np.zeros((np.shape(data)[1],1))
for i in range(np.shape(data)[1]):
    scores[i]=likelihood(data[:,i])
#plt.hist(scores,50)
#plt.show()
d = np.histogram(scores)
print 'min of scores',min(scores)
print 'average= %s, std = %s, median = %s ' % (np.average(scores),np.std(scores),np.median(scores))
#scores = scores[scores[<np.average(scores)]
above_avg,throw = np.where( scores < d[1][1] )
copy_data = data[:,above_avg]
print np.shape(data)
print np.shape(copy_data)
copy_scores = scores[above_avg]
print np.shape(scores)
print np.shape(copy_scores)
#avg=np.mean(data,axis=1)
#[pcas,pcab] = np.linalg.eig(np.cov(data.T));
#si = np.argsort(-pcas.ravel());# print si;
#pcas = np.diag(pcas);
#pcab = pcab[si,:];
#
#pcacoffs = np.dot(pcab.conj().T, data.T-avg);
#pcacoffs =np.real(pcacoffs)
#
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#cNorm  = plt.matplotlib.colors.Normalize(vmin=np.min(scores), vmax=np.max(scores))
#sc = ax.scatter(pcacoffs[:,0],pcacoffs[:,1],pcacoffs[:,2],c=scores,cmap=cm,norm=cNorm)
#plt.colorbar(sc)
#plt.show()

#for i in range(len(k_ids)):
#    print np.std(data[i,:]),np.abs(np.average(data[i,:]))
#    plt.hist(data[i,:],50)
#    plt.show()


pca = mlabPCA(data.T)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cNorm  = plt.matplotlib.colors.Normalize(vmin=0, vmax=np.max(scores))
sc = ax.scatter(pca.Y[:,0],pca.Y[:,1],pca.Y[:,2],c=scores,cmap=cm,norm=cNorm)
plt.colorbar(sc)
plt.show()

pca = mlabPCA(copy_data.T)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cNorm  = plt.matplotlib.colors.Normalize(vmin=0, vmax=np.max(copy_scores))
sc = ax.scatter(pca.Y[:,0],pca.Y[:,1],pca.Y[:,2],c=copy_scores,cmap=cm,norm=cNorm)
plt.colorbar(sc)
plt.show()
