from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
from necroptosismodule import model
from tropical.dynamic_signatures_range import run_tropical_multi
# from tropical.examples.double_enzymatic.mm_two_paths_model import model
# from tropical import clustering
from necro_pso_kd_plot import display
from operator import itemgetter

# model.enable_synth_deg()

pos = np.load('position_pso.npy')


pars =  np.load('optimizer_best_5000_all_new_2.npy')
print(pars)
pars = pars.remove(pars[26])
print(len(pars))
print(pars)
# print(pars)
quit()
# ex = np.load('examples.npy')
# # print(ex.shape)
# print(ex[:,:])
# quit()

print(len(model.parameters))
print(len(pars))
print(len(model.parameters_rules()))
print(len(model.initial_conditions))
quit()

initials1 = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000]
params = 10 ** pars
print(len(model.parameters_rules()))
print(params)
# print(len(params))
params1 = np.concatenate([[1.0], params])
print(params1)
quit()


params = np.concatenate([initials1, params])
# print(params)

#
# print(len(params))
# print(len(model.parameters))
# quit()

tspan = np.linspace(0, 480, 481)
# sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=params).run()
#
#
# tro = run_tropical_multi(model, simulations=sim, passengers_by='imp_nodes', diff_par=1.0, cpu_cores=25)
#
#
# with open('necro_signatures.pickle', 'wb') as handle:
#     pickle.dump(tro, handle, protocol=pickle.HIGHEST_PROTOCOL)
#
# signatures_23 = tro[23]['production']
# clus = clustering.ClusterSequences(seqdata=signatures_23, unique_sequences=False)
#
# #property of the sequences of the class
# clus.diss_matrix(n_jobs=25)
# sil_df = clus.silhouette_score_spectral_range(cluster_range=range(2, 20), n_jobs=cpus, random_state=1234)
# print (sil_df) #dataframe w/ # of clusters and sil score
# best_silh, n_clus = sil_df.loc[sil_df['cluster_silhouette'].idxmax()]
# n_clus = int(n_clus)
#
#
# clus.spectral_clustering(n_clusters=n_clus, n_jobs=cpus, random_state=1234)
# np.save('labels_rip1.npy', clus.labels)
# b = clustering.PlotSequences(clus)
# b.modal_plot(title='RIP1')
# b.all_trajectories_plot(title='RIP1')


#
#
#
#
sim2 = ScipyOdeSimulator(model, tspan=tspan)
L4 = sim2.run(param_values=params)
a = [0, 1*60, 2*60, 3*60, 4*60, 5*60, 6*60, 7*60, 8*60]
print(L4.observables['MLKLa_obs'][a[1]])
print(L4.observables['MLKLa_obs'][a[2]])
print(L4.observables['MLKLa_obs'][a[3]])
print(L4.observables['MLKLa_obs'][a[4]])
print(L4.observables['MLKLa_obs'][a[5]])
print(L4.observables['MLKLa_obs'][a[6]])
print(L4.observables['MLKLa_obs'][a[7]])
print(L4.observables['MLKLa_obs'][a[8]])
# L3 = sim2.run()
#
#
plt.figure()
# # plt.plot(tspan/60, L3.observables['MLKLa_obs'],color = 'purple',label = 'MLKLp')
# plt.plot(tspan/60, sim_result1.observables['MLKLa_obs'],label = '100 ng/ml TNF')
# plt.plot(tspan/60, sim_result2.observables['MLKLa_obs'],label = '10 ng/ml TNF')
# plt.plot(tspan/60, sim_result3.observables['MLKLa_obs'],label = '1 ng/ml TNF')
# plt.plot(tspan/60, sim_result4.observables['MLKLa_obs'],label = '0.1 ng/ml TNF')
# # plt.plot(tspan/60, L3.observables['MLKL_obs'],color = 'red',label = 'MLKL')
plt.plot(tspan/60, L4.observables['MLKLa_obs'],color = 'red',label = 'MLKLp_cal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("MLKLp [Molecules/Cell]", fontsize=15)
# plt.title('MLKLp Trajectories')
# plt.legend(loc ='best')x
plt.show()

#
# print(L3.observables['C8a_obs'][:])
# print(L4.observables['RIP13_obs'][:])
#
# quit()
#
# plt.figure(figsize = (18,7))
# # # plt.figure()
# plt.subplot(231)
# # # plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
# # # plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
# plt.plot(tspan/60, L3.observables['TNF_obs'],label = 'TNF')
# plt.plot(tspan/60, L4.observables['TNF_obs'],label = 'TNFcal')
# #
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# #
# # # plt.figure()
# plt.subplot(232)
# # # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # # plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
# # # plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
# plt.plot(tspan/60, L3.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, L4.observables['CI_k63_obs'],label = 'CI_k63cal')
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # #
# # # plt.figure()
# plt.subplot(233)
# # plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# # plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po4_obs')
# plt.plot(tspan/60, L3.observables['RIP13_obs'],label = 'RIP13')
# plt.plot(tspan/60, L4.observables['RIP13_obs'],label = 'RIP13cal')
#
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(234)
# # plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub_obs')
# # plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k63_obs')
# plt.plot(tspan/60, L3.observables['CI'],label = 'CI')
# plt.plot(tspan/60, L4.observables['CI'],label = 'CIcal')
# # plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# plt.legend(loc = 0)
#
# plt.subplot(235)
# # plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# # plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# plt.plot(tspan/60, L3.observables['C8a_obs'],label = 'C8a')
# plt.plot(tspan/60, L4.observables['C8a_obs'],label = 'C8acal')
# # plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# plt.legend(loc = 0)
#
#
# plt.subplot(236)
# # plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# # plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLp')
# # plt.plot(tspan/60, L4.observables['CIIa_obs'],label = 'CIIa')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# plt.legend(loc = 0)
#
#
#
# plt.figure(figsize = (18,7))
# # # plt.figure()
# plt.subplot(221)
# # # plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
# # # plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
# plt.plot(tspan/60, L3.observables['A20_obs'],label = 'A20')
# plt.plot(tspan/60, L4.observables['A20_obs'],label = 'A20cal')
# #
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# #
# # # plt.figure()
# plt.subplot(222)
# # # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # # plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
# # # plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
# plt.plot(tspan/60, L3.observables['Fadd_obs'],label = 'Fadd')
# plt.plot(tspan/60, L4.observables['Fadd_obs'],label = 'Faddcal')
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # #
# # # plt.figure()
# plt.subplot(223)
# # plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# # plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po4_obs')
# plt.plot(tspan/60, L3.observables['Tradd_obs'],label = 'Tradd')
# plt.plot(tspan/60, L4.observables['Tradd_obs'],label = 'Traddcal')
#
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(224)
# # plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub_obs')
# # plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k63_obs')
# plt.plot(tspan/60, L3.observables['C8i_obs'],label = 'C8a')
# plt.plot(tspan/60, L4.observables['C8i_obs'],label = 'C8acal')
# # plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# plt.legend(loc = 0)
#
#
# plt.tight_layout()
# plt.show()
# #

