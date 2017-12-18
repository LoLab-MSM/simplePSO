from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.backends.backend_pdf import PdfPages
# from QQSB.numtools import simulatortools as st
from pysb.util import alias_model_components
from new_necro_pso import model


# pars = np.array([2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 4800,
#           2.26216506e-06,  5.03882772e-04,   1.71088800e-03,   9.25954827e-07,
#            4.10179461e-04,  1.21120342e-06,   1.04431157e-04,   2.36918819e-06,
#            5.99145732e-03,  1.18073726e-06,   5.86948599e-04,   9.70876569e-03,
#            2.57629434e-06,  4.27412324e-03,   2.68048512e-06,   3.12504841e-03,
#            9.15814868e-03,  5.46825036e-07,   8.04641467e-05,   2.12995771e-02,
#            1.04942303e+00,  7.47419785e-07,   1.08122555e-04,   6.66971920e-03,
#            1.82947728e-05,  2.26036706e-02,   4.90003351e-05,   1.29291276e-02,
#            1.61970726e+00,  1.35343038e-06,   7.77429364e-03,   1.80132522e-02,
#            1.44832257e-02,  1.05449953e-02,   3.60566804e-02,   7.40500933e-06,
#            4.96307811e+01, 1.0])

model.enable_synth_deg()
# print(len((model.parameters)))
# print(len(pars))
# print(pars.shape)
# quit()

tspan = np.linspace(0, 720, 721)
sim1 = ScipyOdeSimulator(model, tspan=tspan)
# sim2 = ScipyOdeSimulator(model, tspan=tspan)
L4 = sim1.run()
# L5 = sim2.run()

plt.figure(figsize = (18,7))
# # plt.figure()
plt.subplot(231)
# # plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
# # plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
# plt.plot(tspan/60, L3.observables['TNF_obs'],label = 'TNF10d')
plt.plot(tspan/60, L4.observables['TNF_obs'],label = 'TNF100s')
#
# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
#
# # plt.figure()
plt.subplot(232)
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
# # plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
# plt.plot(tspan/60, L3.observables['CI_k63_obs'],label = 'CI_k6310d')
plt.plot(tspan/60, L4.observables['CI_k63_obs'],label = 'CI_k63_obs')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
# #
# # plt.figure()
plt.subplot(233)
# plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# plt.plot(tspan/60, L2.observables['RIP13po4_obs'],label = 'RIP13po41')
# plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po410d')
plt.plot(tspan/60, L4.observables['RIP13_obs'],label = 'RIP13_obs')

# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
# plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa10d')
# plt.plot(tspan/60, L4.observables['CIIa_obs'],label = 'CIIa')
plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa100s')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
plt.legend(loc = 0)
#
plt.tight_layout()
plt.show()
#
# plt.figure()
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],color = 'black',label = 'MLKLp')
# plt.plot(tspan/60, L5.observables['MLKLa_obs'],color = 'black',label = 'MLKLp_cal')
# plt.plot(tspan/60, L4.observables['MLKL_obs'],color = 'red',label = 'MLKL')
# plt.plot(tspan/60, L5.observables['MLKL_obs'],color = 'red',label = 'MLKL_cal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# plt.legend(loc = 0)
# plt.show()
