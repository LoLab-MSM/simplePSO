from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.backends.backend_pdf import PdfPages
# from QQSB.numtools import simulatortools as st
from pysb.util import alias_model_components
from correct_necro_molecules import model


# pars = np.array([-5.46774398, -2.60184267, -3.12073463, -4.31653794, -2.23669114, -5.28107178,
#                  -3.79688244, -4.27828926, -1.95684989, -6.52415156, -4.03994837, -0.7192964 ,
#                  -5.96853617, -3.77692131, -6.34104252, -3.69201844, -1.49511334, -4.20121292,
#                  -2.63041475, -1.24124118, -0.58964676, -7.3424675 , -6.16308612, -2.64858185,
#                  -5.35010387, -1.36125674, -2.21808402, -1.99465669, -0.79578398, -4.69971761,
#                  -2.99007416, -1.75841718, -1.90863638, -2.03532575, -3.99256982, -6.20617955,
#                   0.24089193])

pars = np.array([  2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 10000,
                    3.40608922e-06,   2.50125132e-03,   7.57295489e-04,   4.82460830e-05,
                   5.79840919e-03,   5.23513903e-06,   1.59631120e-04,   5.26878819e-05,
                   1.10446030e-02,   2.99122058e-07,   9.12119268e-05,   1.90855025e-01,
                   1.07513705e-06,   1.67139343e-04,   4.55992269e-07,   2.03227072e-04,
                   3.19806039e-02,   6.29197633e-05,   2.34199115e-03,   5.73797722e-02,
                   2.57248730e-01,   4.54498548e-08,   6.86932209e-07,   2.24604343e-05,
                   4.46576772e-06,   4.35254490e-03,   6.05223775e-03,   3.01237942e-02,
                   1.60035385e-01,   1.99656011e-02,   1.02311827e-03,   1.74414593e-02,
                   1.23413770e-02,   9.21879696e-03,   1.01725581e-04,   6.22043061e-07,
                   1.74137350e+00, 1.0])

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

# pars = np.array([-5.05470813, -3.1580303 , -3.11056345, -5.03792988, -4.07031751, -5.86592676,
#                  -2.41246963, -6.47791499, -2.01955062, -5.76756868, -3.62378031, -1.42306076,
#                  -5.68659694, -2.70665402, -5.44713892, -2.33097527, -1.83837442, -5.07432414,
#                  -2.88651406, -0.28211213, -0.63682123, -6.96521154, -5.46102854, -3.06642727,
#                  -6.21093252, -2.17821161, -2.21935119, -1.58879521,  0.17259997, -5.6894568 ,
#                  -3.52767876, -1.39768111, -1.4625396 , -2.56810944, -3.56284713, -6.95863975,
#                   1.74242706])

# pars = np.array([2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 4800,
#                    8.81641186e-06,   6.94975829e-04,   7.75240673e-04,   9.16368433e-0,
#                    8.50516004e-05,   1.36167430e-06,   3.86839105e-03,   3.32724675e-0,
#                    9.55981262e-03,   1.70777763e-06,   2.37804293e-04,   3.77519370e-0,
#                    2.05779951e-06,   1.96492501e-03,   3.57158574e-06,   4.66685954e-0,
#                    1.45086024e-02,   8.427055AZ61e-06,   1.29863152e-03,   5.22261330e-0,
#                    2.30769692e-01,   1.08339907e-07,   3.45916645e-06,   8.58168818e-0,
#                    6.15272465e-07,   6.63419740e-03,   6.03460447e-03,   2.57753630e-0,
#                    1.48798985e+00,   2.04429328e-06,   2.96702523e-04,   4.00238526e-0,
#                    3.44715173e-02,   2.70327707e-03,   2.73623170e-04,   1.09991785e-0,
#                    5.52620587e+01, 1.0])
# print(list(pars).T)
# quit()


# pars = np.array([ 2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 20000,
#                    6.12461255e-04,   2.56036122e-05,   2.87222481e-03,   9.99480050e-04,
#                    2.45549486e-02,   3.44124632e-04,   1.63991153e-04,   2.26485264e-04,
#                    1.77349133e-03,   1.20006406e-05,   9.03749209e-05,   6.18265524e-01,
#                    3.05752763e-05,   2.70116659e-02,   1.12851945e-05,   5.96351986e-05,
#                    1.88904037e+00,   4.61354680e-06,   3.97950001e-05,   3.46268679e-02,
#                    7.04525359e-06,   2.66858598e-04,   2.30427352e-01,   1.11064373e-03,
#                    2.20415931e-08,   1.27541711e-04,   2.22669026e-03,   6.95847514e-03,
#                    7.55415219e-03,   2.10637834e-07,   5.39741523e-03,   2.37668842e-04,
#                    6.71727069e-02,   7.51875355e-04,   1.14508504e+00,   1.00000000e+01])
# #
# pars = np.array([-3.21292138, -4.59169876, -2.54178157, -3.00022587, -1.60986097, -3.46328424,
#                  -3.78517958, -3.64496005, -2.75117093, -4.92079557, -4.04395207, -0.20882497,
#                  -4.51462961, -1.56844863, -4.94749095, -4.22449733,  0.27624124, -5.33596507,
#                  -4.40017149, -1.46058679, -5.15210337, -3.5737188 , -0.63746597, -2.95442523,
#                  -7.65675702, -3.89434776, -2.65234019, -2.15748592, -2.12181427, -6.67646362,
#                  -2.26781417, -3.62402775, -1.17280715, -3.12385415,  0.05883774, 1.0])

# new_params = 10**pars
# print(new_params)
# quit()


model.enable_synth_deg()
# print(len((model.parameters)))
# print(len(pars))
# # print(pars.shape)
# quit()

tspan = np.linspace(0, 1200, 1201)
sim1 = ScipyOdeSimulator(model, tspan=tspan)
sim2 = ScipyOdeSimulator(model, tspan=tspan)
L4 = sim1.run(param_values=pars)
L3 = sim2.run()
# sim_df = L3.dataframe

# print(L3.observables['C8a_obs'][:])
# print(sim_df.observables['C8i_obs'].iloc[:])

# quit()

plt.figure(figsize = (18,7))
# # plt.figure()
plt.subplot(231)
# # plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
# # plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
plt.plot(tspan/60, L3.observables['TNF_obs'],label = 'TNF')
plt.plot(tspan/60, L4.observables['TNF_obs'],label = 'TNFcal')
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
plt.plot(tspan/60, L3.observables['CI_k63_obs'],label = 'CI_k63')
plt.plot(tspan/60, L4.observables['CI_k63_obs'],label = 'CI_k63cal')
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
# plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po4_obs')
plt.plot(tspan/60, L3.observables['RIP13_obs'],label = 'RIP13')
plt.plot(tspan/60, L4.observables['RIP13_obs'],label = 'RIP13cal')

# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
# plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub_obs')
# plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k63_obs')
plt.plot(tspan/60, L3.observables['CI'],label = 'CI')
plt.plot(tspan/60, L4.observables['CI'],label = 'CIcal')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
plt.legend(loc = 0)

plt.subplot(235)
# plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# plt.plot(tspan/60, L3.observables['C8a_obs'],label = 'C8a')
plt.plot(tspan/60, L4.observables['C8a_obs'],label = 'C8acal')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
plt.legend(loc = 0)


plt.subplot(236)
# plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLp')
# plt.plot(tspan/60, L4.observables['CIIa_obs'],label = 'CIIa')
plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
plt.legend(loc = 0)



plt.figure(figsize = (18,7))
# # plt.figure()
plt.subplot(221)
# # plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
# # plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
plt.plot(tspan/60, L3.observables['A20_obs'],label = 'A20')
plt.plot(tspan/60, L4.observables['A20_obs'],label = 'A20cal')
#
# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
#
# # plt.figure()
plt.subplot(222)
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
# # plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
plt.plot(tspan/60, L3.observables['Fadd_obs'],label = 'Fadd')
plt.plot(tspan/60, L4.observables['Fadd_obs'],label = 'Faddcal')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
# #
# # plt.figure()
plt.subplot(223)
# plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po4_obs')
plt.plot(tspan/60, L3.observables['Tradd_obs'],label = 'Tradd')
plt.plot(tspan/60, L4.observables['Tradd_obs'],label = 'Traddcal')

# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(224)
# plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub_obs')
# plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k63_obs')
plt.plot(tspan/60, L3.observables['C8i_obs'],label = 'C8a')
plt.plot(tspan/60, L4.observables['C8i_obs'],label = 'C8acal')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLpcal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("Molecules/Cell", fontsize=15)
plt.legend(loc = 0)


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
