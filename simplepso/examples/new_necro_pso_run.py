from pysb.core import *
# from pysb.integrate import Solver
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.backends.backend_pdf import PdfPages
# from QQSB.numtools import simulatortools as st
from pysb.util import alias_model_components
from correct_necro_molecules import model

from necro_pso_kd import display


# pars =  np.load('optimizer_best_10000.npy')
# print(10 ** pars)
# quit()

# pars = np.array([
# -3.000000000000000000e+00,
# -3.090842869377300683e+00,
# -4.535733029499747992e+00,
# -3.001779214881634683e+00,
# -4.571363867887905386e+00,
# -3.000000000000000000e+00,
# -6.000000000000000000e+00,
# -3.000000000000000000e+00,
# -5.982958244298056449e+00,
# -3.000063563793294641e+00,
# -3.595271259531686026e-01,
# 1.785839987372437765e+00 ,
# -3.000007556571619105e+00,
# -3.333304995414421157e+00,
# -3.000000000000000000e+00,
# -5.993713416826150109e+00,
# -4.000000000000000000e+00,
# -5.379370708464467299e+00,
# -4.258262351082532149e+00,
# 2.000000000000000000e+00 ,
# -7.747894247036030135e-02,
# -6.408648675403814110e+00,
# -4.237107406033835844e+00,
# -2.090564604636510726e+00,
# -6.780371356918183601e+00,
# 1.255272505103306013e+00 ,
# -1.901892602467097459e+00,
# -2.185575801983683064e+00,
# 1.243192385359765284e+00 ,
# -5.288698622615229716e+00,
# -5.249132838452612582e+00,
# -9.740084403570974381e-01,
# -2.238894808492042365e+00,
# 0.000000000000000000e+00 ,
# -2.464375752635486805e+00,
# -6.033114126716795589e+00,
# 3.717180930857924004e-01])
# pars = np.array([
# -3.000000000000000000e+00 ,
# -1.252317339283543607e+00 ,
# -5.082789313249298502e+00 ,
# -3.003189290936787437e+00 ,
# -3.780149772707062983e+00 ,
# -3.309059289555452033e+00 ,
# -6.000000000000000000e+00 ,
# -3.000072740009938066e+00 ,
# -4.136077970611998111e+00 ,
# -3.034915454986220507e+00 ,
# -3.038551754556173101e+00 ,
# 1.563556082282690163e+00  ,
# -7.266248364566838092e+00 ,
# -3.997605834854700113e+00 ,
# -3.004439232863266973e+00 ,
# -3.974786922462485439e+00 ,
# -3.999999999718301780e+00 ,
# -4.817018895571956172e+00 ,
# -2.109962334946007267e+00 ,
# 1.418008368291111543e+00  ,
# -4.311788145896281454e-01 ,
# -5.877428636327199918e+00 ,
# -3.311136310426499474e+00 ,
# -1.796026494183759370e+00 ,
# -8.485134358575241009e+00 ,
# -2.422510443920192991e+00 ,
# -2.062986040532480914e+00 ,
# -3.019018272317760765e+00 ,
# -1.819827882296911836e+00 ,
# -5.151351274439876171e+00 ,
# -3.351476569402714212e+00 ,
# 1.768521394307117189e+00  ,
# 9.990576123462765468e-01  ,
# 0.000000000000000000e+00  ,
# -1.717319136311646899e+00 ,
# -5.396518489284537701e+00 ,
# 1.889321233519191789e+00
# ])

# pars = np.array([
# -3.000000000000000000e+00,
# -5.190463909336157755e+00,
# -6.000000000000000000e+00,
# -3.033085665199291014e+00,
# -4.616053825158105361e+00,
# -6.356697247845557897e+00,
# -2.168782791959026124e+00,
# -3.058696627754011299e+00,
# -3.205201740640312824e+00,
# -3.000000000000000000e+00,
# -2.098242182818893742e+00,
# 1.251276896008965078e+00 ,
# -3.360171520761877861e+00,
# -4.443801585952002142e+00,
# -3.002469573770719879e+00,
# -1.516034355984396509e+00,
# -3.999999762801250913e+00,
# -4.253324882630691306e+00,
# -3.307450965106528251e+00,
# 1.542167033113567953e+00 ,
# 1.958004696407286138e+00 ,
# -8.009639890723171618e+00,
# -2.892524316001871654e+00,
# -1.688715707452753456e+00,
# -8.332899618496059091e+00,
# 5.613842426630625271e-01 ,
# -3.393485506662298756e-01,
# -2.384168876963753902e+00,
# -1.841615835076017449e+00,
# -5.609768825686679783e+00,
# -1.580335592915386167e+00,
# 8.935245594185757811e-01 ,
# 1.000000000000000000e+00 ,
# 0.000000000000000000e+00 ,
# -3.454070976893416578e+00,
# -5.942968322097673450e+00,
# 3.362295758042641580e-01])

# pars = np.array([-3.90543513, -3.84631005, -2.93469037, -4.97554492, -5.5154541 , -6.08242418,
#                  -2.86229483, -6.80280754, -3.89629432, -4.91123185, -2.8854796 , -0.33340022,
#                  -4.45343191, -2.17690192, -4.50299611, -4.03538748, -4.        , -5.58059923,
#                  -2.81792633,  0.86244067, -2.06393558, -8.55032924, -5.30963198, -0.61639341,
#                  -6.98793546, -1.10961786, -1.75126626, -2.69282528, -2.24177105, -5.90212771,
#                  -4.58622347, -0.40257703, -1.09440659, -2.53660987, -2.62267299, -5.83944408,
#                   0.49786203])

# params = np.array([-3.00000000e+00,  -2.78247317e+00,  -3.43246928e+00,  -5.05877799e+00,
#                   -3.79144230e+00 , -6.24937055e+00 , -2.48675468e-07 , -5.92891152e+00 ,
#                   -2.68086333e+00 , -3.00000000e+00 , -4.25001012e+00 ,  1.48470006e+00 ,
#                   -3.00301837e+00 , -1.66052896e-03 , -5.10845808e+00 , -1.93227937e+00 ,
#                    1.22585685e+00 , -3.00036514e+00 , -4.86410582e+00 ,  1.54364832e+00 ,
#                    1.73184292e+00 , -7.83377682e+00 , -3.36377726e+00 , -2.52208917e+00 ,
#                   -7.55153153e+00 , -2.21524249e+00 , -1.43141578e+00 , -1.85217468e+00 ,
#                   -5.72547202e-01 , -5.26414183e+00 , -2.51130609e+00 ,  2.00000000e+00 ,
#                   -2.96841789e+00 , -1.57596579e-05 , -2.00187299e+00 , -6.05891821e+00 ,
#                    2.76071207e+00])

# params = np.array([2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 10000,
#                    1.00000000e-03,   1.65016294e-03,  3.69428775e-04,   8.73417743e-06,
#                    1.61643297e-04,   5.63156952e-07,  9.99999427e-01,   1.17784591e-06,
#                    2.08514696e-03,   1.00000000e-03,  5.62328222e-05,   3.05281200e+01,
#                    9.93074042e-04,   9.96183791e-01,  7.79008004e-06,   1.16874733e-02,
#                    1.68211952e+01,   9.99159587e-04,  1.36739561e-05,   3.49661906e+01,
#                    5.39315522e+01,   1.46630117e-08,  4.32735715e-04,   3.00545915e-03,
#                    2.80846147e-08,   6.09196655e-03,  3.70326013e-02,   1.40548210e-02,
#                    2.67579475e-01,   5.44324861e-06,  3.08101569e-03,   1.00000000e+02,
#                    1.07542991e-03,   9.99963713e-01,  9.95696568e-03,   8.73135789e-07,
#                    5.76384204e+02,   1.0])

# pars = 10**params
# print(pars)
# quit()
display(pars)
quit()

# params = np.array([  2326, 4800, 9000, 40000, 9000, 9000, 2070, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 10000,
#                      1.39827142e-05,   7.84751944e-03,   9.51836183e-05,   5.88487293e-05,
#                        1.60614371e-02,   9.40268203e-06,   8.08999401e-02,   8.44124820e-06,
#                        2.33046519e-02,   9.93187511e-06,   2.10323727e-02,   3.36544037e-03,
#                        4.92553920e-07,   4.83968928e-06,   2.85503982e-04,   2.09928228e-02,
#                        4.91236966e-03,   4.96880920e-05,   3.65153006e-04,   2.59328125e-01,
#                        1.04958786e-02,   2.04207190e-07,   4.20427572e-06,   6.13573748e-01,
#                        9.58225778e-06,   1.24798053e-03,   1.80946068e-02,   4.63841711e+00,
#                        5.83766488e-01,   1.92165265e-06,   6.91656591e-04,   7.05279272e-01,
#                        1.18786671e-03,   6.57210828e-02,   4.39462791e-04,   4.20074810e-08,
#                        7.12748086e+01, 1.0])

# pars = np.array([-5.46774398, -2.60184267, -3.12073463, -4.31653794, -2.23669114, -5.28107178,
#                  -3.79688244, -4.27828926, -1.95684989, -6.52415156, -4.03994837, -0.7192964 ,
#                  -5.96853617, -3.77692131, -6.34104252, -3.69201844, -1.49511334, -4.20121292,
#                  -2.63041475, -1.24124118, -0.58964676, -7.3424675 , -6.16308612, -2.64858185,
#                  -5.35010387, -1.36125674, -2.21808402, -1.99465669, -0.79578398, -4.69971761,
#                  -2.99007416, -1.75841718, -1.90863638, -2.03532575, -3.99256982, -6.20617955,
#                   0.24089193])

#NEWEST
# pars = np.array([  2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 10000,
#                     3.40608922e-06,   2.50125132e-03,   7.57295489e-04,   4.82460830e-05,
#                    5.79840919e-03,   5.23513903e-06,   1.59631120e-04,   5.26878819e-05,
#                    1.10446030e-02,   2.99122058e-07,   9.12119268e-05,   1.90855025e-01,
#                    1.07513705e-06,   1.67139343e-04,   4.55992269e-07,   2.03227072e-04,
#                    3.19806039e-02,   6.29197633e-05,   2.34199115e-03,   5.73797722e-02,
#                    2.57248730e-01,   4.54498548e-08,   6.86932209e-07,   2.24604343e-05,
#                    4.46576772e-06,   4.35254490e-03,   6.05223775e-03,   3.01237942e-02,
#                    1.60035385e-01,   1.99656011e-02,   1.02311827e-03,   1.74414593e-02,
#                    1.23413770e-02,   9.21879696e-03,   1.01725581e-04,   6.22043061e-07,
#                    1.74137350e+00, 1.0])

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
L4 = sim1.run(param_values=params)
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
