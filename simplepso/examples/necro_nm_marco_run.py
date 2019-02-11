from pylab import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
from pysb.simulator import CupSodaSimulator
from pysb.simulator.scipyode_test import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator
import necro_mn_marco as model

tspan = np.linspace(0, 1440, 1441)
sim = ScipyOdeSimulator(model, tspan=tspan, log_level=False)
result = sim.run()

# # for p in model.parameters:
# #     print('{},{:e}'.format(p.name,p.value))

# #
# # with open('params_necro.txt', 'w') as f:
# #     for p, v in zip(model.parameters, result.param_values[0]):
# #         f.write('{},{:e}\n'.format(p.name, v))
# #
# # quit()
# #
# #
# # print(result.param_values)
# # quit()
# # for p in result.param_values:
# #     print('{},{:e}'.format(p.name, p.value))
# pars =  np.load('optimizer_best_50_50.npy')
# rate_params = model.parameters_rules() # these are only the parameters involved in the rules
# param_values = np.array([p.value for p in model.parameters]) # these are all the parameters
# rate_mask = np.array([p in rate_params for p in model.parameters])
# param_values[rate_mask] = 10 ** pars


# tspan = np.linspace(0, 1400, 1441)
# sim = BngSimulator(model, tspan=tspan)
# result = sim.run(method='ode')
#
#
# plt.figure()
# # plt.plot(tspan/60, result.observables['MLKLu'][:])
# plt.plot(tspan/60, result.observables['MLKLa_obs'][:])
# # plt.scatter(x, mlkl)
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.ylim(ymax = )
# plt.show()
#


# tnf = [2326, 232, 23,2]
# color = ['r', 'm', 'g', 'b']
# lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
#
# tspan = np.linspace(0, 1440, 1441)
# sim = BngSimulator(model, tspan=tspan)
# result = sim.run(method = 'ode',param_values=param_values,initials= {TNF(brec=None):tnf})
# df = result.dataframe
#
# # plt.figure()
# # plt.plot(tspan/60, result.observables['MLKLa_obs'], color = 'tab:gray', marker = 's')
# # plt.xlabel('Time [minutes]', fontsize=16)
# # plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.legend('100 ng/ml ' , title = 'TNF FP', loc=0, fontsize = 5)
# # plt.show()
# #
# #
# # quit()
# # mlkl = [0, 170, 900, 4880, 9940, 10000]
# # x = [0, 60, 120, 240, 360, 480]
#
# plt.figure()
# for n in range(0,4):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan/60, df.loc[n]['MLKLa_obs'].iloc[:], '--', c = color[n],lw = 1.5) #fppf
# # plt.scatter(x, mlkl, color = 'tab:gray', marker = 's')
# plt.xlabel('Time [hours]', fontsize=14)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=14)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF Doses', loc=0, fontsize = 10)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# plt.show()

# plt.figure()
# plt.plot(tspan/60, result.observables['MLKLu'][:])
# plt.plot(tspan/60, result.observables['MLKLa_obs'][:])
# # plt.scatter(x, mlkl)
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.ylim(ymax = )
# plt.show()
