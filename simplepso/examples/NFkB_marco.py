__author__ = 'geena'

from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from sympy import sympify
from pysb.bng import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator


Model ()

#added monomers
Monomer('TNF_ext')
Monomer('TNFR1', ['state'], {'state': ['i', 'a']})
Monomer('IKKK', ['state'], {'state': ['a', 'n']})
Monomer('IKK', ['state'], {'state': ['n', 'a', 'i', 'ii']})
Monomer('NFkB', ['b','loc'], {'loc': ['n', 'c']})
Monomer('IkBa_gene', ['state'], {'state': ['on', 'off']})
Monomer('IkBa_mRNA')
Monomer('IkBa', ['b', 'phos', 'loc'], {'phos': ['u', 'p'], 'loc': ['n', 'c']})
Monomer('A20_gene', ['state'], {'state': ['on', 'off']})
Monomer('A20_mRNA')
Monomer('A20')

Parameter('KNN', 2.0e5)
Parameter('KN', 1.0e5)
Parameter('M', 2000)
Parameter('AN', 2)
Parameter('ANa', 2)


#Declaring initial conditions
Initial(IKKK(state = 'a'), Parameter('IKKKa_0'))
Initial(IKK(state = 'n'), Parameter('IKKn_0', 2e5))
Initial(IKK(state = 'a'), Parameter('IKKa_0'))
Initial(IKK(state = 'i'), Parameter('IKKi_0'))
Initial(IkBa(b = None, phos = 'p', loc = 'c'), Parameter('IkBap_0'))
Initial(NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c'), Parameter('IkBap_NFkB'))
Initial(NFkB(b = None, loc = 'c'), Parameter('NFkBc'))
Initial(NFkB(b = None, loc = 'n'), Parameter('NFkBn_0', 1))
Initial(A20(), Parameter('A20_0', 10000))
Initial(A20_mRNA(), Parameter('A20t_0', 10))
Initial(IkBa(b = None, phos = 'u', loc = 'c'), Parameter('IkBa_0', 0.14*100000))
Initial(IkBa(b = None, phos = 'u', loc = 'n'), Parameter('IkBan_0', 0.06*100000))
Initial(IkBa_mRNA(), Parameter('IkBat_0', 10))
Initial(NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c'), Parameter('IkBa_NFkB_0', 100000))
Initial(NFkB(b = 1, loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n'), Parameter('IkBan_NFkBn'))
Initial(TNFR1(state = 'a'), Parameter('TNFR1a'))
Initial(A20_gene(state = 'on'), Parameter('A20_on_0'))
Initial(IkBa_gene(state = 'on'), Parameter('IkBa_on_0'))
Initial(TNF_ext(), Parameter('TNF_ext_0', 0.1))
Initial(IKKK(state = 'n'), Parameter('IKKKn_0', KN.value))
Initial(IKK(state = 'ii'), Parameter('IKKii'))
Initial(TNFR1(state = 'i'), Parameter('TNFR1i_0', M.value))
Initial(A20_gene(state = 'off'), Parameter('A20_off_0', AN.value))
Initial(IkBa_gene(state = 'off'), Parameter('ANa_0', ANa.value))

#Declaring Parameters
Parameter('kv', 5.0) #Nucealr to cytoplasm volume
Parameter('kb', 1.2e-5) #receptor activation rate
Parameter('kf', 1.2e-3) #receptor inactivation rate
Parameter('Tdeg', 7.7e-4) #TNF loss same as cdeg
Parameter('ka',1e-5) #IKKK kinase activation rate
Parameter('ka20', 1e5) #A20 TNFR1 block
Parameter('ki', 0.01) #IKKK kinase inactivation rate
Parameter('k1', 2.0*6e-10) #IKKn activation by IKKK
Parameter('k3', 0.002) #IKKa inactivation by A20
Parameter('k2', 10000) #IKKa inactivation
Parameter('k4', 0.001) #IKKii transfer rate
Parameter('c1a', 0.1) #inducible IkBa mRNA synthesis
Parameter('a2', 1e-7) #IkBa phosphorylation b/c IKKa
Parameter('tp', 0.01) #degradation of phosph-IkBa complex with NFkB
Parameter('a3', 5e-7) #IkBa_NFkB phosphorylation b/c IKKa
Parameter('a1', 5e-7) #IkBa*NFkB association
Parameter('c6a', 0.00002) #spontaneous IkBa_NFkB defg of IkBa complex to NFkB
Parameter('c5a', 0.0001) #IkBa deg rate
Parameter('i1a', 0.002) #IkBa nuclear import
Parameter('e1a', 0.005) #IkBa nuclear export
Parameter('e2a', 0.05) #IkBa_NFkB nuclear export
Parameter('i1', 0.01) #NFkB nuclear import
Parameter('c4', 0.5) #A20, IkBa transformation rate
Parameter('c5', 0.0005) #A20 degredation rate
Parameter('c1', 0.1) #inducible A20 mRNA synthesis
Parameter('c3', 0.00075) #A20 and IkBa mRNA deg rate
Parameter('q1', 4e-7) #NFkB attaching @ A20 and IkBa site
Parameter('q2', 1e-6) #IkBa inducible detaching from A20, IkBa site
Parameter('k3_div_k2', k3.value/k2.value) #for IKKa and A20
Parameter('a1_mult_kv', a1.value*kv.value) #for volume IkBa association NFkB

#Declaring observables
Observable('IKKKa_obs', IKKK(state = 'a'))
Observable('IKKKn_obs', IKKK(state = 'n'))
Observable('IKKa_obs', IKK(state = 'a'))
Observable('IKKn_obs', IKK(state = 'n'))
Observable('IKKi_obs', IKK(state = 'i'))
Observable('IKKii_obs', IKK(state = 'ii'))
Observable('TNF_ext_obs', TNF_ext())
Observable('IkBap_obs', IkBa(b = None, phos = 'p', loc = 'c'))
Observable('IkBan_NFkBn_obs', NFkB(b = 1,loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n'))
Observable('IkBap_NFkBc_obs', NFkB(b = 1,loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c'))
Observable('NFkBc_obs', NFkB(b = None, loc = 'c'))
Observable('NFkBn_obs', NFkB(b = None, loc = 'n'))
Observable('Nuclear_NFkBn', NFkB(b = 1,loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n') + NFkB(b = None, loc = 'n'))
Observable('A20_obs', A20())
Observable('TNFR1a_obs', TNFR1(state = 'a'))
Observable('TNFR1i_obs', TNFR1(state = 'i'))
Observable('A20_off_obs', A20_gene(state = 'off'))
Observable('IkBa_off_obs', IkBa_gene(state = 'off'))
Observable('A20_on_obs', A20_gene(state = 'on'))
Observable('IkBa_on_obs', IkBa_gene(state = 'on'))
Observable('A20t_obs', A20_mRNA())
Observable('IkBac_obs', IkBa(b = None, phos = 'u', loc = 'c'))
Observable('IkBan_obs', IkBa(b = None, phos = 'u', loc = 'n'))
Observable('IkBat_obs', IkBa_mRNA())
Observable('IkBa_NFkB_obs', NFkB(b = 1,loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c'))

#Declaring expression
Expression('keff', ka*ka20/(ka20+A20_obs)) #10000 #michaelis menten

#Declaring rules
Rule('TNF_activate_TNFR1', TNF_ext() + TNFR1(state = 'i') >> TNF_ext() + TNFR1(state = 'a'), kb) #TNFa binding to TNFR1
Rule('TNFR1a_creates_IKKKa', TNFR1(state = 'a') + IKKK(state = 'n') >> TNFR1(state = 'a') + IKKK(state = 'a'), keff) #TNFR1 activating IKKK
Rule('inactivate_TNFR1', TNFR1(state = 'a') >> TNFR1(state = 'i'), kf) #Inactivation of TNFR1
Rule('IKKKa_to_IKKKn', IKKK(state = 'a') >> IKKK(state = 'n'), ki) #IKKKa creates KN (total amount of IKKK kinase molecules)
Rule('IKKKa_activates_IKK', IKKK(state = 'a') + IKKK(state = 'a') + IKK(state = 'n') >> IKKK(state = 'a') + IKKK(state = 'a') + IKK(state = 'a'), k1) #IKKK activates IKK
Rule('IKKa_to_IKKi', IKK(state = 'a') >> IKK(state = 'i'), k3) #IKK inactivation
Rule('IKKa_and_A20', IKK(state = 'a') + A20() >> IKK(state = 'i') + A20(), k3_div_k2) #Exp1) #A20 mediated IKKa to IKKi
Rule('IKKi_to_IKKii', IKK(state = 'i') >> IKK(state = 'ii'), k4) #inactive IKK inactivates immediate IKK
Rule('TNF_ext_deg', TNF_ext() >> None, Tdeg) #degredation of TNF
Rule('IkBapc_NFkBpc_to_NFkB_c', IkBa(b = 1, phos = 'p', loc = 'c') % NFkB(b = 1,  loc = 'c') >> NFkB(b = None, loc = 'c'), tp) #IkBa|NFkB phos complex cyto to free NFkB in cyto
Rule('IkBa_p_deg', IkBa(b = None, phos = 'p', loc = 'c') >> None, tp) #degredation of phos IkBa
Rule('A20_deg', A20() >> None, c5) #degrade A20
Rule('A20_mRNA_deg', A20_mRNA() >> None, c3) #degrade A20t
Rule('IkBa_cyt_deg', IkBa(b = None, phos = 'u', loc = 'c') >> None, c5a) #IkBa cyto degredation
Rule('IkBa_NFkB_to_NFkB_c', IkBa(b = 1, phos = 'u', loc = 'c') % NFkB(b = 1, loc = 'c') >> NFkB(b = None, loc = 'c'), c6a) #release of NFkB in cyto
Rule('IkBa_mRNA_deg', IkBa_mRNA() >> None, c3) #degredation of IkBa_mRNA
Rule('NFkB_c_to_NFkBn', NFkB(b = None, loc = 'c') >> NFkB(b = None, loc = 'n'), i1) #Transport of NFkB from cyto to nuc
Rule('NFkB_c_and_IkBa', NFkB(b = None, loc = 'c') + IkBa(b = None, phos = 'u', loc = 'c') >> NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c') , a1) #NFkB in cyto and IkBa in cyto create bind in cyto
Rule('NFkBn_and_IkBan', NFkB(b = None, loc = 'n') + IkBa(b = None, phos = 'u', loc = 'n') >> NFkB(b = 1, loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n'), a1_mult_kv) #NFkBn and IkBan create IkBan|NFkBn complex in nucleus
Rule('IKKa_and_IkBa_NFkB', IKK(state = 'a') + NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c') >> IKK(state = 'a') + NFkB(b = 1,  loc = 'c') % IkBa(b = 1, phos = 'p', loc = 'c'), a3) #IKKa creates IkBa|NFkB phos complex in cyto
Rule('IKKa_and_IkBap', IKK(state = 'a') + IkBa(b = None, phos = 'u', loc = 'c') >> IKK(state  = 'a') + IkBa(b = None, phos = 'p', loc = 'c'), a2) #IKKa phos IkBa
Rule('A20_mRNA_to_A20', A20_mRNA() >> A20_mRNA() + A20(), c4) #A20 translation
Rule('A20_gene_to_A20_mRNA', A20_gene(state = 'on') >> A20_gene(state ='on') + A20_mRNA(), c1) #A20 transcription
Rule('IkBa_mRNA_to_IkBa', IkBa_mRNA() >> IkBa_mRNA() + IkBa(b = None, phos = 'u', loc = 'c'), c4) #IkBa translation
Rule('IkBac_to_IkBan', IkBa(b = None, phos = 'u', loc = 'c') >> IkBa(b = None, phos = 'u', loc = 'n'), i1a) #IkBa cyto to nuc
Rule('IkBan_to_IkBac', IkBa(b = None, phos = 'u', loc = 'n') >> IkBa(b = None, phos = 'u', loc = 'c'), e1a) #IkBa nuc to cyto
Rule('IkBa_on_to_IkBat', IkBa_gene(state = 'on') >> IkBa_gene(state = 'on') + IkBa_mRNA(), c1a) #IkBa Transcription
Rule('IkBan_NFkBn_to_IkBac_NFkBc', NFkB(b = 1, loc = 'n') % IkBa(b = 1, phos = 'u', loc = 'n') >> NFkB(b = 1, loc = 'c') % IkBa(b = 1, phos = 'u', loc = 'c'), e2a) #IkBa\NFkB from nuc to cyto
Rule('IkBan_and_A20_on', IkBa(b = None, phos = 'u', loc = 'n') + A20_gene(state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n') + A20_gene(state = 'off'), q2) #turning off A20 gene state
Rule('IkBan_and_IkBa_on', IkBa(b = None, phos = 'u', loc = 'n') + IkBa_gene(state = 'on') >> IkBa(b = None, phos = 'u', loc = 'n') + IkBa_gene(state = 'off'), q2) #turning off IkBa gene state
Rule('IKKii_to_IKKn', IKK(state = 'ii') >> IKK(state = 'n'), k4) #IKKii creates IKKn
Rule('NFkBn_and_A20_off', NFkB(b = None, loc = 'n') + A20_gene(state = 'off') >> NFkB(b = None, loc = 'n') + A20_gene(state = 'on'), q1) #NFkB turning off A20 gene state
Rule('NFkBn_and_IkBa_off', NFkB(b = None, loc = 'n') + IkBa_gene(state = 'off') >> NFkB(b = None, loc = 'n') + IkBa_gene(state = 'on'), q1) #NFkB turning off IkBa gene state

# if __name__ == "__main__":
generate_equations(model)

time = np.linspace(0, 1800, 1801)
sim = BngSimulator(model, time) #integrator='lsoda',
result = sim.run(method='ode')

plt.figure()
plt.plot(time, result.observables['NFkBn_obs'], label='NFkBn_obs')
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.legend(loc = 'best')

plt.figure()
plt.plot(time, result.observables['IkBa_on_obs'], label='Ikba_gene')
plt.xlabel("Time (in minutes)", fontsize=16)
plt.ylabel("Concentration", fontsize=16)
plt.legend(loc = 'best')

#
# plt.figure(5)
# plt.plot(time/60, x["IkBa_off_obs"], label=IkBa_off_obs)
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend()
#
# plt.figure(6)
# plt.plot(time/60, x["IkBa_on_obs"], label=IkBa_on_obs)
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend()

plt.show()

    # for obs in ["IkBap_NFkBc_obs"]:
    #     plt.figure(1)
    #     plt.plot(time/60, x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     # plt.subplot(2,1,1)
    #     # plt.plot(time/60, x["TNF_ext_obs"], label = 'TNF')
    #     # plt.subplot(2,1,2)
    #     # plt.plot(time/60, x["IKKKa_obs"], label = 'IKKKa')
    #     # plt.figure(2)
    #     # plt.subplot(2,1,1)
    #     # plt.plot(time/60, x["IKKii_obs"], label = 'IKKii')
    #     # plt.subplot(2,1,2)
    #     # plt.plot(time/60, x["IKKn_obs"], label = 'IKKn')
    #     # plt.figure(3)
    #     # plt.subplot(2,1,1)
    #     # plt.plot(time/60, x["IKKa_obs"], label = 'IKKa')
    #     # plt.subplot(2,1,2)
    #     # plt.plot(time/60, x["IKKi_obs"], label = 'IKKi')
    #     # plt.figure(4)
    #     # plt.subplot(2,1,1)
    #     # plt.plot(time/60, x["IkBap_obs"], label = 'IkBap')
    #     # plt.subplot(2,1,2)
    #     # plt.plot(time/60, x["IkBap_NFkBc_obs"], label = 'IkBap_NFkBc')
    #     plt.legend(loc=0, prop={'size': 16})
    #     plt.xlabel("Time (in minutes)", fontsize=16)
    #     plt.ylabel("Concentrations", fontsize=16)
    #     plt.xlim(xmax = 300)
    #     plt.ylim(ymin = 0, ymax = 1000)
    #
    #
    # plt.figure(2)
    # for obs in ["IKKKa_obs", "IKKii_obs", "IKKa_obs", "IKKi_obs", "IkBap_obs", "IkBan_NFkBn_obs"]:
    #     # plt.subplot(3,2,1)
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     # re.match(r"(\w+)_obs", obs).group()
    #     #plt.plot(time/60, x[obs.name], label=obs.name)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # # plt.figure()
    # # for i in range(len(model.species)):
    # #     plt.plot(time/60, x["__s%d" %i])
    # plt.figure(3)
    # for obs in ["TNFR1a_obs", "A20_on_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # plt.figure(4)
    # for obs in ["A20_obs", "A20_off_obs", "IkBac_obs", "IkBan_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # plt.figure(5)
    # for obs in ["A20_on_obs", "A20_off_obs", "IkBa_off_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)

    # plt.figure(6)
    # for obs in ["A20t_obs"]: #, "IkBat_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # plt.figure(7)
    # for obs in ["IkBat_obs"]: #, "IkBat_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # plt.figure(8)
    # for obs in ["NFkBn_obs", "NFkBc_obs"]: #, "IkBat_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # plt.figure(9)
    # for obs in ["IkBac_obs"]: #, "IkBat_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)

    # plt.figure(10)
    # for obs in ["A20_obs"]: #, "IkBat_obs"]:
    #     plt.plot(time/60., x[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
    #     plt.legend(loc=0, prop={'size': 16})
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentrations", fontsize=16)
    # plt.xticks(fontsize=10)
    # plt.yticks(fontsize=10)
    #
    # plt.show()


    # x = odesolve(model, time, verbose=True)
    # plt.figure(4)
    # plt.plot(time/60, x["NFkBn_obs"], label='NFkBn')
    # plt.xlabel("Time (in minutes)", fontsize=16)
    # plt.ylabel("Concentration", fontsize=16)
    # plt.legend()
    # #
    # # plt.figure(5)
    # # plt.plot(time/60, x["IkBa_off_obs"], label=IkBa_off_obs)
    # # plt.xlabel("Time (in minutes)", fontsize=16)
    # # plt.ylabel("Concentration", fontsize=16)
    # # plt.legend()
    # #
    # # plt.figure(6)
    # # plt.plot(time/60, x["IkBa_on_obs"], label=IkBa_on_obs)
    # # plt.xlabel("Time (in minutes)", fontsize=16)
    # # plt.ylabel("Concentration", fontsize=16)
    # # plt.legend()
    #
    # plt.show()

    # generate_equations(model, verbose = True)
    # for i,ode in enumerate(model.odes):
    #     print i,":",ode



    # yobs = pd.DataFrame(array,columns=[str(i) for i in model.observables])
# yobs = pd.DataFrame(array,columns=[str(i) for i in model.observables])
# print yobs
# yobs = pd.DataFrame(array,columns=['observable'+str(a) for a in range(cols)]) #creating panda data frame for labeling of array

# print (yobs.shape)
# yobs = pd.DataFrame(index = nsims, columns = len(model.observables))
# model.parameters['IKK2_0'].value = pchip_obj1[0]
    # last_conc[1] = ikk2[i+1] #updating IKK2 values using matlab IKK2 array values (pchip interpolation)
    # last_conc[1] = pchip_obj1[i] #updating IKK2 interpolated values over each simulation
    # if i == 0:
    #     array[0,:] = solver.y[0]
    # array[i+1 ,:] = solver.y[1] #taking each simulation index and all obs and keeping the last obs in last time point of simulation
        # array[0,:] = solver.yobs_view[0]
    # array[i+1 ,:] = solver.yobs_view[1] #taking each simulation index and all obs and keeping the last obs in last time point of simulation
    # print(solver.yobs_view[1])
    # ikk2_vals.append(solver.yobs_view[0][1])



# plt.figure('IkBa_P2')
# for obs in ["IkBa_mRNA_obs"]:
#     plt.plot(tspan/60., array[:,4], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure('IkBb_P2')
# for obs in ["IkBb_mRNA_obs"]:
#     plt.plot(tspan/60., array[:,5], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure('IkBe_P2')
# for obs in ["IkBe_mRNA_obs"]:
#     plt.plot(tspan/60., array[:,6], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure('IkBd_P2')
# for obs in ["IkBd_mRNA_obs"]:
#     plt.plot(tspan/60., array[:,7], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
# print(ikk2_vals)
# print(myarray)

# plt.figure('IkB_P2')
# for obs in ["IkBa_obs", "IkBb_obs", "IkBe_obs", "IkBd_obs"]:
#     plt.plot(tspan/60., array[:,[]], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#




# for i,ode in enumerate(model.odes):
#      print i,":",ode

# for i,obs in enumerate(model.observables):
#      print i,":",obs
#
# for i,sp in enumerate(model.species):
#      print i,":",sp
# #
# # for i,reactions in enumerate(model.reactions):
# #      print i,":",reactions
# #
# for i,parameters in enumerate(model.parameters):
#      print i,":",parameters
# if __name__ == "__main__":
    # tspan = np.linspace(0, 720, 721)
    # # print len(tspan)
    # solver = Solver(model,tspan,verbose=True)
    # generate_equations(model, verbose = True)


# plt.figure()
# plt.plot( time_equil, equil["NFkBn_free"], label = NFkBn_free.name)
# plt.legend(loc = 0)
#
# # time = np.linspace(0,720, 721)
# nsims = len(tspan) - 1 #simulation time is 720
# # phase2 = odesolve(model, time, verbose=True)
#
# last_conc = [equil[x][-1] for x in ['__s%d' %i for i in np.arange(len(model.species))]]
# assert model.observables[1].name == 'IKK2_obs'
# assert str(model.species[1]) == "IKK2(ikb=None, S='C')"
# for i in range(nsims):
#     print(i)
#     model.parameters['IKK2_0'].value = myarray[i]
#     solver.run(y0 = last_conc) #setting initial concentrations to previous simulation value
#     last_conc = solver.y[1,:] # running solver taking first row and all simulation observables
#     last_conc[1] = myarray[i]
# print(last_conc)
    # print(i) #printing each simulation iteration



# plt.figure()
# plt.plot(time, phase2["IKK2_obs"], label = IKK2_obs.name)
# plt.legend(loc = 0)
#
# plt.figure()
# plt.plot(time, phase2["NFkBn_free"], label = NFkBn_free.name)
# plt.legend(loc = 0)
#

# plt.figure('IkBa_equil')
# for obs in ["IkBa_mRNA_obs"]:
#     plt.plot(time_equil, equil[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure('IkBb_equil')
# for obs in ["IkBb_mRNA_obs"]:
#     plt.plot(time_equil, equil[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure('IkBe_equil')
# for obs in ["IkBe_mRNA_obs"]:
#     plt.plot(time_equil, equil[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
# plt.figure('IkBd_equil')
# for obs in ["IkBd_mRNA_obs"]:
#     plt.plot(time_equil, equil[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)

#
# plt.figure('IkB_P2')
# for obs in ["IkBa_mRNA_obs", "IkBb_mRNA_obs", "IkBe_mRNA_obs", "IkBd_mRNA_obs"]:
#     plt.plot(time/60., phase2[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)

# plt.figure('IkB_total')
# for obs in ["IkBa_obs", "IkBb_obs", "IkBe_obs", "IkBd_obs"]:
#     plt.plot(time/60., phase2[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
#
#
# plt.figure('IkB_total_equil')
# for obs in ["IkBa_obs", "IkBb_obs", "IkBe_obs", "IkBd_obs"]:
#     plt.plot(np.linspace(-4000, 1, 4001), equil[obs], label=re.match(r"(\w+)_obs", obs).group(), linewidth=3)
#     plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time (in minutes)", fontsize=16)
# plt.ylabel("Concentrations", fontsize=16)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)

# plt.figure(2)
# plt.plot( np.linspace(-4000, 1, 4001), equil["IKK2_obs"], label = IKK2_obs.name)
# plt.legend(loc = 0)



# x = np.array([0, 5, 10, 15, 20, 25, 30, 260, 480, 2880], dtype = float)
# y1 = np.array([0.01, 0.6, 1, 0.65, 0.5, 0.36, 0.2, 0.1, 0.03, 0.01], dtype = float)
# # pchip_obj1 = interp1d(x,y1, kind='cubic')(tspan)
# pchip_obj1 = scipy.interpolate.PchipInterpolator(x, y1)(tspan)
# print(pchip_obj1)
# print(result)
# quit()
# exit()
# xi = np.arange(x[0],x[-1],0.01)
# yi1 = pchip_obj1(xi)
#
# plt.figure(1)
# plt.plot(tspan/60., array[:,2], label= NFkBn_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)
#
# plt.figure(1)
# plt.plot(tspan/60., array[:,3], label= NFkBn_bound.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)

# plt.figure()
# plt.plot(tspan/60., array[:,6], label= IkBa_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)
#
#
# plt.figure()
# plt.plot(tspan/60., array[:,1], label= IKK2_obs.name)
# plt.xlabel("Time (in hours)", fontsize=16)
# plt.ylabel("Concentration", fontsize=16)
# plt.legend(loc=0)


# print str(yi1)[0:-1]
# plt.plot(xi,yi1)
# plt.show()

# if __name__ == "__main__":
#     solver = Solver(model,time,verbose=True)
#     generate_equations(model, verbose = True)
#
#     # time = np.linspace(0, 4000 + 12*60, 4721)
#     time = np.linspace(0, 600, 601)
#     x = solver.run()
#     # x = odesolve(model, time, verbose=True) #integrator='lsoda'
#
#

#     # plt.ylim(0, .125)
#
#     for i in x["NFkBn_free"]:
#         print i
#     print
#     for i in x["NFkBn_bound"]:
#         print i
    # quit()
#     plt.figure(2)
#     plt.plot(time/60, x["IkBa_mRNA_obs"], label= IkBa_mRNA_obs)
#     plt.xlabel("Time (in minutes)", fontsize=16)
#     plt.ylabel("Concentration", fontsize=16)
#     plt.legend()
# #
# #     plt.figure(2)
# #     plt.plot(time/60, x["IkBa_obs"], label=IkBa_obs)
# #     plt.xlabel("Time (in minutes)", fontsize=16)
# #     plt.ylabel("Concentration", fontsize=16)
# #     plt.legend()
# # #
# plt.show()