# import the pysb module and all its methods and functions
from pysb import *
from pysb.simulator import ScipyOdeSimulator
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import *
from pysb.bng import generate_equations, generate_network
from pysb.core import ComponentSet
from pysb.export import export

# instantiate a model
Model ()

#Monomers: Name, binding site (named according to what it binds to or if the site gets modified=mod), optional (the states the binding site can be in)
Monomer('TNF', ['tnfr'])
Monomer('TNFR', ['tnf', 'tradd', 'rip1'])
Monomer('TRADD', ['tnfr', 'rip1'])
#ub=ubiquitinated, db=deubiquitinated
Monomer('RIP1', ['tnfr', 'tradd', 'a20', 'fadd', 'rip3c8', 'mod'], {'mod': ['ub', 'db']})
Monomer('A20', ['rip1'])
Monomer('FADD', ['rip1'])
Monomer('RIP3', ['rip1'])
#i=inactive, a=active
Monomer('C8', ['rip1', 'flip', 'mod'], {'mod' : ['i', 'a']})
Monomer('FLIP', ['c8'])

#Parameters: Rates of reactions
Parameter('kf_TNF_TNFR_bind', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_bind', 1.00e-03) #1.00e-03
Parameter('kc_degrade_TNF', 1.00e-03) #1.00e-03
Parameter('kf_TNF_TNFR_bind_TRADD', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_bind_TRADD', 1.00e-03) #1.00e-03
Parameter('kf_TNF_TNFR_TRADD_bind_RIP1', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_TRADD_bind_RIP1', 1.00e-03) #1.00e-03
Parameter('kf_TNF_TNFR_TRADD_RIP1_bind_A20', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_TRADD_RIP1_bind_A20', 1.00e-03) #1.00e-03
Parameter('kc_TNF_TNFR_TRADD_RIP1_A20_disassociate', 1.00e-01) #1.00e-01
Parameter('kf_RIP1_TRADD_bind_FADD', 1.00e-01) #1.00e-01
Parameter('kr_RIP1_TRADD_bind_FADD', 3.11e-07) #3.11e-07
Parameter('kf_RIP1_TRADD_FADD_bind_RIP3', 3.27e-06) #same binding rate as C8
Parameter('kr_RIP1_TRADD_FADD_bind_RIP3', 0.018) #same binding rates as C8
Parameter('kf_RIP1_TRADD_FADD_bind_C8', 3.27e-06) #3.27e-06
Parameter('kr_RIP1_TRADD_FADD_bind_C8', 0.018) #0.018
Parameter('kf_RIP1_TRADD_FADD_C8_bind_FLIP', 3.27e-06) #3.27e-06
Parameter('kr_RIP1_TRADD_FADD_C8_bind_FLIP', 0.018) #0.018
Parameter('kc_FLIP_activates_C8', 3.27e-02) #3.27e-02

#Rules     ****Is the last reaction labeled well?*****
#Unbound TNF and unbound TNFR reversibly bind to create the TNF/TNFR complex
Rule('TNF_TNFR_bind', TNF(tnfr=None) + TNFR(tnf=None, tradd=None, rip1=None) | TNF(tnfr=1) % TNFR(tnf=1, tradd=None, rip1=None), kf_TNF_TNFR_bind, kr_TNF_TNFR_bind)

#Unbound TNF is degraded
Rule('degrade_TNF', TNF(tnfr=None) >> None, kc_degrade_TNF)

#TNF/TNFR complex and unbound TRADD reversibly bind to create the TNF/TNFR/TRADD complex
Rule('TNF_TNFR_bind_TRADD', TNF(tnfr=1) % TNFR(tnf=1, tradd=None, rip1=None) + TRADD(tnfr=None, rip1=None) | TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=None) % TRADD(tnfr=2, rip1=None),
     kf_TNF_TNFR_bind_TRADD, kr_TNF_TNFR_bind_TRADD)

#TNF/TNFR/TRADD complex and RIP1 reversibly bind to create the TNF/TNFR/TRADD/RIP1 complex
Rule('TNF_TNFR_TRADD_bind_RIP1', TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=None) % TRADD(tnfr=2, rip1=None) + RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub')
     | TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=3) % TRADD(tnfr=2, rip1=4) % RIP1(tnfr=3, tradd=4, a20=None, fadd=None, rip3c8=None, mod='ub'),
     kf_TNF_TNFR_TRADD_bind_RIP1, kr_TNF_TNFR_TRADD_bind_RIP1)

#TNF/TNFR/TRADD/RIP1 complex and A20 reversibly bind to create the TNF/TNFR/TRADD/RIP1/A20 complex
Rule('TNF_TNFR_TRADD_RIP1_bind_A20', TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=3) % TRADD(tnfr=2, rip1=4) % RIP1(tnfr=3, tradd=4, a20=None, fadd=None, rip3c8=None, mod='ub') + A20(rip1=None)
    | TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=3) % TRADD(tnfr=2, rip1=4) % RIP1(tnfr=3, tradd=4, a20=5, fadd=None, rip3c8=None, mod='ub') % A20(rip1=5),
     kf_TNF_TNFR_TRADD_RIP1_bind_A20, kr_TNF_TNFR_TRADD_RIP1_bind_A20)

#TNF/TNFR/TRADD/RIP1/A20 complex nonreversibly disassociates leaving only the TRADD/RIP1 complex. RIP1 is deubiquitinated by A20 in the process.
Rule('TNF_TNFR_TRADD_RIP1_A20_disassociate', TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=3) % TRADD(tnfr=2, rip1=4) % RIP1(tnfr=3, tradd=4, a20=5, fadd=None, rip3c8=None, mod='ub') % A20(rip1=5)
     >> TNF(tnfr=None) + TNFR(tnf=None, tradd=None, rip1=None) + TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=None, rip3c8=None, mod='db') + A20(rip1=None),
     kc_TNF_TNFR_TRADD_RIP1_A20_disassociate)

#RIP1/TRADD and FADD reversibly bind to create the RIP1/TRADD/FADD complex
Rule('RIP1_TRADD_bind_FADD', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=None, rip3c8=None, mod='db') + FADD(rip1=None)
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=None, mod='db') % FADD(rip1=6),
     kf_RIP1_TRADD_bind_FADD, kr_RIP1_TRADD_bind_FADD)

#RIP1/TRADD/FADD complex and RIP3 reversibly bind to create the RIP1/TRADD/FADD/RIP3 complex
Rule('RIP1_TRADD_FADD_bind_RIP3', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=None, mod='db') % FADD(rip1=6) + RIP3(rip1=None)
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % RIP3(rip1=7),
     kf_RIP1_TRADD_FADD_bind_RIP3, kr_RIP1_TRADD_FADD_bind_RIP3)

#IRP1/TRADD/FADD complex and C8 reversibly bind to create the RIP1/TRADD/FADD/C8 complex
Rule('RIP1_TRADD_FADD_bind_C8', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=None, mod='db') % FADD(rip1=6) + C8(rip1=None, flip=None, mod='i')
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % C8(rip1=7, flip=None, mod='i'),
     kf_RIP1_TRADD_FADD_bind_C8, kr_RIP1_TRADD_FADD_bind_C8)

#RIP1/TRADD/FADD/C8 complex and FLIP reversibly bind to create the RIP1/TRADD/FADD/C8/FLIP complex
Rule('RIP1_TRADD_FADD_C8_bind_FLIP', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % C8(rip1=7, flip=None, mod='i') + FLIP(c8=None)
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % C8(rip1=7, flip=8, mod='i') % FLIP(c8=8),
     kf_RIP1_TRADD_FADD_C8_bind_FLIP, kr_RIP1_TRADD_FADD_C8_bind_FLIP)

#In the RIP1/TRADD/FADD/C8/FLIP complex, FLIP nonreversibly activates C8
Rule('FLIP_activates_C8', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % C8(rip1=7, flip=8, mod='i') % FLIP(c8=8)
     >> TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % C8(rip1=7, flip=8, mod='a') % FLIP(c8=8),
     kc_FLIP_activates_C8)

#Initial Conditions
Parameter('TNF_0', 2326) #698 is 30ng/ml of TNF
Parameter('TNFR_0', 4800) #0.00246
Parameter('TRADD_0', 9000) #
Parameter('RIP1_0', 40000) #47000 0.04
Parameter('A20_0', 9000) #2256
Parameter('FADD_0', 8030) # 0.0033
Parameter('RIP3_0', 40000) #20000
Parameter('C8_0', 9000) # 0.004 # 0.09
Parameter('FLIP_0', 3900) # 0.004 # 0.09
# Parameter('TNF_0', 2326) #698 is 30ng/ml of TNF
# Parameter('TNFR_0', 4800) #0.00246
# Parameter('TRADD_0', 9000) #
# Parameter('RIP1_0', 40000) #47000 0.04
# Parameter('A20_0', 9000) #2256
# Parameter('FADD_0', 8030) # 0.0033
# Parameter('RIP3_0', 40000) #20000
# Parameter('C8_0', 9000) # 0.004 # 0.09
# Parameter('FLIP_0', 3900) # 0.004 # 0.09
Initial(TNF(tnfr=None), TNF_0)
Initial(TNFR(tnf=None, tradd=None, rip1=None), TNFR_0)
Initial(TRADD(tnfr=None, rip1=None), TRADD_0)
Initial(RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub'), RIP1_0)
Initial(A20(rip1=None), A20_0)
Initial(FADD(rip1=None), FADD_0)
Initial(RIP3(rip1=None), RIP3_0)
Initial(C8(rip1=None, flip=None, mod='i'), C8_0)
Initial(FLIP(c8=None), FLIP_0)

#Observables
Observable('obsComplexI', TNF(tnfr=1) % TNFR(tnf=1, tradd=2, rip1=3) % TRADD(tnfr=2, rip1=4) % RIP1(tnfr=3, tradd=4, a20=None, fadd=None, rip3c8=None, mod='ub'))
Observable('obsComplexIIa', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % C8(rip1=7, flip=None, mod='i'))
Observable('obsComplexIIb', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=6, rip3c8=7, mod='db') % FADD(rip1=6) % RIP3(rip1=7))
# Observable('obsRIP1', RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub'))

print(len(model.rules()))
quit()

tspan = np.linspace(0, 1200, 1201) #put time of sim you want
sim = ScipyOdeSimulator(model, tspan=tspan)
sim_result = sim.run()

plt.figure(figsize=(6,3))
plt.subplot(3, 1, 1)
plt.plot(tspan/60, sim_result.observables['obsComplexI'],label = '100 ng/ml TNF')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("obsComplexI [Molecules/Cell]", fontsize=15)
plt.title('obsComplexI Trajectories')
plt.legend(loc ='best')

plt.figure(figsize=(6,3))
plt.subplot(3, 1, 2)
plt.plot(tspan/60, sim_result.observables['obsComplexIIa'],label = '100 ng/ml TNF')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("obsComplexIIa [Molecules/Cell]", fontsize=15)
plt.title('obsComplexIIa Trajectories')
plt.legend(loc ='best')

plt.figure(figsize=(6,3))
plt.subplot(3, 1, 3)
plt.plot(tspan/60, sim_result.observables['obsComplexIIb'],label = '100 ng/ml TNF')
plt.xlabel("Time (in hr)", fontsize=15)
plt.ylabel("obsComplexIIb [Molecules/Cell]", fontsize=15)
plt.title('obsComplexIIb Trajectories')
plt.legend(loc ='best')
plt.show()

generate_network(model)
generate_equations(model)
math_output = export(model, 'mathematica')
with open('mathematica', 'w') as f:
    f.write(math_output)

# plt.figure()
# plt.plot(tspan/60, sim_result.observables['obsComplexIIa'],label = '100 ng/ml TNF')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("obsComplexI [Molecules/Cell]", fontsize=15)
# plt.title('obsComplexI Trajectories')
# plt.legend(loc ='best')
# plt.show()