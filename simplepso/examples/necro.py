from pylab import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
from pysb.simulator import CupSodaSimulator
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator

Model()

# par_list = np.load('allchains5000.npy')
# # print(len(par_list))
# # quit()
# rate_params = model.parameters_rules()
# param_values = np.array([p.value for p in model.parameters])
# rate_mask = np.array([p in rate_params for p in model.parameters])
# unmasked_pars = [0] * 25000
# # ode_list = []
#
# for i,par in enumerate(par_list):
#     param_values[rate_mask] = 10 ** np.asarray(par)
#     unmasked_pars[i] = param_values
# unmasked_pars = np.array(unmasked_pars)
# print(unmasked_pars)

# fstpso = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 10000,
# 3.304257485026848768e-05,
# 9.791215971368776028e-03,
# 6.110068937548310437e-03,
# 4.319219461882335177e-05,
# 4.212644667502709987e-03,
# 1.164332269398710173e-05,
# 2.404257105292715788e-02,
# 3.311085510054400207e-05,
# 4.280399320119887552e-02,
# 2.645814637693955063e-05,
# 1.437707485722601597e-02,
# 2.303744002126363599e-01,
# 2.980688423948379739e-05,
# 4.879773212139151134e-02,
# 1.121503480627013705e-05,
# 1.866712857727401229e-03,
# 7.572177664736867708e-01,
# 1.591282730353343167e-05,
# 3.897146248650381478e-02,
# 3.076363174012029411e+00,
# 3.734859661130401243e+00,
# 3.216200470463125876e-06,
# 8.782429548444865440e-05,
# 2.906341314225960662e-02,
# 5.663104188870970508e-05,
# 2.110469405515222677e-02,
# 1.294086380531199176e-01,
# 3.127598126760888220e-01,
# 4.298489868360909627e-01,
# 2.332910188537793332e-06,
# 7.077504621536276526e-03,
# 6.294061533948092091e-01,
# 6.419313355304218094e-02,
# 8.584653667640911989e-04,
# 8.160445062706172072e-05,
# 4.354383618147421691e-06,
# 4.278903092658225660e+00
# ]


Monomer('TNF', ['brec'])
Monomer('TNFR', ['blig', 'brip', 'bDD'])
Monomer('TRADD', ['brec', 'brip', 'state','bDD1', 'bDD2'], {'state': ['unmod', 'K63ub']})
Monomer('RIP1', ['bscf', 'bub1', 'bub2', 'bub3','bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state': ['unmod', 'K63ub', 'deub','ub', 'po4', 'trunc']})
Monomer('TRAF', ['brip', 'bciap', 'bcyld', 'state'], {'state': ['unmod', 'K63ub']})
Monomer('cIAP', ['btraf'])
Monomer('A20', ['brip'])
Monomer('A20t')
Monomer('CYLD', ['brip','btraf'])
Monomer('FADD', ['bDD', 'bDED1', 'bDED2'])
Monomer('proC8', ['bDED'])
Monomer('C8', ['bf', 'flip', 'state'], {'state': ['A', 'I']})
Monomer('flip_L', ['bDED', 'state'], {'state': ['A', 'I']})
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('LUBAC', ['brip'])

Parameter('TNF_0', 2326)
Parameter('TNFR_0', 4800)
Parameter('TRADD_0', 4696) #4696
Parameter('RIP1_0', 40000)
Parameter('TRAF_0', 11776) #11776
Parameter('cIAP_0', 9000) #9000
Parameter('A20_0', 9000)
Parameter('CYLD_0', 9000)
Parameter('FADD_0', 3109) #3109
Parameter('flip_L_0', 3900) #3900
Parameter('Lubac_0', 7226)
Parameter('C8_0', 3799) #3799
Parameter('RIP3_0', 10654) #10654
Parameter('MLKLa_0', 25544) #5544


Initial(TNF(brec=None), TNF_0)
Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
Initial(cIAP(btraf=None), cIAP_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)
Initial(A20(brip = None), A20_0)
Initial(CYLD(brip=None, btraf = None), CYLD_0)
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
Initial(flip_L(bDED=None, state = 'A'), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
Initial(C8(bf=None, flip = None, state='I'), C8_0)


#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
Parameter('p1f', 3.304257e-05)
Parameter('p1r',9.791216e-03)
Parameter('p2f', 6.110069e-03)
Parameter('p3f', 4.319219e-05)
Parameter('p3r',4.212645e-03)
Parameter('p4f', 1.164332e-05)
Parameter('p4r', 2.404257e-02)
Parameter('p5f', 3.311086e-05)
Parameter('p5r', 4.280399e-02)
Parameter('p6f', 2.645815e-05)
Parameter('p6r',1.437707e-02)

Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) | TNF(brec=1) % TNFR(blig=1, brip=None), p1f, p1r)

Rule('TNF_deg', TNF(brec = None) >> None, p2f)

Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) |
     TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), p3f,p3r)

Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
     | TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), p4f,p4r)

Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
     | TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), p5f,p5r)


Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
     | TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), p6f,p6r)


Parameter('p7f', 2.303744e-01)
Parameter('p8f', 2.980688e-05)
Parameter('p8r', 4.879773e-02)
Parameter('p9f', 1.121503e-05)
Parameter('p9r',1.866713e-03)
Parameter('p10f', 7.572178e-01)
Parameter('p11f', 1.591283e-05)
Parameter('p11r', 3.897146e-02)
Parameter('p12f', 3.076363e+00)

Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), p7f)


Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     | TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     p8f,p8r)


#RIP1 K63ub to be deub by A20 or CYLD

Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
     | TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
     p9f,p9r)


Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
     >>  TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None), p10f)



Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
     | TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
     p11f,p11r)

Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
     >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     p12f)

#Initiating Necroptosis
Parameter('p13f', 3.734860e+00)
Parameter('p13r', 3.216200e-06)
Parameter('p14f', 8.782430e-05)
Parameter('p14r', 2.906341e-02)
Parameter('p15f', 5.663104e-05)
Parameter('p15r', 2.110469e-02)
Parameter('p16f', 1.294086e-01)
Parameter('p16r', 3.127598e-01)

# Parameter('p17f', 4.298490e-01)
Parameter('p18f', 2.332910e-06)
Parameter('p18r', 7.077505e-03)
Parameter('p19f', 6.294062e-01)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     | TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), p13f, p13r)


Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, flip = None, state='I')
     | TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = None, state='I'), p14f,p14r)


Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, flip = None,state='I') + flip_L(bDED=None, state = 'A')
     | TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A'), p15f,p15r)

Rule('activate_C8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A')
     | TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A'), p16f,p16r)

# Rule('catalyze_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A') >>
#      TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,flip = None, state='I') + flip_L(bDED=None, state = 'I') , p17f)

#@TODOAddedReactions
#added params
Parameter('p17f', 4.298490e-01)
Parameter('p24f', 2.332910e-06)
Parameter('p24r', 7.077505e-03)
Parameter('p25f', 6.294062e-01)

Rule('catalyze_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A') >>
     TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,flip = 4, state='A') % flip_L(bDED=4, state = 'A') , p17f)
#
#
Rule('bind_C8Flip_RIP1RIP3', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod') +
     C8(bf=None,flip = 4, state='A') % flip_L(bDED=4, state = 'A') | TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=6,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None)
     % RIP3(bRHIM=5, bDD = None, state='unmod')% C8(bf=6,flip = 4, state='A') % flip_L(bDED=4, state = 'A'), p24f, p24r)

Rule('RIP1truc_C8Flip', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=6,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None)
     % RIP3(bRHIM=5, bDD = None, state='unmod')% C8(bf=6,flip = 4, state='A') % flip_L(bDED=4, state = 'A') >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') +
     FADD(bDD=None,bDED1 = None, bDED2 = None) +  RIP3(bRHIM=None, bDD = None, state='unmod') + C8(bf=None,flip = 4, state='A') % flip_L(bDED=4, state = 'A'), p25f)
# END ADDED REACTIONS





#RIP3 reactions to MLKL

Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
     | TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), p18f,p18r)

Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod')
     >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), p19f)

Parameter('p20f', 6.419313e-02)
Parameter('p21f', 8.584654e-04)
Parameter('p22f', 8.160445e-05)
Parameter('p22r', 4.354384e-06)
Parameter('p23f', 4.278903e+00)

Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), p20f)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), p21f)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
     | RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod'), p22f,p22r)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') + RIP1(bscf = None, bub1 = None, bub2= None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') , p23f)


Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('CII_C8a_obs', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') %
           FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A'))
Observable('CII_RIP3_obs', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub')
           % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'))
Observable('RIP1RIP3unmod_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod'))
Observable('CI_k63_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub')
           % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
Observable('CII_RIP1deub_obs', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub')
           % FADD(bDD=1,bDED1 = None, bDED2 = None))
Observable('A20_CI_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub')
           % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7))
Observable('CYLD_CI_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub')
           % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7))
# Observable('C8_obs', )
Observable('C8Flip_obs', C8(bf=None,flip = 4, state='A') % flip_L(bDED=4, state = 'A'))
Observable('RIP1TRADD_obs',  TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub'))
Observable('RIP1k63_obs', RIP1(state='K63ub'))
Observable('RIP1deub_obs', RIP1(state='deub'))
Observable('RIP1trunc_obs', RIP1(state='trunc'))
Observable('Flip_obs', flip_L(bDED=None, state = 'A'))
generate_equations(model)

# #
# #
# # print('done')
# #
# # print(len(model.parameters_rules()))
# # quit()
# # print(len(model.species))
#
# # par_files = os.listdir('10000pso')[:500]
# # n_pars = len(par_files)
# # all_pars = np.zeros((n_pars, 51))
#
# # rate_params = model.parameters_rules() # these are only the parameters involved in the rules
# # param_values = np.array([p.value for p in model.parameters]) # these are all the parameters
# # rate_mask = np.array([p in rate_params for p in model.parameters])
#
# # for i in range(n_pars):
# #     par = np.load('10000pso/{0}'.format(par_files[i]))
# #     param_values[rate_mask] = 10 ** par
# #     all_pars[i] = param_values
# # obs_y_range = ['MLKLa_obs']
# # print(len(all_pars))
#
# # tspan = np.linspace(0, 1440, 1441)
# # result = ScipyOdeSimulator(model, tspan=tspan).run(param_values=all_pars[:])
# # df = result.dataframe
# # print(df)
# # quit()
# #
# #
# # pars = np.array([
# # -4.         ,-4.56486775 ,-1.98677211 ,-4.57835751 ,-1.93305594 ,-6.94977109  ,
# #  -1.55078776, -6.58662292, -1.9162339 , -4.50644771, -2.69849836, -0.93036868 ,
# #  -4.64730879, -3.8032061 , -5.00570374, -2.55052241, -1.10643895, -5.241637   ,
# #  -3.38531962, -2.73739373, -0.17986558, -8.34729917, -4.28342644, -2.81749845 ,
# #  -4.96165901, -3.12204622,  0.51454775, -1.0972641 , -1.53217513, -5.16454949 ,
# #  -4.16652382, -0.81268578, -0.5795986 , -1.22351651, -4.49737028, -4.99775591 ,
# #   1.01013847])
# # # #
# # # # # #
# # # #
#
# pars = np.array([-2.48093845e+00, -1.54667397e-01, -2.24894369e+00 -6.36459478e+00 ,
# #  -6.10173330e-01, -3.56792093e+00 ,-2.43588076e+00,-5.91784104e+00                  ,
# #  -1.74601460e+00, -6.06447883e+00 , 6.26210176e-03,-1.48488691e-01                  ,
# #  -5.87884291e+00, -2.91476836e+00 ,-6.95019956e+00,-2.77208784e+00                  ,
# #  -1.21736663e+00, -6.79823713e+00 ,-2.34519663e+00,-1.51186077e+00                  ,
# #  -6.72442508e-02, -6.94976629e+00 ,-5.08751300e+00,-3.03126227e+00                  ,
# #  -3.95838232e+00, -2.25731733e+00 ,-1.12673197e+00,-1.73810601e+00                  ,
# #  -3.69051414e-01, -4.82778547e+00 ,-4.15011982e+00,-2.20106898e+00                  ,
# #  -3.03698641e+00, -2.58057532e+00 ,-2.08828616e+00,-5.72086803e+00                  ,
# #   1.14955541e+00])
# # print(len(pars))
# # quit()
# #
# pars = np.load('optimizer_best_1000_14.npy')
# # # # # pars =  np.load('optimizer_best_1000_mar7.npy')
pars =  np.load('optimizer_best_75_50_100TNF.npy')
rate_params = model.parameters_rules() # these are only the parameters involved in the rules
param_values = np.array([p.value for p in model.parameters]) # these are all the parameters
rate_mask = np.array([p in rate_params for p in model.parameters])
param_values[rate_mask] = 10 ** pars
print(len(param_values))
print(param_values)
# quit()
# #
# # data10 = np.array([0.0096*5544, 0.048*5544, 0.178*5544, 0.287*5544, 0.497*5544, 0.547*5544, 0.770*5544, 0.808*5544, 0.953*5544, 1.0*5544])
# # stdev10 = np.array([.05, .02, .08, .11, .11, .12, .16, .09, .06, .01])
tspan = np.linspace(0, 1440, 1441)
solver1 = BngSimulator(model, tspan=tspan)
result = solver1.run(param_values=param_values)
# result = solver1.run(method ='ode', initials = {TNF(brec=None): [2326]},  param_values= param_values)
# # print(result.observables['MLKLa_obs'].max())
# # plt.figure()
# # plt.plot(tspan/60, result.observables['MLKLa_obs'])
# # plt.show()
# # quit()
df = result.dataframe
# # #with open('params_cal_necromlkl14neww.txt', 'w') as f:
# # #    for p, v in zip(model.parameters, result.param_values[0]):
# # #        f.write('{},{:e}\n'.format(p.name, v))
# # # t = np.array([0., 30,  60,   120,  180, 270,  480,  960, 1440])
# # # data = np.array([0., 0., 0., 0., 0.01*5544, 0.05*5544, 0.5*5544, 0.99*5544, 1.*5544])
# # #
# # #

# x100 = np.array([0., .5, 1.5, 4.5, 8, 12, 16])
# y100 = np.array([
# 0.00885691708746097*5544,0.0161886154261265*5544,
# 0.0373005242261882*5544,
# 0.0498939020159581*5544,
# 0.639729406776*5544,
# 1*5544])
# y10 = np.array([0.0106013664572332*5544,
# 0.00519576571714913*5544,
# 0.02967443048221*5544,
# 0.040022163974868*5544,
# 0.198128107774737*5544,
# 0.540551401148672*5544])

y10_long = np.array([0., 0.001, 0.02, 0.03, 0.04, 0.06, .09, .21, .40, .65, .81])
t = np.array([0, 2,4,6,8,10,12,14,16,18, 20])
# y1 = np.array([0.0060632926030143*5544,
# 0.00942691670200054*5544,
# 0.0113342231044983*5544,
# 0.0312716821493274*5544,
# 0.10956134602783*5544,
# 0.361058859808041*5544])
y0 = np.array([0.00777642331264696*5544,
0.00919829048298192*5544,
0.00177263616389778*5544,
0.00637961107042362*5544,
0.0136330127587403*5544,
0.0852337335697491*5544])

x1 = np.array([0, 0.5, 1.5, 4.5, 8., 12., 16.])
# y1 = np.array([0, 0.006, 0.009, 0.011, 0.0313, 0.110, 0.261])

x10 = np.array([0, .5, 1.5, 4.5, 8, 12, 16])
# y10 = np.array([0, 0.001, 0.02, .10, .15, .41, 0.70])

y10 = np.array([0., 0.0106013664572332,
0.00519576571714913,
0.02967443048221,
0.050022163974868,
0.198128107774737,
0.640551401148672])
y1 = np.array([0., 0.0060632926030143,
0.00942691670200054,
0.0113342231044983,
0.0312716821493274,
0.10956134602783,
0.36105885980804])
# y100 = np.array([0., 0.008, 0.016, 0.037, 0.079, 0.639, 1.0])



x100 = np.array([0., 0.5, 1.5, 4.5, 8., 12., 16.])
# y100 = np.array([0., 0.008, 0.016, 0.037, 0.079, 0.639, 1.0])
y100 = np.array([0., 0.01, 0.0152, 0.03, 0.337, 0.693, 0.99])
plt.figure()
# for n in range(0,2):

plt.scatter(x100,y100, label = '100 ng/ml')
plt.plot(x100, y100)
# plt.scatter(x10, y100, label = '10 ng/ml')
# plt.scatter(x1, y1, label ='1')
# plt.scatter(x100,y1/5544, label = '1 ng/ml')
# plt.scatter(x100,y0/5544, label = '0.1 ng/ml')
# plt.plot(tspan/60, df.loc[0]['MLKLa_obs'].iloc[:]/5544,lw =1.5, label = '100') #fst-pso
# plt.plot(tspan/60, df.loc[0]['MLKLa_obs'].iloc[:]/5544,lw =1.5, label = '10') #fst-pso
# plt.plot(tspan/60, df.loc[1]['MLKLa_obs'].iloc[:]/5544,lw =1.5, label = '1') #fst-pso
# plt.plot(tspan/60, df.loc[3]['MLKLa_obs'].iloc[:]/5544,lw =1.5, label = '01') #fst-pso
plt.plot(tspan/60, result.observables['MLKLa_obs'][:]/25544,lw = 1.5) #fppf
# plt.scatter(t/60, data, lw = 2.5, color = 'k',marker = '.', zorder = 1)
# plt.plot(t/60, data)
plt.xlabel('Time [hours]', fontsize=14)
plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=14)
plt.title('pMLKL to 100 ng/ml TNFa (param set 14)')
# plt.title('Sensitivity of pMLKL to varying TNFa doses (WTKD cal)')
# plt.legend([ 'FDKD', 'TDKD', 'A20KD', 'C8KD','WT '] , title = 'KD Conditions', loc='best', fontsize = 8)
plt.legend(loc = 'best', fontsize = 12)
# plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# plt.ylim(ymax = 6000)
# plt.savefig('pMLKL to paramset14')
plt.show()
quit()
