from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
from pysb.simulator.bng import BngSimulator

Model()

# model.enable_synth_deg()

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
# 4.278903092658225660e+00,1.0
# ]

#fst-pso best fit so far 300
fstpso = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 10000,
4.297407813099238380e-05,
4.562857778542710359e-03,
4.083583806088330465e-02,
3.268921932083641298e-05,
4.420584558572939718e-04,
1.485042551633140245e-05,
1.826361893397378602e-02,
3.507084838987814003e-05,
4.170264612462448945e-04,
8.959852502456857637e-07,
3.776291959790734065e-02,
1.314285237600542500e+00,
1.093149358778181988e-05,
3.403706736201566121e-03,
2.140740462939950052e-06,
1.929684455430751358e-03,
4.911556238775962036e-01,
1.942601055868577861e-05,
2.600881539371887458e-02,
2.025760859508057621e+00,
1.807053446596232626e-01,
9.630081674576946899e-06,
1.087005537066853601e-04,
5.560331800852477691e-02,
4.672696328534567471e-05,
8.792909418558267354e-02,
1.089562750437665477e-02,
1.331030111195359855e-01,
3.199511015473860409e-01,
2.455343173071242920e-05,
3.665659601263383804e-02,
2.930188791118327085e+00,
2.079368801801380340e-01,
3.026297571627492231e-02,
4.498574912060039448e-03,
6.729428782950642980e-07,
4.485287992136575974e-02,
1.0
]

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
Parameter('TRADD_0', 9000)
Parameter('RIP1_0', 40000)
Parameter('TRAF_0', 9000)
Parameter('cIAP_0', 9000)
Parameter('A20_0', 9000)
Parameter('CYLD_0', 9000)
Parameter('FADD_0', 8030)
Parameter('flip_L_0', 3900)
Parameter('Lubac_0', 7226)
Parameter('C8_0', 9000)
Parameter('RIP3_0', 40000)
Parameter('MLKLa_0', 10000)


Initial(TNF(brec=None), TNF_0)
Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
Initial(cIAP(btraf=None), cIAP_0)
Initial(A20(brip = None), A20_0)
Initial(CYLD(brip=None, btraf = None), CYLD_0)
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
Initial(flip_L(bDED=None, state = 'A'), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
Initial(C8(bf=None, flip = None, state='I'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)


#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
Parameter('bind_TNF_TNFR_kf', 1e-6)
Parameter('bind_TNF_TNFR_kr',1e-3)
Parameter('TNF_deg1', 0.001)
Parameter('bind_TNFRANY_TRADD_kf', 1e-6)
Parameter('bind_TNFRANY_TRADD_kr',1e-3)
Parameter('bind_TNFRANY_RIP1unmod_kf', 1e-6)
Parameter('bind_TNFRANY_RIP1unmod_kr', 1e-3)
Parameter('bind_RIP1ANY_TRAFunmod_kf', 1e-6)
Parameter('bind_RIP1ANY_TRAFunmod_kr', 1e-3)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf', 1e-6)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr',1e-3)



Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)

Rule('TNF_deg', TNF(brec = None) >> None, TNF_deg1)

Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <>
     TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)

Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)

Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)


Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)


Parameter('CompI_UB2', 1e-1)
Parameter('bind_LUBAC_kf', 1e-6)
Parameter('bind_LUBAC_kr', 1e-3)
Parameter('bind_RIP1K63ubANY_A20_kf', 1e-6)
Parameter('bind_RIP1K63ubANY_A20_kr',1e-3)
Parameter('k_A20_1', 1e-1)
Parameter('bind_RIP1K63ubANY_CYLD_kf', 1e-6)
Parameter('bind_RIP1K63ubANY_CYLD_kr', 1e-3)
Parameter('k_CYLD_1', 1e-1)

Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)


Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     bind_LUBAC_kf, bind_LUBAC_kr)


#RIP1 K63ub to be deub by A20 or CYLD

Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
     bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)


Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
     >>  TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None), k_A20_1)



Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
     bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)

Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
     >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     k_CYLD_1)

#Initiating Necroptosis
Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 1e-1)
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 3.11e-7)

Parameter('bind_FADD_proC8_2_kf', 3.27e-06)
Parameter('bind_FADD_proC8_2_kr', 0.018)

Parameter('bind_FADDANY_flip_L_kf',3.27e-06)
Parameter('bind_FADDANY_flip_L_kr', 0.018)

Parameter('bind_C8_flip_L_kf',3.27e-2)
Parameter('bind_C8_flip_L_kr', 0.018) #not used

Parameter('kc_c8_1', 1e-1)
Parameter('bind_FADDANY_RIP3_kf', 1e-6)
Parameter('bind_FADDANY_RIP3_kr', 1e-3)
Parameter('kc_c8_2', 1e-1)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)


Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, flip = None, state='I')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = None, state='I'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, flip = None,state='I') + flip_L(bDED=None, state = 'A')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A'), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)

Rule('activate_C8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A'), bind_C8_flip_L_kf,bind_C8_flip_L_kr)

Rule('catalyze_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A') >>
     TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,flip = None, state='I') + flip_L(bDED=None, state = 'I') , kc_c8_1)


#RIP3 reactions to MLKL

Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_RIP3_kf, bind_FADDANY_RIP3_kr)

Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod')
     >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)

Parameter('bind_RIP1_RIP3po4_kf', 1e-2)
Parameter('RIP1po4_RIP3po4_kf', 1e-3)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1e-3)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 1e-6)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 1)

Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_RIP1_RIP3po4_kf)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), RIP1po4_RIP3po4_kf)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
     <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)

generate_equations(model)

# print(model.parameters)

print(len(model.parameters))
print(len(fstpso))
# quit()

Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))

tnf = [2326, 233, 23, 2]
color = ['r', 'm', 'g', 'b']
lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']

tspan = np.linspace(0, 1440, 1441)
sim = BngSimulator(model, tspan=tspan)
result = sim.run(method='ode', param_values=fstpso, initials={TNF(brec=None): tnf})
df = result.dataframe

mlkl = [0, 170, 900, 4880, 9940, 10000]
x = [0, 60, 120, 240, 360, 480]

plt.figure()
for n in range(0,4):
    # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
    plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], '--', c = color[n],lw = 1.5) #fppf
plt.scatter(x, mlkl, color = 'tab:gray', marker = 's')
plt.xlabel('Time [minutes]', fontsize=16)
plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FP', loc=0, fontsize = 5)
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
plt.show()



# sim = BngSimulator(model, tspan=tspan)
# result = sim.run(method='ode', param_values=fstpso)
#
# mlkl = [0, 170, 900, 4880, 9940, 10000]
# x = [0, 60, 120, 240, 360, 480]
#
# plt.figure()
# plt.plot(tspan, result.observables['MLKLa_obs'][:])
# plt.scatter(x, mlkl, color = 'r')
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.ylim(ymax = 10000)
# plt.show()
