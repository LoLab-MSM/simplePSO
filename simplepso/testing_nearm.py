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
4.485287992136575974e-02]

#necroptosis monomers
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
Monomer('C8', ['bf', 'state'], {'state': ['pro', 'A']})
Monomer('flip', ['bf'])
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('LUBAC', ['brip'])

#apoptosis monomers
Monomer('BAR', ['bf'])
Monomer('Bid', ['bf', 'state'], {'state': ['U', 'T', 'M']})
Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state': ['C', 'M', 'A']})
Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state': ['M', 'A']})
Monomer('Bcl2', ['bf', 'state'], {'state': ['C', 'M']})
Monomer('BclxL', ['bf', 'state'], {'state': ['C', 'M']})
Monomer('Mcl1', ['bf', 'state'], {'state': ['C', 'M']})
Monomer('Bad', ['bf', 'state'], {'state': ['C', 'M']})
Monomer('Noxa', ['bf', 'state'], {'state': ['C', 'M']})
Monomer('CytoC', ['bf', 'state'], {'state': ['M', 'C', 'A']})
Monomer('Smac', ['bf', 'state'], {'state': ['M', 'C', 'A']})
Monomer('Apaf', ['bf', 'state'], {'state': ['I', 'A']})
Monomer('Apop', ['bf'])
Monomer('C3', ['bf', 'state'], {'state': ['pro', 'A', 'ub']})
Monomer('C6', ['bf', 'state'], {'state': ['pro', 'A']})
Monomer('C9', ['bf'])
Monomer('PARP', ['bf', 'state'], {'state': ['U', 'C']})
Monomer('XIAP', ['bf'])

#necroptosis initials
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
Parameter('C8_0', 20000)
Parameter('RIP3_0', 40000)
Parameter('MLKLa_0', 10000)

#apoptosis initials
Parameter('BAR_0', 1000.0)
Parameter('Apaf_0', 100000.0)
Parameter('C3_0', 10000.0)
Parameter('C6_0', 10000.0)
Parameter('C9_0', 100000.0)
Parameter('PARP_0', 1000000.0)
Parameter('XIAP_0', 100000.0)
Parameter('Bid_0', 40000.0)
Parameter('Bax_0', 80000.0)
Parameter('Bak_0', 20000.0)
Parameter('Bcl2_0', 20000.0)
Parameter('BclxL_0', 20000.0)
Parameter('Mcl1_0', 20000.0)
Parameter('Bad_0', 1000.0)
Parameter('Noxa_0', 1000.0)
Parameter('CytoC_0', 500000.0)
Parameter('Smac_0', 100000.0)

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
Initial(flip(bf = None), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
Initial(C8(bf=None, state='pro'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)

Initial(BAR(bf=None), BAR_0)
Initial(Apaf(bf=None, state='I'), Apaf_0)
Initial(C3(bf=None, state='pro'), C3_0)
Initial(C6(bf=None, state='pro'), C6_0)
Initial(C9(bf=None), C9_0)
Initial(PARP(bf=None, state='U'), PARP_0)
Initial(XIAP(bf=None), XIAP_0)
Initial(Bid(bf=None, state='U'), Bid_0)
Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
Initial(Bcl2(bf=None, state='M'), Bcl2_0)
Initial(BclxL(bf=None, state='M'), BclxL_0)
Initial(Mcl1(bf=None, state='M'), Mcl1_0)
Initial(Bad(bf=None, state='M'), Bad_0)
Initial(Noxa(bf=None, state='M'), Noxa_0)
Initial(CytoC(bf=None, state='M'), CytoC_0)
Initial(Smac(bf=None, state='M'), Smac_0)

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

Parameter('bind_proC8_C8a_kf', 3.27e-06)
Parameter('bind_proC8_C8a_kr', 0.018)

Parameter('bind_FADDANY_flip_L_kf',3.27e-06)
Parameter('bind_FADDANY_flip_L_kr', 0.018)

# Parameter('bind_C8_flip_L_kf',3.27e-2)
# Parameter('bind_C8_flip_L_kr', 0.018) #not used

# Parameter('kc_c8_1', 1e-1)
Parameter('bind_FADDANY_RIP3_kf', 1e-6)
Parameter('bind_FADDANY_RIP3_kr', 1e-3)
Parameter('kc_c8_2', 1e-1)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)


Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, state='pro')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='pro'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


Rule('C8_apop', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,state='pro')
     >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,state='A'), bind_proC8_C8a_kf, bind_proC8_C8a_kr)



Rule('bind_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + flip(bf=None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % flip(bf=2), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)


Rule('bind_flip_L_RIP3', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % flip(bf=2) + RIP3(bRHIM=None, bDD = None, state='unmod')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % flip(bf=2) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_RIP3_kf, bind_FADDANY_RIP3_kr)


# Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
     # <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_RIP3_kf, bind_FADDANY_RIP3_kr)


Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % flip(bf=2) % RIP3(bRHIM=5, bDD = None, state='unmod')
     >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + FADD(bDD=None,bDED1 = None, bDED2 = None) + flip(bf=None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)


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

#RULES FOR APOP
# Parameter('bind_L_R_to_LR_kf', 4e-07)
# Parameter('bind_L_R_to_LR_kr', 0.001)
# Parameter('convert_LR_to_DISC_kc', 1e-05)
# Parameter('bind_DISC_C8pro_to_DISCC8pro_kf', 1e-06)
# Parameter('bind_DISC_C8pro_to_DISCC8pro_kr', 0.001)
# Parameter('catalyze_DISCC8pro_to_DISC_C8A_kc', 1.0)
Parameter('bind_C8A_BidU_to_C8ABidU_kf', 1e-06)
Parameter('bind_C8A_BidU_to_C8ABidU_kr', 0.001)
Parameter('catalyze_C8ABidU_to_C8A_BidT_kc', 1.0)
# Parameter('bind_DISC_flip_kf', 1e-06)
# Parameter('bind_DISC_flip_kr', 0.001)
Parameter('bind_BAR_C8A_kf', 1e-06)
Parameter('bind_BAR_C8A_kr', 0.001)
# Parameter('Apaf_0', 100000.0)
# Parameter('C3_0', 10000.0)
# Parameter('C6_0', 10000.0)
# Parameter('C9_0', 100000.0)
# Parameter('XIAP_0', 100000.0)
# Parameter('PARP_0', 1000000.0)
#
# Rule('bind_L_R_to_LR', L(bf=None) + R(bf=None) <> L(bf=1) % R(bf=1), bind_L_R_to_LR_kf, bind_L_R_to_LR_kr)
# Rule('convert_LR_to_DISC', L(bf=1) % R(bf=1) >> DISC(bf=None), convert_LR_to_DISC_kc)
# Rule('bind_DISC_C8pro_to_DISCC8pro', DISC(bf=None) + C8(bf=None, state='pro') <>  DISC(bf=1) % C8(bf=1, state='pro'), bind_DISC_C8pro_to_DISCC8pro_kf, bind_DISC_C8pro_to_DISCC8pro_kr)
# Rule('catalyze_DISCC8pro_to_DISC_C8A', DISC(bf=1) % C8(bf=1, state='pro') >> DISC(bf=None) + C8(bf=None, state='A'), catalyze_DISCC8pro_to_DISC_C8A_kc)
Rule('bind_C8A_BidU_to_C8ABidU', C8(bf=None, state='A') + Bid(bf=None, state='U') <> C8(bf=1, state='A') % Bid(bf=1, state='U'), bind_C8A_BidU_to_C8ABidU_kf, bind_C8A_BidU_to_C8ABidU_kr)
Rule('catalyze_C8ABidU_to_C8A_BidT', C8(bf=1, state='A') % Bid(bf=1, state='U') >> C8(bf=None, state='A') + Bid(bf=None, state='T'), catalyze_C8ABidU_to_C8A_BidT_kc)
# Rule('bind_DISC_flip', DISC(bf=None) + flip(bf=None) <> DISC(bf=1) % flip(bf=1), bind_DISC_flip_kf, bind_DISC_flip_kr)
Rule('bind_BAR_C8A', BAR(bf=None) + C8(bf=None, state='A') <> BAR(bf=1) % C8(bf=1, state='A'), bind_BAR_C8A_kf, bind_BAR_C8A_kr)


Parameter('equilibrate_SmacC_to_SmacA_kf', 0.01)
Parameter('equilibrate_SmacC_to_SmacA_kr', 0.01)
Parameter('equilibrate_CytoCC_to_CytoCA_kf', 0.01)
Parameter('equilibrate_CytoCC_to_CytoCA_kr', 0.01)
Parameter('bind_CytoCA_ApafI_to_CytoCAApafI_kf', 5e-07)
Parameter('bind_CytoCA_ApafI_to_CytoCAApafI_kr', 0.001)
Parameter('catalyze_CytoCAApafI_to_CytoCA_ApafA_kc', 1.0)
Parameter('convert_ApafA_C9_to_Apop_kf', 5e-08)
Parameter('convert_ApafA_C9_to_Apop_kr', 0.001)
Parameter('bind_Apop_C3pro_to_ApopC3pro_kf', 5e-09)
Parameter('bind_Apop_C3pro_to_ApopC3pro_kr', 0.001)
Parameter('catalyze_ApopC3pro_to_Apop_C3A_kc', 1.0)
Parameter('bind_Apop_XIAP_kf', 2e-06)
Parameter('bind_Apop_XIAP_kr', 0.001)
Parameter('bind_SmacA_XIAP_kf', 7e-06)
Parameter('bind_SmacA_XIAP_kr', 0.001)
Rule('equilibrate_SmacC_to_SmacA', Smac(bf=None, state='C') <> Smac(bf=None, state='A'), equilibrate_SmacC_to_SmacA_kf, equilibrate_SmacC_to_SmacA_kr)
Rule('equilibrate_CytoCC_to_CytoCA', CytoC(bf=None, state='C') <> CytoC(bf=None, state='A'), equilibrate_CytoCC_to_CytoCA_kf, equilibrate_CytoCC_to_CytoCA_kr)
Rule('bind_CytoCA_ApafI_to_CytoCAApafI', CytoC(bf=None, state='A') + Apaf(bf=None, state='I') <> CytoC(bf=1, state='A') % Apaf(bf=1, state='I'), bind_CytoCA_ApafI_to_CytoCAApafI_kf, bind_CytoCA_ApafI_to_CytoCAApafI_kr)
Rule('catalyze_CytoCAApafI_to_CytoCA_ApafA', CytoC(bf=1, state='A') % Apaf(bf=1, state='I') >> CytoC(bf=None, state='A') + Apaf(bf=None, state='A'), catalyze_CytoCAApafI_to_CytoCA_ApafA_kc)
Rule('convert_ApafA_C9_to_Apop', Apaf(bf=None, state='A') + C9(bf=None) <> Apop(bf=None), convert_ApafA_C9_to_Apop_kf, convert_ApafA_C9_to_Apop_kr)
Rule('bind_Apop_C3pro_to_ApopC3pro', Apop(bf=None) + C3(bf=None, state='pro') <> Apop(bf=1) % C3(bf=1, state='pro'), bind_Apop_C3pro_to_ApopC3pro_kf, bind_Apop_C3pro_to_ApopC3pro_kr)
Rule('catalyze_ApopC3pro_to_Apop_C3A', Apop(bf=1) % C3(bf=1, state='pro') >> Apop(bf=None) + C3(bf=None, state='A'), catalyze_ApopC3pro_to_Apop_C3A_kc)
Rule('bind_Apop_XIAP', Apop(bf=None) + XIAP(bf=None) <> Apop(bf=1) % XIAP(bf=1), bind_Apop_XIAP_kf, bind_Apop_XIAP_kr)
Rule('bind_SmacA_XIAP', Smac(bf=None, state='A') + XIAP(bf=None) <> Smac(bf=1, state='A') % XIAP(bf=1), bind_SmacA_XIAP_kf, bind_SmacA_XIAP_kr)

Parameter('bind_C8A_C3pro_to_C8AC3pro_kf', 1e-07)
Parameter('bind_C8A_C3pro_to_C8AC3pro_kr', 0.001)
Parameter('catalyze_C8AC3pro_to_C8A_C3A_kc', 1.0)
Parameter('bind_XIAP_C3A_to_XIAPC3A_kf', 2e-06)
Parameter('bind_XIAP_C3A_to_XIAPC3A_kr', 0.001)
Parameter('catalyze_XIAPC3A_to_XIAP_C3ub_kc', 0.1)
Parameter('bind_C3A_PARPU_to_C3APARPU_kf', 1e-06)
Parameter('bind_C3A_PARPU_to_C3APARPU_kr', 0.01)
Parameter('catalyze_C3APARPU_to_C3A_PARPC_kc', 1.0)
Parameter('bind_C3A_C6pro_to_C3AC6pro_kf', 1e-06)
Parameter('bind_C3A_C6pro_to_C3AC6pro_kr', 0.001)
Parameter('catalyze_C3AC6pro_to_C3A_C6A_kc', 1.0)
Parameter('bind_C6A_C8pro_to_C6AC8pro_kf', 3e-08)
Parameter('bind_C6A_C8pro_to_C6AC8pro_kr', 0.001)
Parameter('catalyze_C6AC8pro_to_C6A_C8A_kc', 1.0)

Rule('bind_C8A_C3pro_to_C8AC3pro', C8(bf=None, state='A') + C3(bf=None, state='pro') <> C8(bf=1, state='A') % C3(bf=1, state='pro'), bind_C8A_C3pro_to_C8AC3pro_kf, bind_C8A_C3pro_to_C8AC3pro_kr)
Rule('catalyze_C8AC3pro_to_C8A_C3A', C8(bf=1, state='A') % C3(bf=1, state='pro') >> C8(bf=None, state='A') + C3(bf=None, state='A'), catalyze_C8AC3pro_to_C8A_C3A_kc)
Rule('bind_XIAP_C3A_to_XIAPC3A', XIAP(bf=None) + C3(bf=None, state='A') <> XIAP(bf=1) % C3(bf=1, state='A'), bind_XIAP_C3A_to_XIAPC3A_kf, bind_XIAP_C3A_to_XIAPC3A_kr)
Rule('catalyze_XIAPC3A_to_XIAP_C3ub', XIAP(bf=1) % C3(bf=1, state='A') >> XIAP(bf=None) + C3(bf=None, state='ub'), catalyze_XIAPC3A_to_XIAP_C3ub_kc)
Rule('bind_C3A_PARPU_to_C3APARPU', C3(bf=None, state='A') + PARP(bf=None, state='U') <> C3(bf=1, state='A') % PARP(bf=1, state='U'), bind_C3A_PARPU_to_C3APARPU_kf, bind_C3A_PARPU_to_C3APARPU_kr)
Rule('catalyze_C3APARPU_to_C3A_PARPC', C3(bf=1, state='A') % PARP(bf=1, state='U') >> C3(bf=None, state='A') + PARP(bf=None, state='C'), catalyze_C3APARPU_to_C3A_PARPC_kc)
Rule('bind_C3A_C6pro_to_C3AC6pro', C3(bf=None, state='A') + C6(bf=None, state='pro') <> C3(bf=1, state='A') % C6(bf=1, state='pro'), bind_C3A_C6pro_to_C3AC6pro_kf, bind_C3A_C6pro_to_C3AC6pro_kr)
Rule('catalyze_C3AC6pro_to_C3A_C6A', C3(bf=1, state='A') % C6(bf=1, state='pro') >> C3(bf=None, state='A') + C6(bf=None, state='A'), catalyze_C3AC6pro_to_C3A_C6A_kc)
Rule('bind_C6A_C8pro_to_C6AC8pro', C6(bf=None, state='A') + C8(bf=None, state='pro') <> C6(bf=1, state='A') % C8(bf=1, state='pro'), bind_C6A_C8pro_to_C6AC8pro_kf, bind_C6A_C8pro_to_C6AC8pro_kr)
Rule('catalyze_C6AC8pro_to_C6A_C8A', C6(bf=1, state='A') % C8(bf=1, state='pro') >> C6(bf=None, state='A') + C8(bf=None, state='A'), catalyze_C6AC8pro_to_C6A_C8A_kc)


Parameter('equilibrate_BidT_to_BidM_kf', 0.1)
Parameter('equilibrate_BidT_to_BidM_kr', 0.001)
Parameter('bind_BidM_BaxC_to_BidMBaxC_kf', 1e-07)
Parameter('bind_BidM_BaxC_to_BidMBaxC_kr', 0.001)
Parameter('catalyze_BidMBaxC_to_BidM_BaxM_kc', 1.0)
Parameter('bind_BidM_BaxM_to_BidMBaxM_kf', 1e-07)
Parameter('bind_BidM_BaxM_to_BidMBaxM_kr', 0.001)
Parameter('catalyze_BidMBaxM_to_BidM_BaxA_kc', 1.0)
Parameter('bind_BidM_BakM_to_BidMBakM_kf', 1e-07)
Parameter('bind_BidM_BakM_to_BidMBakM_kr', 0.001)
Parameter('catalyze_BidMBakM_to_BidM_BakA_kc', 1.0)
Parameter('bind_BaxA_BaxM_to_BaxABaxM_kf', 1e-07)
Parameter('bind_BaxA_BaxM_to_BaxABaxM_kr', 0.001)
Parameter('catalyze_BaxABaxM_to_BaxA_BaxA_kc', 1.0)
Parameter('bind_BakA_BakM_to_BakABakM_kf', 1e-07)
Parameter('bind_BakA_BakM_to_BakABakM_kr', 0.001)
Parameter('catalyze_BakABakM_to_BakA_BakA_kc', 1.0)
Rule('equilibrate_BidT_to_BidM', Bid(bf=None, state='T') <> Bid(bf=None, state='M'), equilibrate_BidT_to_BidM_kf, equilibrate_BidT_to_BidM_kr)
Rule('bind_BidM_BaxC_to_BidMBaxC', Bid(bf=None, state='M') + Bax(bf=None, state='C') <> Bid(bf=1, state='M') % Bax(bf=1, state='C'), bind_BidM_BaxC_to_BidMBaxC_kf, bind_BidM_BaxC_to_BidMBaxC_kr)
Rule('catalyze_BidMBaxC_to_BidM_BaxM', Bid(bf=1, state='M') % Bax(bf=1, state='C') >> Bid(bf=None, state='M') + Bax(bf=None, state='M'), catalyze_BidMBaxC_to_BidM_BaxM_kc)
Rule('bind_BidM_BaxM_to_BidMBaxM', Bid(bf=None, state='M') + Bax(bf=None, state='M') <> Bid(bf=1, state='M') % Bax(bf=1, state='M'), bind_BidM_BaxM_to_BidMBaxM_kf, bind_BidM_BaxM_to_BidMBaxM_kr)
Rule('catalyze_BidMBaxM_to_BidM_BaxA', Bid(bf=1, state='M') % Bax(bf=1, state='M') >> Bid(bf=None, state='M') + Bax(bf=None, state='A'), catalyze_BidMBaxM_to_BidM_BaxA_kc)
Rule('bind_BidM_BakM_to_BidMBakM', Bid(bf=None, state='M') + Bak(bf=None, state='M') <> Bid(bf=1, state='M') % Bak(bf=1, state='M'), bind_BidM_BakM_to_BidMBakM_kf, bind_BidM_BakM_to_BidMBakM_kr)
Rule('catalyze_BidMBakM_to_BidM_BakA', Bid(bf=1, state='M') % Bak(bf=1, state='M') >> Bid(bf=None, state='M') + Bak(bf=None, state='A'), catalyze_BidMBakM_to_BidM_BakA_kc)
Rule('bind_BaxA_BaxM_to_BaxABaxM', Bax(bf=None, s1=None, s2=None, state='A') + Bax(bf=None, state='M') <> Bax(bf=1, s1=None, s2=None, state='A') % Bax(bf=1, state='M'), bind_BaxA_BaxM_to_BaxABaxM_kf, bind_BaxA_BaxM_to_BaxABaxM_kr)
Rule('catalyze_BaxABaxM_to_BaxA_BaxA', Bax(bf=1, s1=None, s2=None, state='A') % Bax(bf=1, state='M') >> Bax(bf=None, s1=None, s2=None, state='A') + Bax(bf=None, state='A'), catalyze_BaxABaxM_to_BaxA_BaxA_kc)
Rule('bind_BakA_BakM_to_BakABakM', Bak(bf=None, s1=None, s2=None, state='A') + Bak(bf=None, state='M') <> Bak(bf=1, s1=None, s2=None, state='A') % Bak(bf=1, state='M'), bind_BakA_BakM_to_BakABakM_kf, bind_BakA_BakM_to_BakABakM_kr)
Rule('catalyze_BakABakM_to_BakA_BakA', Bak(bf=1, s1=None, s2=None, state='A') % Bak(bf=1, state='M') >> Bak(bf=None, s1=None, s2=None, state='A') + Bak(bf=None, state='A'), catalyze_BakABakM_to_BakA_BakA_kc)

Parameter('bind_BidM_Bcl2M_kf', 1e-06)
Parameter('bind_BidM_Bcl2M_kr', 0.06600000000000002)
Parameter('bind_BidM_BclxLM_kf', 1e-06)
Parameter('bind_BidM_BclxLM_kr', 0.012000000000000002)
Parameter('bind_BidM_Mcl1M_kf', 1e-06)
Parameter('bind_BidM_Mcl1M_kr', 0.010000000000000002)
Parameter('bind_BaxA_Bcl2_kf', 1e-06)
Parameter('bind_BaxA_Bcl2_kr', 0.010000000000000002)
Parameter('bind_BaxA_BclxLM_kf', 1e-06)
Parameter('bind_BaxA_BclxLM_kr', 0.010000000000000002)
Parameter('bind_BakA_BclxLM_kf', 1e-06)
Parameter('bind_BakA_BclxLM_kr', 0.05)
Parameter('bind_BakA_Mcl1M_kf', 1e-06)
Parameter('bind_BakA_Mcl1M_kr', 0.010000000000000002)
Parameter('bind_BadM_Bcl2_kf', 1e-06)
Parameter('bind_BadM_Bcl2_kr', 0.011000000000000001)
Parameter('bind_BadM_BclxLM_kf', 1e-06)
Parameter('bind_BadM_BclxLM_kr', 0.010000000000000002)
Parameter('bind_NoxaM_Mcl1M_kf', 1e-06)
Parameter('bind_NoxaM_Mcl1M_kr', 0.019000000000000006)

Rule('bind_BidM_Bcl2M', Bid(bf=None, state='M') + Bcl2(bf=None, state='M') <> Bid(bf=1, state='M') % Bcl2(bf=1, state='M'), bind_BidM_Bcl2M_kf, bind_BidM_Bcl2M_kr)
Rule('bind_BidM_BclxLM', Bid(bf=None, state='M') + BclxL(bf=None, state='M') <> Bid(bf=1, state='M') % BclxL(bf=1, state='M'), bind_BidM_BclxLM_kf, bind_BidM_BclxLM_kr)
Rule('bind_BidM_Mcl1M', Bid(bf=None, state='M') + Mcl1(bf=None, state='M') <> Bid(bf=1, state='M') % Mcl1(bf=1, state='M'), bind_BidM_Mcl1M_kf, bind_BidM_Mcl1M_kr)
Rule('bind_BaxA_Bcl2', Bax(bf=None, s1=None, s2=None, state='A') + Bcl2(bf=None) <> Bax(bf=1, s1=None, s2=None, state='A') % Bcl2(bf=1), bind_BaxA_Bcl2_kf, bind_BaxA_Bcl2_kr)
Rule('bind_BaxA_BclxLM', Bax(bf=None, s1=None, s2=None, state='A') + BclxL(bf=None, state='M') <> Bax(bf=1, s1=None, s2=None, state='A') % BclxL(bf=1, state='M'), bind_BaxA_BclxLM_kf, bind_BaxA_BclxLM_kr)
Rule('bind_BakA_BclxLM', Bak(bf=None, s1=None, s2=None, state='A') + BclxL(bf=None, state='M') <> Bak(bf=1, s1=None, s2=None, state='A') % BclxL(bf=1, state='M'), bind_BakA_BclxLM_kf, bind_BakA_BclxLM_kr)
Rule('bind_BakA_Mcl1M', Bak(bf=None, s1=None, s2=None, state='A') + Mcl1(bf=None, state='M') <> Bak(bf=1, s1=None, s2=None, state='A') % Mcl1(bf=1, state='M'), bind_BakA_Mcl1M_kf, bind_BakA_Mcl1M_kr)
Rule('bind_BadM_Bcl2', Bad(bf=None, state='M') + Bcl2(bf=None) <> Bad(bf=1, state='M') % Bcl2(bf=1), bind_BadM_Bcl2_kf, bind_BadM_Bcl2_kr)
Rule('bind_BadM_BclxLM', Bad(bf=None, state='M') + BclxL(bf=None, state='M') <> Bad(bf=1, state='M') % BclxL(bf=1, state='M'), bind_BadM_BclxLM_kf, bind_BadM_BclxLM_kr)
Rule('bind_NoxaM_Mcl1M', Noxa(bf=None, state='M') + Mcl1(bf=None, state='M') <> Noxa(bf=1, state='M') % Mcl1(bf=1, state='M'), bind_NoxaM_Mcl1M_kf, bind_NoxaM_Mcl1M_kr)

Parameter('assemble_pore_sequential_Bax_2_kf', 0.0002040816)
Parameter('assemble_pore_sequential_Bax_2_kr', 0.001)
Parameter('assemble_pore_sequential_Bax_3_kf', 0.0002040816)
Parameter('assemble_pore_sequential_Bax_3_kr', 0.001)
Parameter('assemble_pore_sequential_Bax_4_kf', 0.0002040816)
Parameter('assemble_pore_sequential_Bax_4_kr', 0.001)
Parameter('assemble_pore_sequential_Bak_2_kf', 0.0002040816)
Parameter('assemble_pore_sequential_Bak_2_kr', 0.001)
Parameter('assemble_pore_sequential_Bak_3_kf', 0.0002040816)
Parameter('assemble_pore_sequential_Bak_3_kr', 0.001)
Parameter('assemble_pore_sequential_Bak_4_kf', 0.0002040816)
Parameter('assemble_pore_sequential_Bak_4_kr', 0.001)
Parameter('pore_transport_complex_BaxA_4_CytoCM_kf', 2.857143e-05)
Parameter('pore_transport_complex_BaxA_4_CytoCM_kr', 0.001)
Parameter('pore_transport_dissociate_BaxA_4_CytoCC_kc', 10.0)
Parameter('pore_transport_complex_BaxA_4_SmacM_kf', 2.857143e-05)
Parameter('pore_transport_complex_BaxA_4_SmacM_kr', 0.001)
Parameter('pore_transport_dissociate_BaxA_4_SmacC_kc', 10.0)
Parameter('pore_transport_complex_BakA_4_CytoCM_kf', 2.857143e-05)
Parameter('pore_transport_complex_BakA_4_CytoCM_kr', 0.001)
Parameter('pore_transport_dissociate_BakA_4_CytoCC_kc', 10.0)
Parameter('pore_transport_complex_BakA_4_SmacM_kf', 2.857143e-05)
Parameter('pore_transport_complex_BakA_4_SmacM_kr', 0.001)
Parameter('pore_transport_dissociate_BakA_4_SmacC_kc', 10.0)

Rule('assemble_pore_sequential_Bax_2', Bax(bf=None, s1=None, s2=None, state='A') + Bax(bf=None, s1=None, s2=None, state='A') <> Bax(bf=None, s1=None, s2=1, state='A') % Bax(bf=None, s1=1, s2=None, state='A'), assemble_pore_sequential_Bax_2_kf, assemble_pore_sequential_Bax_2_kr)
Rule('assemble_pore_sequential_Bax_3', Bax(bf=None, s1=None, s2=None, state='A') + Bax(bf=None, s1=None, s2=1, state='A') % Bax(bf=None, s1=1, s2=None, state='A') <> MatchOnce(Bax(bf=None, s1=3, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A')), assemble_pore_sequential_Bax_3_kf, assemble_pore_sequential_Bax_3_kr)
Rule('assemble_pore_sequential_Bax_4', Bax(bf=None, s1=None, s2=None, state='A') + MatchOnce(Bax(bf=None, s1=3, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A')) <> MatchOnce(Bax(bf=None, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A')), assemble_pore_sequential_Bax_4_kf, assemble_pore_sequential_Bax_4_kr)
Rule('assemble_pore_sequential_Bak_2', Bak(bf=None, s1=None, s2=None, state='A') + Bak(bf=None, s1=None, s2=None, state='A') <> Bak(bf=None, s1=None, s2=1, state='A') % Bak(bf=None, s1=1, s2=None, state='A'), assemble_pore_sequential_Bak_2_kf, assemble_pore_sequential_Bak_2_kr)
Rule('assemble_pore_sequential_Bak_3', Bak(bf=None, s1=None, s2=None, state='A') + Bak(bf=None, s1=None, s2=1, state='A') % Bak(bf=None, s1=1, s2=None, state='A') <> MatchOnce(Bak(bf=None, s1=3, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A')), assemble_pore_sequential_Bak_3_kf, assemble_pore_sequential_Bak_3_kr)
Rule('assemble_pore_sequential_Bak_4', Bak(bf=None, s1=None, s2=None, state='A') + MatchOnce(Bak(bf=None, s1=3, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A')) <> MatchOnce(Bak(bf=None, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A')), assemble_pore_sequential_Bak_4_kf, assemble_pore_sequential_Bak_4_kr)
Rule('pore_transport_complex_BaxA_4_CytoCM', MatchOnce(Bax(bf=None, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A')) + CytoC(bf=None, state='M') <> MatchOnce(Bax(bf=5, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A') % CytoC(bf=5, state='M')), pore_transport_complex_BaxA_4_CytoCM_kf, pore_transport_complex_BaxA_4_CytoCM_kr)
Rule('pore_transport_dissociate_BaxA_4_CytoCC', MatchOnce(Bax(bf=5, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A') % CytoC(bf=5, state='M')) >> MatchOnce(Bax(bf=None, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A')) + CytoC(bf=None, state='C'), pore_transport_dissociate_BaxA_4_CytoCC_kc)
Rule('pore_transport_complex_BaxA_4_SmacM', MatchOnce(Bax(bf=None, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A')) + Smac(bf=None, state='M') <> MatchOnce(Bax(bf=5, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A') % Smac(bf=5, state='M')), pore_transport_complex_BaxA_4_SmacM_kf, pore_transport_complex_BaxA_4_SmacM_kr)
Rule('pore_transport_dissociate_BaxA_4_SmacC', MatchOnce(Bax(bf=5, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A') % Smac(bf=5, state='M')) >> MatchOnce(Bax(bf=None, s1=4, s2=1, state='A') % Bax(bf=None, s1=1, s2=2, state='A') % Bax(bf=None, s1=2, s2=3, state='A') % Bax(bf=None, s1=3, s2=4, state='A')) + Smac(bf=None, state='C'), pore_transport_dissociate_BaxA_4_SmacC_kc)
Rule('pore_transport_complex_BakA_4_CytoCM', MatchOnce(Bak(bf=None, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A')) + CytoC(bf=None, state='M') <> MatchOnce(Bak(bf=5, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A') % CytoC(bf=5, state='M')), pore_transport_complex_BakA_4_CytoCM_kf, pore_transport_complex_BakA_4_CytoCM_kr)
Rule('pore_transport_dissociate_BakA_4_CytoCC', MatchOnce(Bak(bf=5, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A') % CytoC(bf=5, state='M')) >> MatchOnce(Bak(bf=None, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A')) + CytoC(bf=None, state='C'), pore_transport_dissociate_BakA_4_CytoCC_kc)
Rule('pore_transport_complex_BakA_4_SmacM', MatchOnce(Bak(bf=None, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A')) + Smac(bf=None, state='M') <> MatchOnce(Bak(bf=5, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A') % Smac(bf=5, state='M')), pore_transport_complex_BakA_4_SmacM_kf, pore_transport_complex_BakA_4_SmacM_kr)
Rule('pore_transport_dissociate_BakA_4_SmacC', MatchOnce(Bak(bf=5, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A') % Smac(bf=5, state='M')) >> MatchOnce(Bak(bf=None, s1=4, s2=1, state='A') % Bak(bf=None, s1=1, s2=2, state='A') % Bak(bf=None, s1=2, s2=3, state='A') % Bak(bf=None, s1=3, s2=4, state='A')) + Smac(bf=None, state='C'), pore_transport_dissociate_BakA_4_SmacC_kc)


generate_equations(model)

Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('C8a_obs', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub')
           % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='pro'))
Observable('RIP1_RIP3_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')
           % RIP3(bRHIM=5, bDD = None, state='unmod'))
Observable('flip_necro',TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub')
           % FADD(bDD=1,bDED1 = 2, bDED2 = None) % flip(bf=2))
Observable('mBid', Bid(state='M'))
Observable('aSmac', Smac(state='A'))
Observable('cPARP', PARP(state='C'))

flip_num = list([3900, 9000, 20000])
c8_num = list([9000, 3900, 20000])


# tnf = [2326]
# color = ['r', 'm', 'g', 'b']
# lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
# initials = {flip_L(bDED=None, state = 'A'): flipnum, C8(bf=None, flip = None, state='I'): c8num}
tspan = np.linspace(0, 72000, 72000)
sim = BngSimulator(model, tspan=tspan)
result = sim.run(method='ode')
# df = result.dataframe

# plt.figure()
# for n in range(0,1):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], '--', c = color[n],lw = 1.5) #fppf
# # plt.scatter(x, mlkl, color = 'tab:gray', marker = 's')
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses'
# plt.legend(flip_n, title = 'flip', loc=0, fontsize = 5)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FP', loc=0, fontsize = 5)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# # plt.show()

plt.figure(figsize = (15,10))
# plt.title('High C8, low Flip')
plt.subplot(221)
plt.plot(tspan, result.observables['MLKLa_obs'],label = 'mlkla', marker = 's')
plt.xlabel('Time [hrs]', fontsize=16)
plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.legend(loc='best', fontsize = 10)

# plt.figure()
plt.subplot(222)
plt.plot(tspan, result.observables['C8a_obs'], label = 'C8a', marker = 's')
plt.xlabel('Time [hr]', fontsize=16)
plt.ylabel('C8a [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.legend(loc='best', fontsize = 10)

# plt.figure()
plt.subplot(223)
plt.plot(tspan, result.observables['RIP1_RIP3_obs'],label = 'RIP1:RIP3', marker = 's')
plt.xlabel('Time [hr]', fontsize=16)
plt.ylabel('RIP1:RIP3 [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.legend(loc='best', fontsize = 10)

# plt.figure()
plt.subplot(224)
plt.plot(tspan, result.observables['flip_necro'],label = 'flip_necro', marker = 's')
plt.xlabel('Time [hr]', fontsize=16)
plt.ylabel('flip_necro [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.legend(loc='best', fontsize = 10)

plt.figure(figsize = (15,8))
plt.subplot(131)
plt.plot(tspan, result.observables['mBid'][:])
# plt.scatter(x, mlkl, color = 'r')
plt.xlabel('Time [seconds]', fontsize=16)
plt.ylabel('mBid amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.ylim(ymax = 10000)

plt.subplot(132)
plt.plot(tspan, result.observables['aSmac'][:])
# plt.scatter(x, mlkl, color = 'r')
plt.xlabel('Time [seconds]', fontsize=16)
plt.ylabel('aSmac amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.ylim(ymax = 10000)

plt.subplot(133)
plt.plot(tspan, result.observables['cPARP'][:])
# plt.scatter(x, mlkl, color = 'r')
plt.xlabel('Time [seconds]', fontsize=16)
plt.ylabel('cParp amount [molecules]', fontsize=16)

plt.show()

# quit()
#
# plt.figure()
# plt.plot(tspan/60, result.observables['MLKLa_obs'], color = 'tab:gray', marker = 's')
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.legend('100 ng/ml ' , title = 'TNF FP', loc=0, fontsize = 5)
# plt.show()
#
#
# quit()
# mlkl = [0, 170, 900, 4880, 9940, 10000]
# x = [0, 60, 120, 240, 360, 480]
#
# plt.figure()
# for n in range(0,1):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], '--', c = color[n],lw = 1.5) #fppf
# plt.scatter(x, mlkl, color = 'tab:gray', marker = 's')
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FP', loc=0, fontsize = 5)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# plt.show()
#
#
#
# # sim = BngSimulator(model, tspan=tspan)
# # result = sim.run(method='ode', param_values=fstpso)
# #
# # mlkl = [0, 170, 900, 4880, 9940, 10000]
# # x = [0, 60, 120, 240, 360, 480]
# #
# # plt.figure()
# # plt.plot(tspan, result.observables['MLKLa_obs'][:])
# # plt.scatter(x, mlkl, color = 'r')
# # plt.xlabel('Time [minutes]', fontsize=16)
# # plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# # # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.ylim(ymax = 10000)
# # plt.show()
