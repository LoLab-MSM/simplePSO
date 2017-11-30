from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
# from QQSB.numtools import simulatortools as st
from pysb.util import alias_model_components

Model()

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
Monomer('flip_L', ['bDED'])
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('TAK1', ['brip', 'bmapk'])
Monomer('NEMO', ['brip', 'btak', 'bikk', 'state'], {'state': ['I', 'A']})
Monomer('LUBAC', ['brip'])

Parameter('TNF_0',698 ) #698 is 30ng/ml of TNF
Parameter('TNFR_0', 4800) #0.00246
Parameter('TRADD_0', 9000) #
Parameter('RIP1_0', 40000) #47000 0.04
Parameter('TRAF_0', 9000) # 8.3e-4
Parameter('cIAP_0', 9000) #10000 8.3e-4
Parameter('A20_0', 2000) #2256
Parameter('CYLD_0', 9000) #50000 # 0.004
Parameter('FADD_0', 8030) # 0.0033
Parameter('flip_L_0', 3900) # 0.004 # 0.09
Parameter('Lubac_0', 7226)
Parameter('C8_0', 9000) #10000 # 0.033 # 0.093 # 0.0107
Parameter('RIP3_0', 40000) #20000
Parameter('NEMO_0', 24000) # 1000000
Parameter('MLKLa_0', 4800) # 100000 #0.0034

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
Initial(flip_L(bDED=None), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
Initial(NEMO(brip = None, btak = None, bikk = None, state = 'I'), NEMO_0)
Initial(C8(bf=None, flip = None, state='I'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)


#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
Parameter('bind_TNF_TNFR_kf', 1e-6)
Parameter('bind_TNF_TNFR_kr',1e-3)
Parameter('TNF_deg1', 0.01)
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

Rule('bind_TNFRANY_TRADD', TNFR(blig=None, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <>
     TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)

Rule('bind_TNFRANY_RIP1unmod', TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
     <> TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)

Rule('Complex_I_ubiquitylation1',TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
     <> TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)


Rule('Complex_I_ubiquitylation', TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
     <> TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)


# Rule('TNF_deg', TNF(brec = None) >> None, TNF_deg1)
#
# Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <>
#      TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)
#
# Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)
#
# Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)
#
#
# Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)


Parameter('CompI_UB2', 1e-1)
Parameter('bind_LUBAC_kf', 1e-6)
Parameter('bind_LUBAC_kr', 1e-3)
Parameter('bind_RIP1K63ubANY_A20_kf', 1e-6)
Parameter('bind_RIP1K63ubANY_A20_kr',1e-3)
Parameter('k_A20_1', 1e-1)
Parameter('bind_RIP1K63ubANY_CYLD_kf', 1e-6)
Parameter('bind_RIP1K63ubANY_CYLD_kr', 1e-3)
Parameter('k_CYLD_1', 1e-1)

Rule('Complex_I_ubiquitylation2',TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)


Rule('ComplexI_Lubac', TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     <> TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     bind_LUBAC_kf, bind_LUBAC_kr)


#RIP1 K63ub to be deub by A20 or CYLD

Rule('bind_RIP1K63ubANY_A20',TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
     <> TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
     bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)


Rule('A20_2',  TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
     >>  TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None), k_A20_1)



Rule('bind_RIP1K63ubANY_CYLD', TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
     <>  TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
     bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)

Rule('A20_1', TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
     >> TNFR(blig=None, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     k_CYLD_1)

#Initiating Necroptosis
Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 1e-1)
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 3.11e-7)

Parameter('bind_C8flip_kf', 7.27e-06)
Parameter('bind_C8flip_kr', 0.018)
Parameter('bind_C8flip_kc', 1e-1)

Parameter('bind_FADD_proC8_2_kf', 7.27e-06)
Parameter('bind_FADD_proC8_2_kr', 0.018)
Parameter('bind_FADDANY_flip_L_kf',7.27e-06)
Parameter('bind_FADDANY_flip_L_kr', 0.018)
Parameter('bind_C8_flip_L_kf',7.27e-06)
Parameter('bind_C8_flip_L_kr', 0.018)
Parameter('kc_c8_1', 1e-1)

Parameter('rip3kf', 1e-6)
Parameter('rip3kr', 1e-3)
Parameter('rip3kc', 1e-1)

Parameter('bind_FADDANY_RIP3_kf', 1e-6)
Parameter('bind_FADDANY_RIP3_kr', 1e-3)
Parameter('kc_c8_2', 1e-1)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     <>  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)


# Rule('C8_bind_Flip', flip_L(bDED=None) + C8(bf=None, flip = None, state='I') <> C8(bf=None,flip = 4, state='I') % flip_L(bDED=4), bind_C8flip_kf, bind_C8flip_kr)
#
# Rule('C8_cat_Flip', C8(bf=None,flip = 4, state='I') % flip_L(bDED=4) >> flip_L(bDED=None) + C8(bf=None, flip = None, state='A'), bind_C8flip_kc)

# Rule('bind_FADD_proC8_2',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, flip = None, state='A')
#      <> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = None, state='A'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


# Rule('bind_FADDANY_flip_L', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2, flip = None,state='A') + flip_L(bDED=None)
#      <> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
#
# Rule('activate_C8', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4)
#      >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4), bind_C8_flip_L_kf, bind_C8_flip_L_kr)

#@TODO ADDED THESE REACTION
# Rule('catalyze_FADDANY_flip_L',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = None, state='A')  >>
#       RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD=None, bRHIM=None, bMLKL=None, state='trunc') + FADD(bDD=None,bDED1=None,bDED2=None) + C8(bf=None,flip = None, state='I'), kc_c8_1)

#
Rule('bind_FADD_proC8_2', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, flip = None, state='I')
     <>  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = None, state='I'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


Rule('bind_FADDANY_flip_L',RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2, flip = None,state='I') + flip_L(bDED=None)
     <>  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)

Rule('activate_C8',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4)
     >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4), bind_C8_flip_L_kf, bind_C8_flip_L_kr)


#@TODO ADDED THESE REACTION
Rule('catalyze_FADDANY_flip_L',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4) >>
     RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD=None, bRHIM=None, bMLKL=None, state='trunc') + FADD(bDD=None,bDED1=None,bDED2=None) + C8(bf=None,flip = None, state='I') + flip_L(bDED=None) , kc_c8_1)
#


Rule('bind_RIP3', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4) + RIP3(bRHIM=None, bDD = None, state='unmod')
     <> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=5,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4) % RIP3(bRHIM=None, bDD =5, state='unmod'), rip3kf, rip3kr)

Rule('deg_C8Flip', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=5,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4) % RIP3(bRHIM=None, bDD =5, state='unmod')
     >>  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=5,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=None, bDD = 5, state='unmod'),rip3kc)


# Rule('bind_FADDANY_proC8', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
#      <> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_RIP3_kf, bind_FADDANY_RIP3_kr)



# Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)

# Rule('C8_activation2',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod') + C8(bf=None, flip = None, state='I')
#      >> FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)
#
#
Rule('C8_activation2',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=5,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=None, bDD = 5, state='unmod')
     >> FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=None, bDD = 5, state='unmod'), kc_c8_2)

#@TODO ORIGINAL RULE
# Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)



Parameter('bind_RIP1_RIP3po4_kf', 1e-2)
Parameter('RIP1po4_RIP3po4_kf', 1e-2)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1e-3)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 1e-6)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 1)

Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=None, bDD =5, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=None, bDD = 5, state='po4'), bind_RIP1_RIP3po4_kf)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'deub')% RIP3(bRHIM=None, bDD = 5, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'po4')% RIP3(bRHIM=None, bDD = 5, state='po4'), RIP1po4_RIP3po4_kf)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'po4')% RIP3(bRHIM=None, bDD = 5, state='po4') + MLKL(bRHIM=None, state='unmod')
     <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = 1, bRHIM = None, state = 'po4')% RIP3(bRHIM=None, bDD = 5, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 =5, bDD=None, btraf=None, bMLKL = 1, bRHIM = None, state = 'po4')% RIP3(bRHIM=None, bDD = 5, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'po4')% RIP3(bRHIM=None, bDD =5, state='po4') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)



Observable('RIP1_obs',RIP1(state='unmod'))
Observable('RIP1k63_obs',RIP1(state='K63ub'))
Observable('RIP1deub_obs',RIP1(state='deub'))
Observable('RIP1po4_obs',RIP1(state='po4'))
Observable('RIP3_obs', RIP3(state='unmod'))
Observable('MLKL_obs', MLKL(bRHIM=None, state='unmod'))
Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('RIP13_obs', RIP1(bDD=ANY, bRHIM=1, state='deub') % RIP3(bRHIM=1, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = 5, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'po4')% RIP3(bRHIM=None, bDD = 5, state='po4'))
Observable('RIP1deub3po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TNFR(blig=None, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
# Observable('CI_k63_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
#            RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))
Observable('A20_obs', A20(brip = None))
Observable('C8i_obs', C8(state = 'I'))
Observable('C8a_obs', C8(state = 'A'))
Observable('flip_obs', flip_L())
Observable('CIIa_obs', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4))

# generate_equations(model)


# params = np.array([-5.64547571, -3.29767049, -2.76677842, -6.0334102,  -3.38702609, -5.91678291,
#                  -3.98116991, -5.62540044, -2.22246753, -5.92784673, -3.23139993, -2.01283598,
#                  -5.58900452, -2.36915296, -5.5717866 , -2.50514325, -2.03819231, -6.26215161,
#                  -4.09439759, -1.67162902,  0.02095059, -6.12643541, -3.9660837 , -2.17589245,
#                  -4.73767298, -1.64582103, -4.30980095, -1.88843078,  0.20943653, -5.86856408,
#                  -2.10933906, -1.74440787, -1.8391347 , -1.97695361, -1.44301426, -5.13047439,
#                   1.69575111])

# params = np.array([
#             -5.35958193, -2.44645455, -3.67496492, -4.03363007, -0.12770523, -5.75531172,
#            -3.82257372, -3.88933494, -0.29141231, -5.50417801, -2.76435863, -1.71358552,
#            -5.64657099, -1.29242888, -3.94131499, -1.73861215,  2.        , -6.3967278 ,
#            -2.51372906,  0.9961207 , -2.30724548, -4.17312932, -5.22992234, -0.55118456,
#            -4.65935054, -3.21356214, -1.33307069, -3.31128076, -1.02698157, -1.88445244,
#            -1.87678634, -1.89851486, -3.        , -3.58667063, -0.09375467])
# pars = np.array([2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 4800,
#                    4.36936243e-06,   3.57721835e-03,   2.11365976e-0,   9.25486163e-05,
#                    7.45237620e-01 ,  1.75666230e-06 ,  1.50461809e-04,  1.29022383e-04 ,
#                    5.11196286e-01 ,  3.13200171e-06 ,  1.72044728e-03,  1.93381302e-02 ,
#                    2.25646712e-06 ,  5.10001108e-02 ,  1.14468241e-04,  1.82552527e-02 ,
#                    1.00000000e+02 ,  4.01118045e-07 ,  3.06387427e-03,  9.91107357e+00 ,
#                    4.92895122e-03 ,  6.71228951e-05 ,  5.88948961e-06,  2.81070612e-01 ,
#                    2.19103573e-05 ,  6.11558293e-04 ,  4.64439672e-02,  4.88336560e-04 ,
#                    9.39763190e-02 ,  1.30481085e-02 ,  1.32804766e-02,  1.26323788e-02 ,
#                    1.00000000e-03 ,  2.59017657e-04 ,  8.05833523e-01])
#
#
# print(len(model.parameters))
# print(len(model.parameters_rules()))
# print(len(model.initial_conditions))
# quit()
# initials = np.array([2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 4800])
# # pars = np.array([2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 2400, 4800,
# #           2.26216506e-06,  5.03882772e-04,   1.71088800e-03,   9.25954827e-07,
# #            4.10179461e-04,  1.21120342e-06,   1.04431157e-04,   2.36918819e-06,
# #            5.99145732e-03,  1.18073726e-06,   5.86948599e-04,   9.70876569e-03,
# #            2.57629434e-06,  4.27412324e-03,   2.68048512e-06,   3.12504841e-03,
# #            9.15814868e-03,  5.46825036e-07,   8.04641467e-05,   2.12995771e-02,
# #            1.04942303e+00,  7.47419785e-07,   1.08122555e-04,   6.66971920e-03,
# #            1.82947728e-05,  2.26036706e-02,   4.90003351e-05,   1.29291276e-02,
# #            1.61970726e+00,  1.35343038e-06,   7.77429364e-03,   1.80132522e-02,
# #            1.44832257e-02,  1.05449953e-02,   3.60566804e-02,   7.40500933e-06,
# #            4.96307811e+01, 1.0])
#
# model.enable_synth_deg()
# # print(len((model.parameters)))
# # print(pars.shape)
# # quit()
#
# tspan = np.linspace(0, 720, 721)
# sim = ScipyOdeSimulator(model, tspan=tspan)
# L4 = sim.run(param_values=pars)
# L5 = sim.run()
#
# plt.figure()
# plt.plot(tspan/60, L5.observables['MLKLa_obs'],color = 'blue',label = 'MLKLp')
# # plt.plot(tspan/60, L4.observables['MLKLa_obs'],color = 'orange',label = 'MLKLp_cal')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Molecules/Cell", fontsize=15)
# plt.legend(loc = 0)
# plt.show()
