from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.backends.backend_pdf import PdfPages
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
Monomer('flip_L', ['bDED', 'state'], {'state': ['A', 'I']})
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('TAK1', ['brip', 'bmapk'])
Monomer('NEMO', ['brip', 'btak', 'bikk', 'state'], {'state': ['I', 'A']})
Monomer('LUBAC', ['brip'])

Parameter('TNF_0', 2326) #698 is 30ng/ml of TNF
Parameter('TNFR_0', 4800) #0.00246
Parameter('TRADD_0', 9000) #
Parameter('RIP1_0', 40000) #47000 0.04
Parameter('TRAF_0', 9000) # 8.3e-4
Parameter('cIAP_0', 9000) #10000 8.3e-4
Parameter('A20_0', 9000) #2256
Parameter('CYLD_0', 9000) #50000 # 0.004
Parameter('FADD_0', 8030) # 0.0033
Parameter('flip_L_0', 3090) # 0.004 # 0.09
Parameter('Lubac_0', 7226)
Parameter('C8_0', 9000) #10000 # 0.033 # 0.093 # 0.0107
Parameter('RIP3_0', 40000) #20000
Parameter('NEMO_0', 24000) # 1000000
Parameter('MLKLa_0', 20000) # 100000 #0.0034

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
Initial(NEMO(brip = None, btak = None, bikk = None, state = 'I'), NEMO_0)
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
     >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     k_CYLD_1)

#Initiating Necroptosis
# Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 1e-1)
# Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 3.11e-7)
Parameter('bind_FADD_proC8_2_kf', 7.27e-06)
Parameter('bind_FADD_proC8_2_kr', 0.018)
Parameter('kc_c8_1', 1e-1)
Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 1e-1)
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 3.11e-7)

Parameter('bind_FADD_proC8_kf', 7.27e-06)
Parameter('bind_FADD_proC8_kr', 0.018)

# Parameter('bind_FADDANY_flip_L_kf',7.27e-06)
# Parameter('bind_FADDANY_flip_L_kr', 0.018)
# Parameter('bind_C8_flip_L_kf',7.27e-06)
# Parameter('bind_C8_flip_L_kr', 0.018)

# Parameter('bind_FADDANY_RIP3_kf', 1e-6)
# Parameter('bind_FADDANY_RIP3_kr', 1e-3)
Parameter('kc_c8_2', 1e-1)

#RIP1 deub and necrosome formation

#FADD + C8i <> FADD : C8i
#FADD:C8i + Flip <> FADD:C8i:Flip
#FADD:C8i:Flip >> FADD:C8a%Flip
#FADD:C8a%Flip + RIP1 <> FADD:C8a%Flip%RIP1db
#FADD:C8a%Flip%RIP1db >> FADD + C8i + Flip + RIP1trunc
#
#RIP1db + RIP3 <> RIP1db:RIP3



# Rule('bind_TRADDANYRIP1ANY_FADD', C8(bf=None, flip = None, state='I') <> FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=1,flip = None, state='I'),bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)

Rule('bind_C8_Flip', C8(bf=None,flip = None, state='I') + flip_L(bDED=None, state = 'A') <>  C8(bf=None,flip = 2, state='I') % flip_L(bDED=2, state = 'A'),bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)

Rule('C8_cat_Flip',  C8(bf=None,flip = 2, state='I') % flip_L(bDED=2, state = 'A') >> C8(bf=None,flip = 2, state='A') % flip_L(bDED=2, state = 'A'), kc_c8_1)



# Rule('bind_Rip1', FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=1,flip = 2, state='A') % flip_L(bDED=2, state = 'A') + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') <>
#      FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=1,flip = 2, state='A') % flip_L(bDED=2, state = 'A') % RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub'), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)





     # <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)

Rule('bind_RIP1_FADD', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None,bDED1 = None, bDED2 = None) <>
     RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None),bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)


Rule('bind_FADD_proC8_2',  RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None,flip = 2, state='A') % flip_L(bDED=2, state = 'A')
     <> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=3,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=3,flip = 2, state='A') % flip_L(bDED=2, state = 'A'), bind_FADD_proC8_kf, bind_FADD_proC8_kr)

Rule('Rip1_trunc_C8Flip', RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=3,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=3,flip = 2, state='A') % flip_L(bDED=2, state = 'A') >>
     RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,flip = None, state='I') + flip_L(bDED=None, state = 'A'), kc_c8_2)




# Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2, flip = None,state='I') + flip_L(bDED=None, state = 'A')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A'), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
#
# Rule('activate_C8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A')
#      >> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A'), bind_C8_flip_L_kf, bind_C8_flip_L_kr)
#
# #@TODO ADDED THESE REACTION
# Rule('catalyze_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A') >>
#      TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,flip = None, state='I') + flip_L(bDED=None, state = 'A') , kc_c8_1)
#
#
#
#
# Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_RIP3_kf, bind_FADDANY_RIP3_kr)



# Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)

# Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)

#@TODO ORIGINAL RULE
# Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)


Parameter('bindrip1rip3kf', 1e-3)
Parameter('bindrip1rip3kr', 1e-6)
Parameter('bind_RIP1_RIP3po4_kf', 1e-3)
Parameter('RIP1po4_RIP3po4_kf', 1e-6)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1e-3)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 1e-6)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 1)

Rule('Rip1_bind_Rip3', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = None, state = 'deub') + RIP3(bRHIM=None, bDD = None, state='unmod') <>
     RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), bindrip1rip3kf, bindrip1rip3kr)

Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_RIP1_RIP3po4_kf)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), RIP1po4_RIP3po4_kf)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
     <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)



Observable('RIP1_obs',RIP1(state='unmod'))
Observable('RIP1k63_obs',RIP1(state='K63ub'))
Observable('RIP1deub_obs',RIP1(state='deub'))
Observable('RIP1po4_obs',RIP1(state='po4'))
Observable('RIP3_obs', RIP3(state='unmod'))
Observable('MLKL_obs', MLKL(state='unmod'))
Observable('MLKLa_obs', MLKL(state='active'))
Observable('RIP13_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('RIP1deub3po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
           RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))
Observable('A20_obs', A20(brip = None))
Observable('C8i_obs', C8(state = 'I'))
Observable('C8a_obs', C8(state = 'A'))
Observable('flip_obs', flip_L())
Observable('CIIa_obs', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=2,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A'))



# params = np.array([-5.30443515, -3.97996083, -4.13834667, -5.56469212, -3.12171887, -5.88636881
#          -2.09441601, -4.02225288, -3.37499851, -5.51352587, -3.17279288, -1.66591288
#          -5.47233099, -3.53222612, -7.9528404 , -1.44918836, -0.60323906, -5.56448962,
#          -4.93576555, -1.56250485, -0.04862364, -7.60935808, -5.99280593,  0.39516239,
#          -6.15132573, -0.12083605,  0.71538332, -6.3165361 , -4.62594487, -1.12642492,
#          -1.63085357, -1.66712655, -3., -2.34609272,  0.17064889])
# new_params = 10**params
# print(new_params)
# initials = [600, 4800, 9000, 20000, 9000, 9000,9000, 9000, 8030, 3900, 7226, 160000, 40000, 24000, 4800]
# params = [-3.26583350782,	-5.90793122272,	-6.91071484206, 7.26683657118,	1.19527275645,	1.15773328597,	-1.62148262759,	-4.14988114453,
#           1.37940222603,	3.81916083957,	9.63025632674,	0.793327429331,	-6.02916249511,	-6.84484359606,	-7.19312555916,	7.84528125124,
#           -2.16827411435,	-3.64756788082,	7.59655331234,	-6.78811071013,	-3.76657563219,	9.3434940999,	-7.36652931984,	6.27924286688,
#           -8.61356215219,	-7.68937560312,	-5.70859996861,	7.45745204635,	-8.29400463637,	9.83045337456,	7.51255180804,	2.92830496104,
#           6.49783223747,	-7.99931100374,	-8.97680701083]

# tspan = np.linspace(0, 720, 721)
#
# sim = ScipyOdeSimulator(model, tspan)
# # sim_result = sim.run()
#
#
#
# # sim_result = sim.run(initials = {RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'): rip1num, C8(bf=None, flip = None, state='A'): c8num, flip_L(bDED=None, state = 'A'): flipnum}, param_values=params)
# sim_result = sim.run(initials = {FADD(bDD=None, bDED1=None, bDED2=None): faddnum}, param_values=params)
# sim_df = sim_result.dataframe
#
# for n in range(0,20):
#     plt.plot(tspan/60, sim_df.loc[n]['MLKLa_obs'].iloc[:], lw =1)
# plt.xlabel('Time in Hr', fontsize=16)
# plt.ylabel('MLKL uM Concentration', fontsize=16)
# plt.title('TNF:0.1 ng/ml')
# plt.legend(faddnum, title = 'Fadd', loc='best', fontsize = 5)
# plt.show()
# quit()

# print(sim_df)
# print(model.parameters)
# print(type(model.parameters_rules().values))
# print(sim_df)
# quit()
#
# print(sim_df.loc[:,:]['C8a_obs]'].iloc[:])
# quit()
#
# print(model.parameters)
# # quit()
# print(len(model.parameters))
# print(len(model.rules))

# print('tnf')
# print(sim_df.shape)
# quit()
# print(sim_df.loc[:]['MLKLa_obs'].iloc[:])
# plt.figure()
# plt.plot(tspan/60, sim_df.loc[:]['MLKLa_obs'].iloc[:], lw =0.5)
# plt.show()
# quit()
# print('mlkl2')
# print(sim_df[1]['MLKLa_obs'].iloc[:])
# quit()


# sim1 = ScipyOdeSimulator(model, tspan)
# sim2 =ScipyOdeSimulator(model, tspan)
# sim3 = ScipyOdeSimulator(model, tspan)
# sim4 =ScipyOdeSimulator(model, tspan=tspan)

# L1 = sim1.run(initials= {TNF(brec=None): 1.96e-5},param_values=params)
# L2 = sim2.run(initials= {TNF(brec=None): 1.96e-4},param_values=params)
# L3 = sim3.run(initials= initials,param_values=params)
# L4 = sim4.run()

# L1 = sim1.run(initials= {TNF(brec=None): 47})
# L2 = sim2.run(initials= {TNF(brec=None): 472})
# L3 = sim3.run(initials= {TNF(brec=None): 4721})
# L4 = sim4.run(initials= {TNF(brec=None): 47212})

# print('mlkl')
# print(L1.observables['MLKLa_obs'])
# print(L2.observables['MLKLa_obs'])
# print(L3.observables['MLKLa_obs'])
# print(L4.observables['MLKLa_obs'])
#
# print('c8i')
# print(L1.observables['C8i_obs'])
# print(L2.observables['C8i_obs'])
# print(L3.observables['C8i_obs'])
# print(L4.observables['C8i_obs'])
#
# print('c8a')
# print(L1.observables['C8a_obs'])
# print(L2.observables['C8a_obs'])
# print(L3.observables['C8a_obs'])
# print(L4.observables['C8a_obs'])

#
# # plt.figure()
# # plt.plot(tspan/60, L1.observables['C8_obs'],label = 'C8.1')
# # plt.plot(tspan/60, L2.observables['C8_obs'],label = 'C81')
# # plt.plot(tspan/60, L3.observables['C8_obs'],label = 'C810')
# # plt.plot(tspan/60, L4.observables['C8_obs'],label = 'C8100')
# # plt.xlabel("Time (in hr)", fontsize=15)
# # plt.ylabel("Concentration uM", fontsize=15)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
# #
# # plt.figure()
# # plt.plot(tspan/60, L1.observables['flip_obs'],label = 'flip.1')
# # plt.plot(tspan/60, L2.observables['flip_obs'],label = 'flip1')
# # plt.plot(tspan/60, L3.observables['flip_obs'],label = 'flip10')
# # plt.plot(tspan/60, L4.observables['flip_obs'],label = 'flip100')
# # plt.xlabel("Time (in hr)", fontsize=15)
# # plt.ylabel("Concentration uM", fontsize=15)
# # # plt.ylim(ymin = -10, ymax =100)
# # plt.legend(loc=0)
#
# plt.figure(figsize = (18,7))
# # # plt.figure()
# plt.subplot(221)
# # # plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
# # # plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
# # plt.plot(tspan/60, L3.observables['TNF_obs'],label = 'TNF10d')
# plt.plot(tspan/60, L4.observables['TNF_obs'],label = 'TNF100s')
# #
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# #
# # # plt.figure()
# plt.subplot(222)
# # # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # # plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
# # # plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
# # plt.plot(tspan/60, L3.observables['CI_k63_obs'],label = 'CI_k6310d')
# plt.plot(tspan/60, L4.observables['CI_k63_obs'],label = 'CI_k63100s')
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# # plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
# # #
# # # plt.figure()
# plt.subplot(223)
# # plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# # plt.plot(tspan/60, L2.observables['RIP13po4_obs'],label = 'RIP13po41')
# # plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po410d')
# plt.plot(tspan/60, L4.observables['RIP13po4_obs'],label = 'RIP13po4100s')
#
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(224)
# # plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# # plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# # plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLa10d')
# plt.plot(tspan/60, L4.observables['MLKL_obs'],label = 'MLKL100s')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa100s')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# plt.legend(loc = 0)

# plt.figure(figsize = (30, 10))
# # plt.figure()
# plt.subplot(141)
# plt.plot(tspan/60, L1.observables['RIP1_obs'],label = 'RIP1.1')
# plt.plot(tspan/60, L2.observables['RIP1_obs'],label = 'RIP11')
# plt.plot(tspan/60, L3.observables['RIP1_obs'],label = 'RIP110')
# plt.plot(tspan/60, L4.observables['RIP1_obs'],label = 'RIP1100')
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymax =9e-4)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(142)
# plt.plot(tspan/60, L1.observables['RIP1k63_obs'],label = 'RIP1k63.1')
# plt.plot(tspan/60, L2.observables['RIP1k63_obs'],label = 'RIP1k631')
# plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k6310')
# plt.plot(tspan/60, L4.observables['RIP1k63_obs'],label = 'RIP1k63100')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(143)
# plt.plot(tspan/60, L1.observables['RIP1deub_obs'],label = 'RIP1deub.1')
# plt.plot(tspan/60, L2.observables['RIP1deub_obs'],label = 'RIP1deub1')
# plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub10')
# plt.plot(tspan/60, L4.observables['RIP1deub_obs'],label = 'RIP1deub100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(144)
# plt.plot(tspan/60, L1.observables['RIP1po4_obs'],label = 'RIP1po4.1')
# plt.plot(tspan/60, L2.observables['RIP1po4_obs'],label = 'RIP1po41')
# plt.plot(tspan/60, L3.observables['RIP1po4_obs'],label = 'RIP1po410')
# plt.plot(tspan/60, L4.observables['RIP1po4_obs'],label = 'RIP1po4100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.figure(figsize = (23,5))
# # plt.figure()
# plt.subplot(141)
# plt.plot(tspan/60, L1.observables['RIP13_obs'],label = 'RIP13.1')
# plt.plot(tspan/60, L2.observables['RIP13_obs'],label = 'RIP131')
# plt.plot(tspan/60, L3.observables['RIP13_obs'],label = 'RIP1310')
# plt.plot(tspan/60, L4.observables['RIP13_obs'],label = 'RIP13100')
# #
# # plt.plot(tspan/60, L1.observables['RIP1k63_obs'],label = 'RIP1k63.1')
# # plt.plot(tspan/60, L2.observables['RIP1k63_obs'],label = 'RIP1k631')
# # plt.plot(tspan/60, L3.observables['RIP1k63_obs'],label = 'RIP1k6310')
# # plt.plot(tspan/60, L4.observables['RIP1k63_obs'],label = 'RIP1k63100')
# # plt.plot(tspan/60, L1.observables['RIP1deub_obs'],label = 'RIP1deub.1')
# # plt.plot(tspan/60, L2.observables['RIP1deub_obs'],label = 'RIP1deub1')
# # plt.plot(tspan/60, L3.observables['RIP1deub_obs'],label = 'RIP1deub10')
# # plt.plot(tspan/60, L4.observables['RIP1deub_obs'],label = 'RIP1deub100')
# # plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymax =9e-4)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(142)
# plt.plot(tspan/60, L1.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po4.1')
# plt.plot(tspan/60, L2.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po41')
# plt.plot(tspan/60, L3.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po410')
# plt.plot(tspan/60, L4.observables['RIP1deub3po4_obs'],label = 'RIP1deub3po4100')
# # plt.plot(tspan/60, L1.observables['RIP1_obs'],label = 'RIP1.1')
# # plt.plot(tspan/60, L2.observables['RIP1_obs'],label = 'RIP11')
# # plt.plot(tspan/60, L3.observables['RIP1_obs'],label = 'RIP110')
# # plt.plot(tspan/60, L4.observables['RIP1_obs'],label = 'RIP1100')
# # plt.plot(tspan/60, L1.observables['RIP3_obs'],label = 'RIP3.1')
# # plt.plot(tspan/60, L2.observables['RIP3_obs'],label = 'RIP31')
# # plt.plot(tspan/60, L3.observables['RIP3_obs'],label = 'RIP310')
# # plt.plot(tspan/60, L4.observables['RIP3_obs'],label = 'RIP3100')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# # plt.figure()
# plt.subplot(143)
# plt.plot(tspan/60, L1.observables['RIP13po4_obs'],label = 'RIP13po4.1')
# plt.plot(tspan/60, L2.observables['RIP13po4_obs'],label = 'RIP13po41')
# plt.plot(tspan/60, L3.observables['RIP13po4_obs'],label = 'RIP13po410')
# plt.plot(tspan/60, L4.observables['RIP13po4_obs'],label = 'RIP13po4100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)
#
# plt.subplot(144)
# plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
# plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
# plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLa10')
# plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa100')
# # plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
# plt.xlabel("Time (in hr)", fontsize=15)
# plt.ylabel("Concentration uM", fontsize=15)
# # plt.ylim(ymin = -10, ymax =100)
# plt.legend(loc=0)

# plt.tight_layout()
# plt.show()