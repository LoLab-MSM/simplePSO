from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

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
Monomer('C8', ['bf', 'state'], {'state': ['A', 'I']})
Monomer('flip_L', ['bDED'])
Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
Monomer('TAK1', ['brip', 'bmapk'])
Monomer('NEMO', ['brip', 'btak', 'bikk', 'state'], {'state': ['I', 'A']})
Monomer('LUBAC', ['brip'])
Monomer('IKK', ['bind', 'bnemo', 'state'], {'state': ['I', 'A', 'AI']})
Monomer('IkBa', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBa_mRNA')

# def ikbb_and_mRNA_monomers():
Monomer('IkBb', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBb_mRNA')

# def ikbe_and_mRNA_monomers():
Monomer('IkBe', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBe_mRNA')

# def ikbd_and_mRNA_monomers():
Monomer('IkBd', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBd_mRNA')
# Monomer('gene', ['type','state'], {'type': ['a', 'b', 'e', 'd'],'state':['off', 'on']})

# Declares NFkB, IKK1, and IKK2 which having a binding site to IkB and NFkB can transport between the Cytoplasm
# and the Nucleus. IKK1 and IKK2 only exist in the Cytoplasm.

# def nfkb_and_ikk_monomers():
Monomer('NFkB', ['ikb', 'S'], {'S': ['C', 'N']})

Parameter('TNF_0', 1.96e-4) #698 is 30ng/ml of TNF
Parameter('TNFR_0', .01)
Parameter('TRADD_0', 8.3e-4)
Parameter('RIP1_0', 8.3e-4) #47000
Parameter('TRAF_0', 8.3e-4)
# Parameter('cIAP_0', 8.3e-4) #10000
# Parameter('A20_0', .0049) #2256
Parameter('CYLD_0', 0.004) #50000
Parameter('FADD_0', 0.0033)
# Parameter('flip_L_0', 0.004)
Parameter('Lubac_0', 0.003)
Parameter('C8_0', 0.0033) #10000
Parameter('RIP3_0', 0.004) #20000
Parameter('NEMO_0', 0.1) # 1000000
Parameter('MLKLa_0', 0.0034) # 100000

Initial(TNF(brec=None), TNF_0)
Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
# Initial(cIAP(btraf=None), cIAP_0)
# Initial(A20(brip = None), A20_0)
Initial(CYLD(brip=None, btraf = None), CYLD_0)
Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
# Initial(flip_L(bDED=None), flip_L_0)
Initial(LUBAC(brip=None), Lubac_0)
Initial(NEMO(brip = None, btak = None, bikk = None, state = 'I'), NEMO_0)
Initial(C8(bf=None, state='A'), C8_0)
Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)
# #IkBa initial conditions
Parameter('IkBa_0', 0.0025) #IkBa
Initial(IkBa(nfkb = None, S = 'C'), IkBa_0)

Parameter('IkBan_0', 0.0013) #IkBa
Initial(IkBa(nfkb = None, S = 'N'), IkBan_0)

Parameter('IkBa_NFkB_0', 0.0621) #IkBa
Initial(IkBa(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBa_NFkB_0)

Parameter('IkBa_NFkBn_0', 0.0208) #IkBa
Initial(IkBa(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBa_NFkBn_0)

Parameter('IkBat_0', 0.0020) #IkBa
Initial(IkBa_mRNA(), IkBat_0)
#
# #IkBb initial conditions
Parameter('IkBb_0', 0.0044) #IkBb
Initial(IkBb(nfkb = None, S = 'C'), IkBb_0)

Parameter('IkBbn_0', 0.0002) #IkBa
Initial(IkBb(nfkb = None, S = 'N'), IkBbn_0)

Parameter('IkBb_NFkB_0', 0.0237) #IkBa
Initial(IkBb(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBb_NFkB_0)

Parameter('IkBb_NFkBn_0', 0.0016) #IkBb
Initial(IkBb(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBb_NFkBn_0)

Parameter('IkBbt_0', 0.0033) #IkBa
Initial(IkBb_mRNA(), IkBbt_0)
#
# #IkBe initial conditions
Parameter('IkBe_0', 0.0003) #IkBa
Initial(IkBe(nfkb = None, S = 'C'), IkBe_0)

Parameter('IkBen_0', 0.0001) #IkBa
Initial(IkBe(nfkb = None, S = 'N'), IkBen_0)

Parameter('IkBeNFkB_0', 0.0045) #IkBa
Initial(IkBe(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBeNFkB_0)

Parameter('IkBeNFkBn_0', 0.0015) #IkBa
Initial(IkBe(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBeNFkBn_0)

Parameter('IkBet_0', 0.0003) #IkBa
Initial(IkBe_mRNA(), IkBet_0)

#IkBd initial conditions
Parameter('IkBd_0', 0.0003) #IkBa
Initial(IkBd(nfkb = None, S = 'C'), IkBd_0)

Parameter('IkBdn_0', 0.0003) #IkBa
Initial(IkBd(nfkb = None, S = 'N'), IkBdn_0)

Parameter('IkBdNFkB_0', 0.0058) #IkBa
Initial(IkBd(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBdNFkB_0)

Parameter('IkBdNFkBn_0', 0.0039) #IkBa
Initial(IkBd(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBdNFkBn_0)

Parameter('IkBdt_0', 0.0001) #IkBa
Initial(IkBd_mRNA(), IkBdt_0)

Parameter('NFkB_0', 0.000) #Nuclear Factor-kappaB
Initial(NFkB(ikb=None, S='C'), NFkB_0)

Parameter('NFkBn_0', 0.0012) #Nuclear Factor-kappaB in nucleus
Initial(NFkB(ikb=None, S='N'), NFkBn_0)

Parameter('IKK_off_0', 0.0985) #TRADD-TRAF-RIP
Initial(IKK(bind = None, bnemo = None, state = 'I'), IKK_off_0)

Parameter('IKK_0', 0.0002) #TRADD-TRAF-RIP
Initial(IKK(bind = None, bnemo = None,state = 'A'), IKK_0)

Parameter('IKKi_0', 0.0013) #TRADD-TRAF-RIP
Initial(IKK(bind = None, bnemo = None,state = 'AI'), IKKi_0)

#COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 100)
Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 0.75)
Parameter('bind_TNF_TNFR_kf', 1100)
Parameter('bind_TNF_TNFR_kr', 0.021)
Parameter('TNF_deg1', 0.9)
Parameter('bind_TNFRANY_TRADD_kf', 100)
Parameter('bind_TNFRANY_TRADD_kr', 0.75)
Parameter('bind_TNFRANY_RIP1unmod_kf', 100)
Parameter('bind_TNFRANY_RIP1unmod_kr', 0.75)
Parameter('bind_RIP1ANY_TRAFunmod_kf', 100)
Parameter('bind_RIP1ANY_TRAFunmod_kr', 0.75)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf', 100)
Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr', 0.75)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf', 100)
Parameter('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr', 0.75)

Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)

Rule('TNF_deg', TNF(brec = None) >> None, TNF_deg1)


Rule('bind_TNFRANY_TRADD', TNF(brec=1) % TNFR(blig=1, brip=None) + TRADD(brec = None, brip = None, bDD1=None, bDD2=None) <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None), bind_TNFRANY_TRADD_kf, bind_TNFRANY_TRADD_kr)

Rule('bind_TNFRANY_RIP1unmod', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM=None,bMLKL=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod'), bind_TNFRANY_RIP1unmod_kf, bind_TNFRANY_RIP1unmod_kr)

# Rule('Complex_II_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod') <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod'), k_ComplexII_1)
Rule('Complex_I_ubiquitylation1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') + TRAF(brip=None, bciap=None, state='unmod')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, state='unmod'), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)

# Rule('Complex_II_2', TNFR(brip=2) % RIP1(bscf=2, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') >> RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='deub') + TNFR(brip=None), k_ComplexII_2)
Rule('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod', cIAP(btraf=None) + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') <> cIAP(btraf=1) % TRAF(brip=None, bciap=1, bcyld = None, state='unmod'), bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf, bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr)

Rule('bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub', CYLD(brip=None) + TRAF(brip=None, bciap=None, bcyld=None, state='unmod') <> CYLD(brip=1) % TRAF(brip=None, bciap=None, bcyld=1, state='unmod'), bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kf, bind_CYLD_TRAFK63ub_to_CYLDTRAFK63ub_kr)

Rule('Complex_I_ubiquitylation', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld = None, state='unmod') % cIAP(btraf = 5), bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)



Parameter('CompI_UB1', 30)
Parameter('CompI_UB2', 30)
Parameter('bind_RIP1K63ubANY_A20_kf', 500)
Parameter('bind_RIP1K63ubANY_A20_kr', 0.75)
Parameter('k_A20_1', 100)
Parameter('bind_RIP1K63ubANY_CYLD_kf', 500)
Parameter('bind_RIP1K63ubANY_CYLD_kr', 0.75)

Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
     >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)


Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
     bind_RIP1ANY_TRAFunmod_kf, bind_RIP1ANY_TRAFunmod_kr)


#RIP1 K63ub to be deub by A20 or CYLD
#
Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
     bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)

# Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
#      >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None),
#      k_A20_1)


Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
     >>  TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub'), k_A20_1)



Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
     bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)

Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
     >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
     k_A20_1)

#Initiating Necroptosis

Parameter('bind_FADD_proC8_2_kf', 500)
Parameter('bind_FADD_proC8_2_kr', 0.018)
Parameter('bind_FADDANY_flip_L_kf', 500)
Parameter('bind_FADDANY_flip_L_kr', 0.018)
Parameter('bind_FADDANY_proC8_kf', 500)
Parameter('bind_FADDANY_proC8_kr', 0.018)
Parameter('kc_c8_2', 1000)
Parameter('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf', 30)
Parameter('k20', 0.001)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 30)
Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 0.75)
Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 30)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf', 30)
Parameter('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr', 0.75)
Parameter('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc', 0.75)

#RIP1 deub and necrosome formation

Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)

#DOES THIS HAPPEN?!?!??! @TODO
# Rule('bind_FADD_flip_L_2', C8(bf=None, state='I') + flip_L(bDED=None) <> C8(bf=1, state='I') % flip_L(bDED=1), bind_FADD_flip_L_2_kf, bind_FADD_flip_L_2_kr)
# Rule('C8_activation1', C8(bf=1, state='I') % flip_L(bDED=1) >> flip_L(bDED=None) + C8(bf=None, state='A'), kc_c8_1)


Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, state='A')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='A'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)


Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, state='A') + flip_L(bDED=None)
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)


TNF
Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) + RIP3(bRHIM=None, bDD = None, state='unmod')
     <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_proC8_kf, bind_FADDANY_proC8_kr)


# Rule('bind_C8A_CYLDU_to_C8ACYLDU', C8(bf=None, state='A') +  <> C8(bf=1, state='A') % CYLD(btraf=1, state='U'), bind_C8A_CYLDU_to_C8ACYLDU_kf, bind_C8A_CYLDU_to_C8ACYLDU_kr)

Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = 4) % C8(bf=2, state='A') % flip_L(bDED=4) % RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)


Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod_kf)

Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
     >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), k20)

Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
     <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)

Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = 1, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') % MLKL(bRHIM=1, state='unmod')
     >>  MLKL(bRHIM=None, state='active') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)

Rule('bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod', MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='unmod') <> MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod'), bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kf, bind_MLKLactive_MLKLunmod_to_MLKLactiveMLKLunmod_kr)

Rule('catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive', MLKL(bRHIM=1, state='active') % MLKL(bRHIM=1, state='unmod') >> MLKL(bRHIM=None, state='active') + MLKL(bRHIM=None, state='active'), catalyze_MLKLactiveMLKLunmod_to_MLKLactive_MLKLactive_kc)




#Initiating Survival
Parameter('bind_RIP1K63ubANY_NEMO_kf', 500)
Parameter('bind_RIP1K63ubANY_NEMO_kr', 0.25)
Parameter('bind_NEMORIP1TAK1_IKKinactive_kf', 500)
Parameter('bind_NEMORIP1TAK1_IKKinactive_kr', 0.25)
Parameter('bind_RIP1K63ubANY_TAK1_kf', 500)
Parameter('bind_RIP1K63ubANY_TAK1_kr', 0.25)
Parameter('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc', 0.1)
Parameter('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kf', 500)
Parameter('bind_TAK1ANY_MAPKinactive_to_TAK1ANYMAPKinactive_kr', 0.25)
Parameter('catalyze_TAK1ANYMAPKinactive_to_TAK1ANY_MAPKactive_kc', .75)
Parameter('IKKi_f_IKK', 5e-5)
Parameter('IKKa_r_IKK', 0.02)
Parameter('IKKi_f_IKKaIKKK', 520.0)
Parameter('IKKa_f_IKKai', 0.15)
Parameter('IKKai_f_IKKi', 0.02)


Rule('bind_RIP1K63ubANY_NEMO', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + NEMO(brip=None, state = 'I')
     <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) +  NEMO(brip=None, state = 'A'), bind_RIP1K63ubANY_NEMO_kf, bind_RIP1K63ubANY_NEMO_kr)



# Rule('bind_RIP1K63ubANY_TAK1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) + TAK1(brip=None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8), bind_RIP1K63ubANY_TAK1_kf, bind_RIP1K63ubANY_TAK1_kr)


# Rule('bind_NEMORIP1TAK1_IKKinactive', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) + IKK(bind=None, bnemo = None, state = 'I')
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = 9,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) % IKK(bind=9, bnemo = None, state = 'I'), bind_NEMORIP1TAK1_IKKinactive_kf, bind_NEMORIP1TAK1_IKKinactive_kr)



# Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = 9,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) % IKK(bind=9, bnemo = None, state = 'I')
#      >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) + IKK(bind=None, bnemo = None, state = 'A'), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)


# Rule('catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=8,bDD = 9,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % NEMO(brip=7) % TAK1(brip=8) % IKK(bind=9, bnemo = None, state = 'I')
#      >> NEMO(brip=None) + IKK(bind=None, bnemo = None, state = 'A'), catalyze_NEMORIP1TAK1IKKinactive_to_NEMORIP1TAK1_IKKactive_kc)

Rule('IKKi_IKKa', IKK(bind=None, bnemo = None, state = 'I') <> IKK(bind=None, bnemo = None, state = 'A'), IKKi_f_IKK,IKKa_r_IKK)
Rule('IKKi_IKKaIKKK', IKK(bind=None, bnemo = None, state = 'I') + NEMO(brip=None, state = 'A') >> IKK(bind=None, bnemo = None, state = 'A') + NEMO(brip=None, state = 'A'), IKKi_f_IKKaIKKK)
Rule('IKKi_IKKai', IKK(bind=None, bnemo = None, state = 'A') >> IKK(bind=None, bnemo = None, state = 'AI'), IKKa_f_IKKai)
Rule('IKKai_IKKi', IKK(bind=None, bnemo = None, state = 'AI') >> IKK(bind=None, bnemo = None, state = 'I'), IKKai_f_IKKi)


#
# def  ikb_mrna_to_ikb():
Parameter('psynth_a', 7e-5)
Parameter('psynth_b', 1e-5)
Parameter('psynth_e', 1e-6)
Parameter('psynth_d', 1e-7)
Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
Rule('d_synth', None >> IkBd_mRNA(), psynth_d)
#
# Rule('a_synth', NFkB(ikb=None, S='N') >> IkBa_mRNA() + NFkB(ikb=None, S='N'), psynth_a)
# Rule('b_synth', NFkB(ikb=None, S='N') >> IkBb_mRNA() + NFkB(ikb=None, S='N'), psynth_b)
# Rule('e_synth', NFkB(ikb=None, S='N') >> IkBe_mRNA() + NFkB(ikb=None, S='N'), psynth_e)
# Rule('d_synth', NFkB(ikb=None, S='N') >> IkBd_mRNA() + NFkB(ikb=None, S='N'), psynth_d)

#
# Parameter('a_delay', 0.1)
# Parameter('be_delay', 0.027)
# Parameter('d_delay', 0.011)
#
# Rule('genea_off_on', gene(type = 'a', state = 'off') >> gene(type = 'a',state = 'on'), a_delay)
# Rule('geneb_off_on', gene(type = 'b', state = 'off') >> gene(type = 'b',state = 'on'), be_delay)
# Rule('genee_off_on', gene(type = 'e', state = 'off') >> gene(type = 'e',state = 'on'), be_delay)
# Rule('gened_off_on', gene(type = 'd', state = 'off') >> gene(type = 'd',state = 'on'), d_delay)


Parameter('hill', 3)
Parameter('a', 8)
Parameter('b', 0.02)
Parameter('e', 0.3)
Parameter('d', 0.025)

Observable('NFkBn_free', NFkB(ikb=None, S='N'))

Expression('a_NFkBn', a*(NFkBn_free)**(hill))
Expression('b_NFkBn', b*(NFkBn_free)**(hill))
Expression('e_NFkBn', e*(NFkBn_free)**(hill))
Expression('d_NFkBn', d*(NFkBn_free)**(hill))
#
# Rule('an_mRNA', gene(type = 'a',state = 'on') >> IkBa_mRNA() + gene(type = 'a',state = 'on'), a_NFkBn)
# Rule('bn_mRNA', gene(type = 'b',state = 'on') >> IkBb_mRNA() + gene(type = 'b',state = 'on'), a_NFkBn)
# Rule('en_mRNA', gene(type = 'e',state = 'on') >> IkBe_mRNA() + gene(type = 'e',state = 'on'), a_NFkBn)
# Rule('dn_mRNA', gene(type = 'd',state = 'on') >> IkBd_mRNA() + gene(type = 'd',state = 'on'), a_NFkBn)


Rule('an_mRNA', None >> IkBa_mRNA(), a_NFkBn)
Rule('bn_mRNA', None >> IkBb_mRNA(), b_NFkBn)
Rule('en_mRNA', None >> IkBe_mRNA(), e_NFkBn)
Rule('dn_mRNA', None >> IkBd_mRNA(), d_NFkBn)


# IkB mRNA and protein synthesis reactions

Parameter('mRNA_a', 0.035)
Rule('a_mRNA', IkBa_mRNA() >> None, mRNA_a)


Parameter('mRNA_b', 3e-3)
Parameter('mRNA_e', 4e-3)
Parameter('mRNA_d', 2e-3)

Rule('b_mRNA', IkBb_mRNA() >> None, mRNA_b)
Rule('e_mRNA', IkBe_mRNA() >> None, mRNA_e)
Rule('d_mRNA', IkBd_mRNA() >> None, mRNA_d)


# Parameter('synth', 0.2448)
# Rule('a_psynth', None >> IkBa(ikk=None, nfkb=None, S='C'), synth)
# Rule('b_psynth', None >> IkBb(ikk=None, nfkb=None, S='C'), synth)
# Rule('e_psynth', None >> IkBe(ikk=None, nfkb=None, S='C'), synth)
# Rule('d_psynth', None >> IkBd(ikk=None, nfkb=None, S='C'), synth)

Parameter('syntha', 0.2448)
Rule('a_psynth', IkBa_mRNA() >> IkBa(nfkb=None, S='C') + IkBa_mRNA(), syntha)

Parameter('synthb', 0.2448)
Rule('b_psynth', IkBb_mRNA() >> IkBb(nfkb=None, S='C') + IkBb_mRNA(), synthb)

Parameter('synthe', 0.2448)
Rule('e_psynth', IkBe_mRNA() >> IkBe(nfkb=None, S='C') + IkBe_mRNA(), synthe)

Parameter('synthd', 0.2448)
Rule('d_psynth', IkBd_mRNA() >> IkBd(nfkb=None, S='C') + IkBd_mRNA(), synthd)

#IkB(a,b,e) association and dissociation from IKK2 and IkBd association and dissociation from IKK2
# def ikb_assoc_diss_nfkb():
Parameter('IkB_IKKf', 30)
Parameter('IkB_IKKr', 6e-5)
Rule('an_adc', IkBa(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('bn_adc', IkBb(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('en_adc', IkBe(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('dn_adc', IkBd(nfkb=None, S='C') + NFkB(ikb=None, S='C') <> IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)

Rule('an_adn', IkBa(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('bn_adn', IkBb(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('en_adn', IkBe(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('dn_adn', IkBd(nfkb=None, S='N') + NFkB(ikb=None, S='N') <> IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)


# #IkB and NFkB cellular localization reactions
# def ikb_nfkb_localization():
Parameter('af', 0.09)
Parameter('bf', 0.009)
Parameter('ef', 0.045)
Parameter('df', 0.045)
Parameter('ancf', 0.012)
Parameter('bncf', 0.012)
Parameter('encf', 0.012)
Parameter('dncf', 0.012)
Rule('a_nc', IkBa(nfkb=None, S='C') <> IkBa(nfkb=None, S='N'), af, ancf)
Rule('b_nc', IkBb(nfkb=None, S='C') <> IkBb(nfkb=None, S='N'), bf, bncf)
Rule('e_nc', IkBe(nfkb=None, S='C') <> IkBe(nfkb=None, S='N'), ef, encf)
Rule('d_nc', IkBd(nfkb=None, S='C') <> IkBd(nfkb=None, S='N'), ef, dncf)

Parameter('anf', 0.276)
Parameter('bnf', 0.0276)
Parameter('enf', 0.138)
Parameter('dnf', 0.276)
Parameter('anr', 0.828)
Parameter('bnr', 0.414)
Parameter('enr', 0.414)
Parameter('dnr', 0.414)
Rule('an_nc', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'), anf, anr)
Rule('bn_nc', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'), bnf, bnr)
Rule('en_nc', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'), enf, enr)
Rule('dn_nc', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') <> IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'), dnf, dnr)

Parameter('nf', 5.4)
Parameter('nr', 0.0048)
Rule('n_nc', NFkB(ikb=None, S='C') <> NFkB(ikb=None, S='N'), nf, nr)

# IkB Protein Degradation Reactions
# def ikb_deg_reactions():
Parameter('a_dc', 0.12)
Parameter('b_dc', 0.18)
Parameter('e_dc', 0.18)
Parameter('d_dc', 0.0014)
Parameter('a_dn', 0.12)
Parameter('b_dn', 0.18)
Parameter('e_dn', 0.18)
Parameter('d_dn', 0.0014)

Rule('ad_c', IkBa(nfkb=None, S='C') >> None, a_dc)
Rule('bd_c', IkBb(nfkb=None, S='C') >> None, b_dc)
Rule('ed_c', IkBe(nfkb=None, S='C') >> None, e_dc)
Rule('dd_c', IkBd(nfkb=None, S='C') >> None, d_dc)

Rule('ad_n', IkBa(nfkb=None, S='N') >> None, a_dn)
Rule('bd_n', IkBb(nfkb=None, S='N') >> None, b_dn)
Rule('ed_n', IkBe(nfkb=None, S='N') >> None, e_dn)
Rule('dd_n', IkBd(nfkb=None, S='N') >> None, d_dn)

Parameter('c_bn', 0.00006)
Parameter('n_bn', 0.00006)

Rule('an_c', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('bn_c', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('en_c', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('dn_c', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)

Rule('an_n', IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('bn_n', IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('en_n', IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('dn_n', IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)

#Added Reactions

#IKK-mediated IkB degradation reactions
# def ikb_ikk_mediated_deg():
# Parameter('a_f_deg', 0.36)
# Parameter('b_f_deg', 0.12)
# Parameter('e_f_deg', 0.18)
# Parameter('d_f_deg', 0.18)
#
# Parameter('and_c_n', 0.36)
# Parameter('bnd_c_n', 0.12)
# Parameter('end_c_n', 0.18)
# Parameter('dnd_c_n', 0.18)
#
# # IkBa >> None
# # IkBb >> None
# # IkBe >> None
# # IkBd >> None
#
# # IkBa : NFkB >> NFkB
# # IkBb : NFkB >> NFkB
# # IkBe : NFkB >> NFkB
# # IkBd : NFkB >> NFkB
#
# Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, a_f_deg)
# Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, b_f_deg)
# Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, e_f_deg)
# Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)
#
# Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), and_c_n)
# Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), bnd_c_n)
# Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), end_c_n)
# Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)

#
Observable('IKKa_obs', IKK(bind=None, bnemo = None, state = 'A'))

Parameter('a_f_deg', 3.6)
Parameter('b_f_deg', 1.2)
Parameter('e_f_deg', 1.8)
Parameter('d_f_deg', 0.0018)
#
# Parameter('and_c_n', 0.36)
# Parameter('bnd_c_n', 0.12)
# Parameter('end_c_n', 0.18)
Parameter('dnd_c_n', 0.0018)
# # Parameter('IKK_boom', 0)
#
# # Expression('IKK_ikba_flux', IKK_boom*a_f_deg)
# # Expression('IKK_ikbb_flux', IKK_boom*b_f_deg)
# # Expression('IKK_ikbe_flux', IKK_boom*e_f_deg)
#
#
#
#
Expression('IKK_ikba_flux', IKKa_obs*a_f_deg)
Expression('IKK_ikbb_flux', IKKa_obs*b_f_deg)
Expression('IKK_ikbe_flux', IKKa_obs*e_f_deg)

#
# # IkBa >> None
# # IkBb >> None
# # IkBe >> None
# # IkBd >> None
#
# # IkBa : NFkB >> NFkB
# # IkBb : NFkB >> NFkB
# # IkBe : NFkB >> NFkB
# # IkBd : NFkB >> NFkB
#
Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, IKK_ikba_flux)
Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, IKK_ikbb_flux)
Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, IKK_ikbe_flux)
Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)
#
Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikba_flux)
Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikbb_flux)
Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikbe_flux)
Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)
#
#
# #A20 mRNA and Protein Synthesis and Degradation Reactions
# # def a20_mrna_to_a20():
Parameter('A20_mRNA', 2e-6)
Parameter('A20n', 0.4)
Parameter('A20_mRNA_c_deg', 0.035)
Parameter('a1d_c_deg', 0.36)
Parameter('A20_synth', 0.25)
Parameter('A20_deg', 0.0029)
Parameter('newhill', 2)

Observable('obs_A20t', A20t())

Expression('A20t_NFkBn', A20n*(NFkBn_free)**(hill))
Expression('A20_synthesis', A20_synth*obs_A20t)

Rule('A20t_synth', None >> A20t(), A20_mRNA)
Rule('A20t_mediated_nfkbn', None >> A20t(), A20t_NFkBn)
Rule('A20t_deg', A20t() >> None, A20_mRNA_c_deg)

Rule('synth_A20', None >> A20(brip = None), A20_synthesis)
Rule('deg_A20', A20(brip = None) >> None, A20_deg)

Expression('flipl_NFkBn', A20n*(NFkBn_free)**(hill))
Expression('ciap_NFkBn', A20n*(NFkBn_free)**(newhill))
Rule('synth_flipl', None >> flip_L(bDED = None), flipl_NFkBn)
Rule('synth_ciap', None >> cIAP(btraf = None), ciap_NFkBn)



Observable('Flip_obs', flip_L(bDED = None))
Observable('cIAP_obs', cIAP(btraf = None))
Observable('RIP1_obs',RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, state='unmod'))
Observable('RIP3_obs', RIP3(bRHIM=None, bDD = None, state='unmod'))
Observable('MLKL_obs', MLKL(bRHIM=None, state='unmod'))
Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
Observable('RIP13_obs', RIP1(bDD=ANY, bRHIM=1, state='unmod') % RIP3(bRHIM=1, bDD = None, state='unmod'))
Observable('RIP13po4_obs', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'))
Observable('TNF_obs', TNF(brec = ANY))
Observable('RIPk63_obs',  RIP1(bscf=None, bub1=None, state='K63ub'))
Observable('TNFR_TRADD', TNFR(blig=None, brip=1) % TRADD(brec=1, brip=None))
Observable('CI_k63_obs', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) %
           RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5))
Observable('CI', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec=2, brip=3) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None, state='K63ub') % TRAF(brip=4, bciap=5, state='unmod') % cIAP(btraf = 5))
Observable('NFkBn_obs', NFkB(ikb=None, S='N'))
Observable('NFkB_obs', NFkB(ikb=None, S='C'))
Observable('IkBa_obs', IkBa(nfkb = None, S='C'))
Observable('IkBan_obs', IkBa(nfkb = None, S='N'))
Observable('IkBaNFkB_obs', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBaNFkBn_obs', IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

Observable('IkBb_obs', IkBb(nfkb = None, S='C'))
Observable('IkBbn_obs', IkBb(nfkb = None, S='N'))
Observable('IkBbNFkB_obs', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBbNFkBn_obs', IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'))

Observable('IkBe_obs', IkBe(nfkb = None, S='C'))
Observable('IkBen_obs', IkBe(nfkb = None, S='N'))
Observable('IkBeNFkB_obs', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBeNFkBn_obs', IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'))


Observable('IkBd_obs', IkBd(nfkb = None, S='C'))
Observable('IkBdn_obs', IkBd(nfkb = None, S='N'))
Observable('IkBdNFkB_obs', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C'))
Observable('IkBdNFkBn_obs', IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'))
Observable('A20_obs', A20(brip = None))

tspan = np.linspace(0, 1440, 14401)

sim1 = ScipyOdeSimulator(model, tspan)
sim2 =ScipyOdeSimulator(model, tspan)
sim3 = ScipyOdeSimulator(model, tspan)
sim4 =ScipyOdeSimulator(model, tspan)

L1 = sim1.run(param_values={'TNF_0':1.96e-5})
L2 = sim2.run(param_values={'TNF_0':1.96e-4})
L3 = sim3.run(param_values={'TNF_0':1.96e-3})
L4 = sim4.run(param_values={'TNF_0':1.96e-2})

print(len(model.rules))
print(len(model.odes))
print(len(model.parameters))


# tnf_dose = [1.96e-5, 1.96e-4, 1.96e-3, 1.96e-2]
#
# tspan = np.linspace(0, 1440, 14401)
# # print(len(tspan))
# # print(tspan)
# sim = ScipyOdeSimulator(model, tspan = tspan)
# # sim2 = ScipyOdeSimulator(model, tspan = tspan)
# # simulation_result = sim.run(param_values={'TNF_0': 1.96e-5})
# # simulation_result2 = sim2.run(param_values={'TNF_0': 1.96e-4})
# simulation_result = sim.run(param_values={'TNF_0': [1.96e-5, 1.96e-4, 1.96e-3, 1.96e-2]})
#
# df = simulation_result.dataframe
# # df2 = simulation_result2.dataframe
#
# print('df1')
# print(df.loc[0]['TNF_obs'].iloc[:])
# print('df2')
# print(df.loc[3]['TNF_obs'].iloc[:])
# quit()
# print(df.loc[2]['TNF_obs'].iloc[-1])
# print(df.loc[3]['TNF_obs'].iloc[-1])
# quit()

#
# for i,sp in enumerate(model.species):
#     print i,":",sp
#
# print(len(model.species))

print('nfkbn1')
print(L1.observables['Flip_obs'])
print('nfkbn2')
print(L2.observables['cIAP_obs'])
print('nfkbn3')
print(L3.observables['NFkBn_obs'])
print('nfkbn4')
print(L4.observables['NFkBn_obs'])
# quit()

plt.figure(figsize = (15,10))
# plt.figure()
plt.subplot(231)
plt.plot(tspan/60, L1.observables['TNF_obs'],label = 'TNF.1')
plt.plot(tspan/60, L2.observables['TNF_obs'],label = 'TNF1')
plt.plot(tspan/60, L3.observables['TNF_obs'],label = 'TNF10')
plt.plot(tspan/60, L4.observables['TNF_obs'],label = 'TNF100')

# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(232)
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, L1.observables['CI_k63_obs'],label = 'CI_k63.1')
# plt.plot(tspan/60, L2.observables['CI_k63_obs'],label = 'CI_k631')
# plt.plot(tspan/60, L3.observables['CI_k63_obs'],label = 'CI_k6310')
# plt.plot(tspan/60, L4.observables['CI_k63_obs'],label = 'CI_k63100')
plt.plot(tspan/60, L1.observables['MLKLa_obs'],label = 'MLKLa.1')
plt.plot(tspan/60, L2.observables['MLKLa_obs'],label = 'MLKLa1')
plt.plot(tspan/60, L3.observables['MLKLa_obs'],label = 'MLKLa10')
plt.plot(tspan/60, L4.observables['MLKLa_obs'],label = 'MLKLa100')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['CI_k63_obs'],label = 'CI_k63')
# plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], color = 'r', label = 'TNFR_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(233)
plt.plot(tspan/60, L1.observables['A20_obs'],label = 'A20.1')
plt.plot(tspan/60, L2.observables['A20_obs'],label = 'A201')
plt.plot(tspan/60, L3.observables['A20_obs'],label = 'A2010')
plt.plot(tspan/60, L4.observables['A20_obs'],label = 'A20100') 
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
plt.plot(tspan/60, L1.observables['IKKa_obs'],label = 'IKKa.1')
plt.plot(tspan/60, L2.observables['IKKa_obs'],label = 'IKKa1')
plt.plot(tspan/60, L3.observables['IKKa_obs'],label = 'IKKa10')
plt.plot(tspan/60, L4.observables['IKKa_obs'],label = 'IKKa100')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(235)
plt.plot(tspan/60, L1.observables['NFkBn_obs'],label = 'NFkBn.1')
plt.plot(tspan/60, L2.observables['NFkBn_obs'],label = 'NFkBn1')
plt.plot(tspan/60, L3.observables['NFkBn_obs'],label = 'NFkBn10')
plt.plot(tspan/60, L4.observables['NFkBn_obs'],label = 'NFkBn100')
# plt.plot(tspan/60, L1.observables['NFkB_obs'],label = 'NFkB.1')
# plt.plot(tspan/60, L2.observables['NFkB_obs'],label = 'NFkB1')
# plt.plot(tspan/60, L3.observables['NFkB_obs'],label = 'NFkB10')
# plt.plot(tspan/60, L4.observables['NFkB_obs'],label = 'NFkB100')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(236)
plt.plot(tspan/60, L1.observables['IkBa_obs'],label = 'IkBa.1')
plt.plot(tspan/60, L2.observables['IkBa_obs'],label = 'IkBa1')
plt.plot(tspan/60, L3.observables['IkBa_obs'],label = 'IkBa10')
plt.plot(tspan/60, L4.observables['IkBa_obs'],label = 'IkBa100')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)


plt.figure(figsize = (20,5))
# plt.figure()
plt.subplot(141)
plt.plot(tspan/60, L1.observables['IkBa_obs'],label = 'IkBa.1')
plt.plot(tspan/60, L2.observables['IkBa_obs'],label = 'IkBa1')
plt.plot(tspan/60, L3.observables['IkBa_obs'],label = 'IkBa10')
plt.plot(tspan/60, L4.observables['IkBa_obs'],label = 'IkBa100')

# plt.plot(tspan/60, simulation_result.observables['TNF_obs'], color = 'r', label = 'TNF_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(142)
plt.plot(tspan/60, L1.observables['IkBb_obs'],label = 'IkBb.1')
plt.plot(tspan/60, L2.observables['IkBb_obs'],label = 'IkBb1')
plt.plot(tspan/60, L3.observables['IkBb_obs'],label = 'IkBb10')
plt.plot(tspan/60, L4.observables['IkBb_obs'],label = 'IkBb100')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(143)
plt.plot(tspan/60, L1.observables['IkBe_obs'],label = 'IkBe.1')
plt.plot(tspan/60, L2.observables['IkBe_obs'],label = 'IkBe1')
plt.plot(tspan/60, L3.observables['IkBe_obs'],label = 'IkBe10')
plt.plot(tspan/60, L4.observables['IkBe_obs'],label = 'IkBe100')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(144)
plt.plot(tspan/60, L1.observables['IkBd_obs'],label = 'IkBd.1')
plt.plot(tspan/60, L2.observables['IkBd_obs'],label = 'IkBd1')
plt.plot(tspan/60, L3.observables['IkBd_obs'],label = 'IkBd10')
plt.plot(tspan/60, L4.observables['IkBd_obs'],label = 'IkBd100')
# plt.plot(tspan/60, simulation_result.observables['IKKa_obs'], color = 'r', label = 'IKKa_mat')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.tight_layout()
plt.show()