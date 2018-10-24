from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
from pysb.simulator.bng import BngSimulator

Model()

model.enable_synth_deg()

#
# # file = open("output_data.txt", 'r')
# # print(file.read())
# #
# # initials_1 = np.array([p.value for p in list(model.initial_conditions)])
# # print(initials_1)
# # quit()
# # initials = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000]
# #
# # for i in range(0,101):
# #     x = file[0]
# #     print(x)
# # quit()
# #     # nvp = initials + file[i]
# pre_cal = [  1.00000000e-06,   1.00000000e-03,   1.00000000e-03 ,  1.00000000e-06,
#    1.00000000e-03,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-06,
#    1.00000000e-03,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
#    1.00000000e-06,   1.00000000e-03,   1.00000000e-06 ,  1.00000000e-03,
#    1.00000000e-01,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
#    1.00000000e-01,   3.11000000e-07,   3.27000000e-06 ,  1.80000000e-02,
#    3.27000000e-06,   1.80000000e-02,   3.27000000e-02 ,  1.80000000e-02,
#    1.00000000e-01,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
#    1.00000000e-02,   1.00000000e-03,   1.00000000e-03 ,  1.00000000e-06,
#    1.00000000e+00]
# params_fst_pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# 3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
# 4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
# 4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
# 3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
# 1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
# 3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
# 5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
# 6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
# 1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
# 9.631340201588598493e-01
# ]
# x = range(1,38)
# params_pf_un = [
# 7.244683452113300153e-07,
# 9.517530155749114040e-02,
# 8.637265039628022611e-02,
# 3.193901331928525024e-05,
# 1.258406237885183361e-02,
# 7.575966125090899991e-05,
# 8.364351479918887000e-02,
# 7.436459484254262205e-05,
# 7.768251964618522187e-02,
# 9.819940099507889082e-06,
# 4.836763065283320590e-02,
# 5.181076854275994403e-01,
# 9.056610726977585267e-05,
# 2.338202234638503667e-03,
# 2.117046267096549337e-05,
# 8.693240554107009577e-02,
# 4.110127984888350738e+00,
# 7.183481154705118965e-05,
# 7.015524132740647012e-02,
# 4.116603440959405447e+00,
# 9.694023692347162324e+00,
# 2.218206981846405481e-05,
# 9.837233297311752931e-05,
# 1.188372739414592205e+00,
# 1.004157754027084113e-04,
# 9.913660402796706794e-01, 1.965754591697516096e+00, 1.769000979470284785e+00, 1.129174042737101313e+00,
# 9.269724481355695611e-05, 9.688283102144903958e-02, 6.520509822693891344e+00, 4.227833110161860475e-01,
# 3.229777219239124333e-03, 1.567066653485880284e-02, 4.805448621305080279e-05, 1.389352872855412979e+00
# ]
#
# params_fst = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
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
#
# # plt.figure(figsize=(20,10))
# # # plt.scatter(x, params, label = 'PFu', color = 'blue', lw =1)
# # plt.scatter(x, params_fst, label = 'FST-PSO', color = 'red', lw =1)
# # plt.scatter(x, params_fst_pf, label ='FP-PF', color = 'orange', lw = 1)
# # plt.scatter(x, pre_cal, label ='preCal', color = 'green', lw = 1)
# # plt.xlabel('kinetic parameters', fontsize=25)
# # plt.ylabel('parameter value', fontsize=25)
# # plt.xticks(x)
# # plt.ylim(ymax = 1)
# # # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.legend(loc = 'best',prop={'size': 20})
# # # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# # # plt.show()
# # # quit()
# #
# # params = [
# # 2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # 1.223347573677009020e-07,
# # 1.847763641629712739e-02,
# # 3.000388537507257813e-02,
# # 2.001077439052020398e-05,
# # 1.219295486333531707e-02,
# # 4.943717334826855457e-05,
# # 2.624725680910927617e-02,
# # 1.991351857905033468e-05,
# # 4.488812697741692559e-03,
# # 1.041577480241394007e-05,
# # 3.473936580568400597e-03,
# # 1.072217409826840662e+00,
# # 1.506747791022632791e-05,
# # 4.360182640163040406e-03,
# # 1.792572991490754334e-05,
# # 2.642504878427706302e-02,
# # 3.558771687312299870e-01,
# # 1.106213891407066866e-06,
# # 2.101813258955530922e-02,
# # 4.020551464879950743e+00,
# # 2.458525808081724495e+00,
# # 3.153989353897445407e-06,
# # 1.676074001871189672e-04,
# # 1.631660435291963120e-03,
# # 3.264925281087857373e-05,
# # 1.555466530711384077e-01,
# # 4.656082056292588645e-01,
# # 7.442659427818243412e-01,
# # 1.663505513158469062e+00,
# # 2.393601347545811427e-05,
# # 5.997125487854754709e-03,
# # 5.270384065069513291e-01,
# # 7.442196406911848194e-02,
# # 4.354575919268037498e-02,
# # 2.466300351383134983e-04,
# # 2.046932131600806510e-05,
# # 2.107029175494326978e+01]
# #
# # new = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # 5.090352718964910340e-05,
# # 8.846681140663302523e-02,
# # 7.688925127286951045e-02,
# # 4.718531991263202608e-05,
# # 2.480645592309011632e-02,
# # 2.634730892196591635e-05,
# # 7.994366563229754474e-02,
# # 8.004684705208197151e-05,
# # 2.787447426045872381e-02,
# # 2.528487936602940290e-05,
# # 5.965677558019213955e-02,
# # 2.605672790260717964e+00,
# # 8.794418149861516561e-05,
# # 7.695471791403299400e-02,
# # 7.898521723337450790e-05,
# # 4.694350552814404581e-02,
# # 6.550932796050894069e+00,
# # 8.186365727139884572e-05,
# # 6.872144915868418080e-02,
# # 1.112414613253134510e+00,
# # 1.497407583427553535e+00,
# # 2.679386619274433570e-06,
# # 1.252840334405051683e-04,
# # 2.142769060164864958e-01,
# # 7.132599412461731482e-05,
# # 1.733224145970182262e+00,
# # 2.345217685036323108e+00,
# # 6.735732655693388304e-01,
# # 5.176940275194190200e+00,
# # 5.364253277853003110e-06,
# # 6.539060406153636429e-02,
# # 7.377928959306133017e+00,
# # 7.512669157064747472e-01,
# # 5.520988157998800439e-02,
# # 8.772775867498763813e-02,
# # 1.123578955588378102e-05,
# # 1.062391316471664356e-02,
# #
# # ]

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

Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))

fstpso = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
3.304257485026848768e-05,
9.791215971368776028e-03,
6.110068937548310437e-03,
4.319219461882335177e-05,
4.212644667502709987e-03,
1.164332269398710173e-05,
2.404257105292715788e-02,
3.311085510054400207e-05,
4.280399320119887552e-02,
2.645814637693955063e-05,
1.437707485722601597e-02,
2.303744002126363599e-01,
2.980688423948379739e-05,
4.879773212139151134e-02,
1.121503480627013705e-05,
1.866712857727401229e-03,
7.572177664736867708e-01,
1.591282730353343167e-05,
3.897146248650381478e-02,
3.076363174012029411e+00,
3.734859661130401243e+00,
3.216200470463125876e-06,
8.782429548444865440e-05,
2.906341314225960662e-02,
5.663104188870970508e-05,
2.110469405515222677e-02,
1.294086380531199176e-01,
3.127598126760888220e-01,
4.298489868360909627e-01,
2.332910188537793332e-06,
7.077504621536276526e-03,
6.294061533948092091e-01,
6.419313355304218094e-02,
8.584653667640911989e-04,
8.160445062706172072e-05,
4.354383618147421691e-06,
4.278903092658225660e+00
]


pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
     3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
     4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
     4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
     3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
     1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
     3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
     5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
     6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
     1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
     9.631340201588598493e-01]

tnf = [2326, 233, 23, 2]
lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
# tspan = np.linspace(0, 1440, 1440)
tspan = np.linspace(0, 480, 481)
sim = BngSimulator(model, tspan=tspan)
# result = sim.run(method='ode', param_values=params_fst, initials={TNF(brec=None): tnf})
res = sim.run(method='ode', param_values=fstpso)
# result2 = sim.run(method='ode', param_values=params_ode_new, initials={TNF(brec=None): tnf})
# df = result.dataframe
df1 = res.dataframe

# mlkl = [0,   170, 900, 4880,    9940,    10000]
# y = [0, 1, 2, 4, 6, 8]
color = ['b', 'g', 'm', 'r']

plt.figure()
# for n in range(0,4):
    # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
plt.plot(tspan, df1['MLKLa_obs'][:])
    # plt.plot(tspan, df1.loc[n]['MLKLa_obs'].iloc[:],'--',c = color[n],lw = 1.5) #fppf
plt.xlabel('Time [minutes]', fontsize=16)
plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.ylim(ymax = 10000)
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
plt.show()




from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
from pysb.simulator.bng import BngSimulator

Model()

model.enable_synth_deg()


# # file = open("output_data.txt", 'r')
# # print(file.read())
# #
# # initials_1 = np.array([p.value for p in list(model.initial_conditions)])
# # print(initials_1)
# # quit()
# # initials = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000]
# #
# # for i in range(0,101):
# #     x = file[0]
# #     print(x)
# # quit()
# #     # nvp = initials + file[i]
# pre_cal = [  1.00000000e-06,   1.00000000e-03,   1.00000000e-03 ,  1.00000000e-06,
#    1.00000000e-03,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-06,
#    1.00000000e-03,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
#    1.00000000e-06,   1.00000000e-03,   1.00000000e-06 ,  1.00000000e-03,
#    1.00000000e-01,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
#    1.00000000e-01,   3.11000000e-07,   3.27000000e-06 ,  1.80000000e-02,
#    3.27000000e-06,   1.80000000e-02,   3.27000000e-02 ,  1.80000000e-02,
#    1.00000000e-01,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
#    1.00000000e-02,   1.00000000e-03,   1.00000000e-03 ,  1.00000000e-06,
#    1.00000000e+00]
# params_fst_pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# 3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
# 4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
# 4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
# 3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
# 1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
# 3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
# 5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
# 6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
# 1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
# 9.631340201588598493e-01
# ]
# x = range(1,38)
# params_pf_un = [
# 7.244683452113300153e-07,
# 9.517530155749114040e-02,
# 8.637265039628022611e-02,
# 3.193901331928525024e-05,
# 1.258406237885183361e-02,
# 7.575966125090899991e-05,
# 8.364351479918887000e-02,
# 7.436459484254262205e-05,
# 7.768251964618522187e-02,
# 9.819940099507889082e-06,
# 4.836763065283320590e-02,
# 5.181076854275994403e-01,
# 9.056610726977585267e-05,
# 2.338202234638503667e-03,
# 2.117046267096549337e-05,
# 8.693240554107009577e-02,
# 4.110127984888350738e+00,
# 7.183481154705118965e-05,
# 7.015524132740647012e-02,
# 4.116603440959405447e+00,
# 9.694023692347162324e+00,
# 2.218206981846405481e-05,
# 9.837233297311752931e-05,
# 1.188372739414592205e+00,
# 1.004157754027084113e-04,
# 9.913660402796706794e-01, 1.965754591697516096e+00, 1.769000979470284785e+00, 1.129174042737101313e+00,
# 9.269724481355695611e-05, 9.688283102144903958e-02, 6.520509822693891344e+00, 4.227833110161860475e-01,
# 3.229777219239124333e-03, 1.567066653485880284e-02, 4.805448621305080279e-05, 1.389352872855412979e+00
# ]
#
# params_fst = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
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
#
# # plt.figure(figsize=(20,10))
# # # plt.scatter(x, params, label = 'PFu', color = 'blue', lw =1)
# # plt.scatter(x, params_fst, label = 'FST-PSO', color = 'red', lw =1)
# # plt.scatter(x, params_fst_pf, label ='FP-PF', color = 'orange', lw = 1)
# # plt.scatter(x, pre_cal, label ='preCal', color = 'green', lw = 1)
# # plt.xlabel('kinetic parameters', fontsize=25)
# # plt.ylabel('parameter value', fontsize=25)
# # plt.xticks(x)
# # plt.ylim(ymax = 1)
# # # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.legend(loc = 'best',prop={'size': 20})
# # # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# # # plt.show()
# # # quit()
# #
# # params = [
# # 2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # 1.223347573677009020e-07,
# # 1.847763641629712739e-02,
# # 3.000388537507257813e-02,
# # 2.001077439052020398e-05,
# # 1.219295486333531707e-02,
# # 4.943717334826855457e-05,
# # 2.624725680910927617e-02,
# # 1.991351857905033468e-05,
# # 4.488812697741692559e-03,
# # 1.041577480241394007e-05,
# # 3.473936580568400597e-03,
# # 1.072217409826840662e+00,
# # 1.506747791022632791e-05,
# # 4.360182640163040406e-03,
# # 1.792572991490754334e-05,
# # 2.642504878427706302e-02,
# # 3.558771687312299870e-01,
# # 1.106213891407066866e-06,
# # 2.101813258955530922e-02,
# # 4.020551464879950743e+00,
# # 2.458525808081724495e+00,
# # 3.153989353897445407e-06,
# # 1.676074001871189672e-04,
# # 1.631660435291963120e-03,
# # 3.264925281087857373e-05,
# # 1.555466530711384077e-01,
# # 4.656082056292588645e-01,
# # 7.442659427818243412e-01,
# # 1.663505513158469062e+00,
# # 2.393601347545811427e-05,
# # 5.997125487854754709e-03,
# # 5.270384065069513291e-01,
# # 7.442196406911848194e-02,
# # 4.354575919268037498e-02,
# # 2.466300351383134983e-04,
# # 2.046932131600806510e-05,
# # 2.107029175494326978e+01]
# #
# # new = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # 5.090352718964910340e-05,
# # 8.846681140663302523e-02,
# # 7.688925127286951045e-02,
# # 4.718531991263202608e-05,
# # 2.480645592309011632e-02,
# # 2.634730892196591635e-05,
# # 7.994366563229754474e-02,
# # 8.004684705208197151e-05,
# # 2.787447426045872381e-02,
# # 2.528487936602940290e-05,
# # 5.965677558019213955e-02,
# # 2.605672790260717964e+00,
# # 8.794418149861516561e-05,
# # 7.695471791403299400e-02,
# # 7.898521723337450790e-05,
# # 4.694350552814404581e-02,
# # 6.550932796050894069e+00,
# # 8.186365727139884572e-05,
# # 6.872144915868418080e-02,
# # 1.112414613253134510e+00,
# # 1.497407583427553535e+00,
# # 2.679386619274433570e-06,
# # 1.252840334405051683e-04,
# # 2.142769060164864958e-01,
# # 7.132599412461731482e-05,
# # 1.733224145970182262e+00,
# # 2.345217685036323108e+00,
# # 6.735732655693388304e-01,
# # 5.176940275194190200e+00,
# # 5.364253277853003110e-06,
# # 6.539060406153636429e-02,
# # 7.377928959306133017e+00,
# # 7.512669157064747472e-01,
# # 5.520988157998800439e-02,
# # 8.772775867498763813e-02,
# # 1.123578955588378102e-05,
# # 1.062391316471664356e-02,
# #
# # ]

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

Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))

fstpso = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
3.304257485026848768e-05,
9.791215971368776028e-03,
6.110068937548310437e-03,
4.319219461882335177e-05,
4.212644667502709987e-03,
1.164332269398710173e-05,
2.404257105292715788e-02,
3.311085510054400207e-05,
4.280399320119887552e-02,
2.645814637693955063e-05,
1.437707485722601597e-02,
2.303744002126363599e-01,
2.980688423948379739e-05,
4.879773212139151134e-02,
1.121503480627013705e-05,
1.866712857727401229e-03,
7.572177664736867708e-01,
1.591282730353343167e-05,
3.897146248650381478e-02,
3.076363174012029411e+00,
3.734859661130401243e+00,
3.216200470463125876e-06,
8.782429548444865440e-05,
2.906341314225960662e-02,
5.663104188870970508e-05,
2.110469405515222677e-02,
1.294086380531199176e-01,
3.127598126760888220e-01,
4.298489868360909627e-01,
2.332910188537793332e-06,
7.077504621536276526e-03,
6.294061533948092091e-01,
6.419313355304218094e-02,
8.584653667640911989e-04,
8.160445062706172072e-05,
4.354383618147421691e-06,
4.278903092658225660e+00
]


# pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
#      3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
#      4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
#      4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
#      3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
#      1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
#      3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
#      5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
#      6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
#      1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
#      9.631340201588598493e-01]

# tnf = [2326, 233, 23, 2]
# lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
# tspan = np.linspace(0, 1440, 1440)
tspan = np.linspace(0, 480, 481)
sim = BngSimulator(model, tspan=tspan)
# result = sim.run(method='ode', param_values=params_fst, initials={TNF(brec=None): tnf})
res = sim.run(method='ode', param_values=fstpso)
# result2 = sim.run(method='ode', param_values=params_ode_new, initials={TNF(brec=None): tnf})
# df = result.dataframe
df1 = res.dataframe

# mlkl = [0,   170, 900, 4880,    9940,    10000]
# y = [0, 1, 2, 4, 6, 8]
color = ['b', 'g', 'm', 'r']

plt.figure()
# for n in range(0,4):
    # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
plt.plot(tspan, res.observables['MLKLa_obs'][:])
    # plt.plot(tspan, df1.loc[n]['MLKLa_obs'].iloc[:],'--',c = color[n],lw = 1.5) #fppf
plt.xlabel('Time [minutes]', fontsize=16)
plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
plt.title('Sensitivity of pMLKL to varying TNFa doses')
plt.ylim(ymax = 10000)
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
plt.show()
# quit()
#
#
#
# pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
#      3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
#      4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
#      4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
#      3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
#      1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
#      3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
#      5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
#      6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
#      1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
#      9.631340201588598493e-01]
#
# tnf = [2326, 233, 23, 2]
# lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
# # tspan = np.linspace(0, 1440, 1440)
# tspan = np.linspace(0, 480, 481)
# sim = BngSimulator(model, tspan=tspan)
# # result = sim.run(method='ode', param_values=params_fst, initials={TNF(brec=None): tnf})
# res = sim.run(method='ode', param_values=pf, initials={TNF(brec=None): tnf})
# # result2 = sim.run(method='ode', param_values=params_ode_new, initials={TNF(brec=None): tnf})
# # df = result.dataframe
# df1 = res.dataframe
#
# # mlkl = [0,   170, 900, 4880,    9940,    10000]
# # y = [0, 1, 2, 4, 6, 8]
# color = ['b', 'g', 'm', 'r']
# plt.figure()
# for n in range(0,4):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan, df1.loc[n]['MLKLa_obs'].iloc[:],'--',c = color[n],lw = 1.5) #fppf
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# plt.show()

#STARTING NECRO_NM FILE

# from pysb.core import *
# from pysb.bng import *
# from pysb.integrate import *
# import matplotlib.pyplot as plt
# import numpy as np
# from pysb.util import alias_model_components
# from pysb.simulator.bng import BngSimulator
#
# Model()
#
# model.enable_synth_deg()
#
#
# # # file = open("output_data.txt", 'r')
# # # print(file.read())
# # #
# # # initials_1 = np.array([p.value for p in list(model.initial_conditions)])
# # # print(initials_1)
# # # quit()
# # # initials = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000]
# # #
# # # for i in range(0,101):
# # #     x = file[0]
# # #     print(x)
# # # quit()
# # #     # nvp = initials + file[i]
# # pre_cal = [  1.00000000e-06,   1.00000000e-03,   1.00000000e-03 ,  1.00000000e-06,
# #    1.00000000e-03,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-06,
# #    1.00000000e-03,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
# #    1.00000000e-06,   1.00000000e-03,   1.00000000e-06 ,  1.00000000e-03,
# #    1.00000000e-01,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
# #    1.00000000e-01,   3.11000000e-07,   3.27000000e-06 ,  1.80000000e-02,
# #    3.27000000e-06,   1.80000000e-02,   3.27000000e-02 ,  1.80000000e-02,
# #    1.00000000e-01,   1.00000000e-06,   1.00000000e-03 ,  1.00000000e-01,
# #    1.00000000e-02,   1.00000000e-03,   1.00000000e-03 ,  1.00000000e-06,
# #    1.00000000e+00]
# # params_fst_pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # 3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
# # 4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
# # 4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
# # 3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
# # 1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
# # 3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
# # 5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
# # 6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
# # 1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
# # 9.631340201588598493e-01
# # ]
# # x = range(1,38)
# # params_pf_un = [
# # 7.244683452113300153e-07,
# # 9.517530155749114040e-02,
# # 8.637265039628022611e-02,
# # 3.193901331928525024e-05,
# # 1.258406237885183361e-02,
# # 7.575966125090899991e-05,
# # 8.364351479918887000e-02,
# # 7.436459484254262205e-05,
# # 7.768251964618522187e-02,
# # 9.819940099507889082e-06,
# # 4.836763065283320590e-02,
# # 5.181076854275994403e-01,
# # 9.056610726977585267e-05,
# # 2.338202234638503667e-03,
# # 2.117046267096549337e-05,
# # 8.693240554107009577e-02,
# # 4.110127984888350738e+00,
# # 7.183481154705118965e-05,
# # 7.015524132740647012e-02,
# # 4.116603440959405447e+00,
# # 9.694023692347162324e+00,
# # 2.218206981846405481e-05,
# # 9.837233297311752931e-05,
# # 1.188372739414592205e+00,
# # 1.004157754027084113e-04,
# # 9.913660402796706794e-01, 1.965754591697516096e+00, 1.769000979470284785e+00, 1.129174042737101313e+00,
# # 9.269724481355695611e-05, 9.688283102144903958e-02, 6.520509822693891344e+00, 4.227833110161860475e-01,
# # 3.229777219239124333e-03, 1.567066653485880284e-02, 4.805448621305080279e-05, 1.389352872855412979e+00
# # ]
# #
# # params_fst = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # 3.304257485026848768e-05,
# # 9.791215971368776028e-03,
# # 6.110068937548310437e-03,
# # 4.319219461882335177e-05,
# # 4.212644667502709987e-03,
# # 1.164332269398710173e-05,
# # 2.404257105292715788e-02,
# # 3.311085510054400207e-05,
# # 4.280399320119887552e-02,
# # 2.645814637693955063e-05,
# # 1.437707485722601597e-02,
# # 2.303744002126363599e-01,
# # 2.980688423948379739e-05,
# # 4.879773212139151134e-02,
# # 1.121503480627013705e-05,
# # 1.866712857727401229e-03,
# # 7.572177664736867708e-01,
# # 1.591282730353343167e-05,
# # 3.897146248650381478e-02,
# # 3.076363174012029411e+00,
# # 3.734859661130401243e+00,
# # 3.216200470463125876e-06,
# # 8.782429548444865440e-05,
# # 2.906341314225960662e-02,
# # 5.663104188870970508e-05,
# # 2.110469405515222677e-02,
# # 1.294086380531199176e-01,
# # 3.127598126760888220e-01,
# # 4.298489868360909627e-01,
# # 2.332910188537793332e-06,
# # 7.077504621536276526e-03,
# # 6.294061533948092091e-01,
# # 6.419313355304218094e-02,
# # 8.584653667640911989e-04,
# # 8.160445062706172072e-05,
# # 4.354383618147421691e-06,
# # 4.278903092658225660e+00
# # ]
# #
# # # plt.figure(figsize=(20,10))
# # # # plt.scatter(x, params, label = 'PFu', color = 'blue', lw =1)
# # # plt.scatter(x, params_fst, label = 'FST-PSO', color = 'red', lw =1)
# # # plt.scatter(x, params_fst_pf, label ='FP-PF', color = 'orange', lw = 1)
# # # plt.scatter(x, pre_cal, label ='preCal', color = 'green', lw = 1)
# # # plt.xlabel('kinetic parameters', fontsize=25)
# # # plt.ylabel('parameter value', fontsize=25)
# # # plt.xticks(x)
# # # plt.ylim(ymax = 1)
# # # # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # # plt.legend(loc = 'best',prop={'size': 20})
# # # # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# # # # plt.show()
# # # # quit()
# # #
# # # params = [
# # # 2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # # 1.223347573677009020e-07,
# # # 1.847763641629712739e-02,
# # # 3.000388537507257813e-02,
# # # 2.001077439052020398e-05,
# # # 1.219295486333531707e-02,
# # # 4.943717334826855457e-05,
# # # 2.624725680910927617e-02,
# # # 1.991351857905033468e-05,
# # # 4.488812697741692559e-03,
# # # 1.041577480241394007e-05,
# # # 3.473936580568400597e-03,
# # # 1.072217409826840662e+00,
# # # 1.506747791022632791e-05,
# # # 4.360182640163040406e-03,
# # # 1.792572991490754334e-05,
# # # 2.642504878427706302e-02,
# # # 3.558771687312299870e-01,
# # # 1.106213891407066866e-06,
# # # 2.101813258955530922e-02,
# # # 4.020551464879950743e+00,
# # # 2.458525808081724495e+00,
# # # 3.153989353897445407e-06,
# # # 1.676074001871189672e-04,
# # # 1.631660435291963120e-03,
# # # 3.264925281087857373e-05,
# # # 1.555466530711384077e-01,
# # # 4.656082056292588645e-01,
# # # 7.442659427818243412e-01,
# # # 1.663505513158469062e+00,
# # # 2.393601347545811427e-05,
# # # 5.997125487854754709e-03,
# # # 5.270384065069513291e-01,
# # # 7.442196406911848194e-02,
# # # 4.354575919268037498e-02,
# # # 2.466300351383134983e-04,
# # # 2.046932131600806510e-05,
# # # 2.107029175494326978e+01]
# # #
# # # new = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# # # 5.090352718964910340e-05,
# # # 8.846681140663302523e-02,
# # # 7.688925127286951045e-02,
# # # 4.718531991263202608e-05,
# # # 2.480645592309011632e-02,
# # # 2.634730892196591635e-05,
# # # 7.994366563229754474e-02,
# # # 8.004684705208197151e-05,
# # # 2.787447426045872381e-02,
# # # 2.528487936602940290e-05,
# # # 5.965677558019213955e-02,
# # # 2.605672790260717964e+00,
# # # 8.794418149861516561e-05,
# # # 7.695471791403299400e-02,
# # # 7.898521723337450790e-05,
# # # 4.694350552814404581e-02,
# # # 6.550932796050894069e+00,
# # # 8.186365727139884572e-05,
# # # 6.872144915868418080e-02,
# # # 1.112414613253134510e+00,
# # # 1.497407583427553535e+00,
# # # 2.679386619274433570e-06,
# # # 1.252840334405051683e-04,
# # # 2.142769060164864958e-01,
# # # 7.132599412461731482e-05,
# # # 1.733224145970182262e+00,
# # # 2.345217685036323108e+00,
# # # 6.735732655693388304e-01,
# # # 5.176940275194190200e+00,
# # # 5.364253277853003110e-06,
# # # 6.539060406153636429e-02,
# # # 7.377928959306133017e+00,
# # # 7.512669157064747472e-01,
# # # 5.520988157998800439e-02,
# # # 8.772775867498763813e-02,
# # # 1.123578955588378102e-05,
# # # 1.062391316471664356e-02,
# # #
# # # ]
#
# Monomer('TNF', ['brec'])
# Monomer('TNFR', ['blig', 'brip', 'bDD'])
# Monomer('TRADD', ['brec', 'brip', 'state','bDD1', 'bDD2'], {'state': ['unmod', 'K63ub']})
# Monomer('RIP1', ['bscf', 'bub1', 'bub2', 'bub3','bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state': ['unmod', 'K63ub', 'deub','ub', 'po4', 'trunc']})
# Monomer('TRAF', ['brip', 'bciap', 'bcyld', 'state'], {'state': ['unmod', 'K63ub']})
# Monomer('cIAP', ['btraf'])
# Monomer('A20', ['brip'])
# Monomer('A20t')
# Monomer('CYLD', ['brip','btraf'])
# Monomer('FADD', ['bDD', 'bDED1', 'bDED2'])
# Monomer('proC8', ['bDED'])
# Monomer('C8', ['bf', 'flip', 'state'], {'state': ['A', 'I']})
# Monomer('flip_L', ['bDED', 'state'], {'state': ['A', 'I']})
# Monomer('RIP3', ['bRHIM', 'bDD', 'state'], {'state': ['unmod', 'po4', 'trunc', 'N']})
# Monomer('MLKL', ['bRHIM', 'state'], {'state': ['unmod', 'active', 'inactive']})
# Monomer('LUBAC', ['brip'])
#
#
# Parameter('TNF_0', 2326)
# Parameter('TNFR_0', 4800)
# Parameter('TRADD_0', 9000)
# Parameter('RIP1_0', 40000)
# Parameter('TRAF_0', 9000)
# Parameter('cIAP_0', 9000)
# Parameter('A20_0', 9000)
# Parameter('CYLD_0', 9000)
# Parameter('FADD_0', 8030)
# Parameter('flip_L_0', 3900)
# Parameter('Lubac_0', 7226)
# Parameter('C8_0', 9000)
# Parameter('RIP3_0', 40000)
# Parameter('MLKLa_0', 10000)
#
#
# Initial(TNF(brec=None), TNF_0)
# Initial(TNFR(blig=None, brip=None, bDD = None), TNFR_0)
# Initial(TRADD(brec=None, brip=None, state='unmod', bDD1 = None, bDD2 = None), TRADD_0)
# Initial(RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None, bDD = None, bRHIM = None, bMLKL = None, state='unmod'), RIP1_0)
# Initial(TRAF(brip=None, bciap=None, bcyld = None, state='unmod'), TRAF_0)
# Initial(cIAP(btraf=None), cIAP_0)
# Initial(A20(brip = None), A20_0)
# Initial(CYLD(brip=None, btraf = None), CYLD_0)
# Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)
# Initial(RIP3(bRHIM=None, bDD = None, state='unmod'), RIP3_0)
# Initial(flip_L(bDED=None, state = 'A'), flip_L_0)
# Initial(LUBAC(brip=None), Lubac_0)
# Initial(C8(bf=None, flip = None, state='I'), C8_0)
# Initial(MLKL(bRHIM=None, state='unmod'), MLKLa_0)
#
#
# #COMPLEX I FORMATION AND RELEASE OF RIP1(K63)
# Parameter('bind_TNF_TNFR_kf', 1e-6)
# Parameter('bind_TNF_TNFR_kr',1e-3)
# Parameter('TNF_deg1', 0.001)
# Parameter('bind_TNFRANY_TRADD_kf', 1e-6)
# Parameter('bind_TNFRANY_TRADD_kr',1e-3)
# Parameter('bind_TNFRANY_RIP1unmod_kf', 1e-6)
# Parameter('bind_TNFRANY_RIP1unmod_kr', 1e-3)
# Parameter('bind_RIP1ANY_TRAFunmod_kf', 1e-6)
# Parameter('bind_RIP1ANY_TRAFunmod_kr', 1e-3)
# Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kf', 1e-6)
# Parameter('bind_cIAP_TRAFunmod_to_cIAPTRAFunmod_kr',1e-3)
#
#
#
# Rule('bind_TNF_TNFR', TNF(brec=None) + TNFR(blig=None, brip=None) <> TNF(brec=1) % TNFR(blig=1, brip=None), bind_TNF_TNFR_kf, bind_TNF_TNFR_kr)
#
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
#
#
# Parameter('CompI_UB2', 1e-1)
# Parameter('bind_LUBAC_kf', 1e-6)
# Parameter('bind_LUBAC_kr', 1e-3)
# Parameter('bind_RIP1K63ubANY_A20_kf', 1e-6)
# Parameter('bind_RIP1K63ubANY_A20_kr',1e-3)
# Parameter('k_A20_1', 1e-1)
# Parameter('bind_RIP1K63ubANY_CYLD_kf', 1e-6)
# Parameter('bind_RIP1K63ubANY_CYLD_kr', 1e-3)
# Parameter('k_CYLD_1', 1e-1)
#
# Rule('Complex_I_ubiquitylation2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='unmod') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5)
#      >> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5), CompI_UB2)
#
#
# Rule('ComplexI_Lubac', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) + LUBAC(brip = None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6),
#      bind_LUBAC_kf, bind_LUBAC_kr)
#
#
# #RIP1 K63ub to be deub by A20 or CYLD
#
# Rule('bind_RIP1K63ubANY_A20', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + A20(brip=None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7),
#      bind_RIP1K63ubANY_A20_kf, bind_RIP1K63ubANY_A20_kr)
#
#
# Rule('A20_2', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % A20(brip=7)
#      >>  TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + A20(brip=None), k_A20_1)
#
#
#
# Rule('bind_RIP1K63ubANY_CYLD', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) + CYLD(brip=None)
#      <> TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7),
#      bind_RIP1K63ubANY_CYLD_kf, bind_RIP1K63ubANY_CYLD_kr)
#
# Rule('A20_1', TNF(brec = 1) % TNFR(blig=1, brip=2) % TRADD(brec = 2, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=4, bub1=6, bub2=7, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='K63ub') % TRAF(brip=4, bciap=5, bcyld =None, state='unmod') % cIAP(btraf = 5) % LUBAC(brip = 6) % CYLD(brip=7)
#      >> TNF(brec = None) + TNFR(blig=None, brip=None) + TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + TRAF(brip=None, bciap=None, bcyld = None, state='unmod') + cIAP(btraf = None) + LUBAC(brip = None) + CYLD(brip=None),
#      k_CYLD_1)
#
# #Initiating Necroptosis
# Parameter('bind_TRADDANYRIP1ANY_FADD_kf', 1e-1)
# Parameter('bind_TRADDANYRIP1ANY_FADD_kr', 3.11e-7)
#
# Parameter('bind_FADD_proC8_2_kf', 3.27e-06)
# Parameter('bind_FADD_proC8_2_kr', 0.018)
#
# Parameter('bind_FADDANY_flip_L_kf',3.27e-06)
# Parameter('bind_FADDANY_flip_L_kr', 0.018)
#
# Parameter('bind_C8_flip_L_kf',3.27e-2)
# Parameter('bind_C8_flip_L_kr', 0.018) #not used
#
# Parameter('kc_c8_1', 1e-1)
# Parameter('bind_FADDANY_RIP3_kf', 1e-6)
# Parameter('bind_FADDANY_RIP3_kr', 1e-3)
# Parameter('kc_c8_2', 1e-1)
#
# #RIP1 deub and necrosome formation
#
# Rule('bind_TRADDANYRIP1ANY_FADD', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='deub') + FADD(bDD=None, bDED1 = None, bDED2 = None)
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None), bind_TRADDANYRIP1ANY_FADD_kf, bind_TRADDANYRIP1ANY_FADD_kr)
#
#
# Rule('bind_FADD_proC8_2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + C8(bf=None, flip = None, state='I')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = None, state='I'), bind_FADD_proC8_2_kf, bind_FADD_proC8_2_kr)
#
#
# Rule('bind_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2, flip = None,state='I') + flip_L(bDED=None, state = 'A')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A'), bind_FADDANY_flip_L_kf, bind_FADDANY_flip_L_kr)
#
# Rule('activate_C8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='I') % flip_L(bDED=4, state = 'A')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A'), bind_C8_flip_L_kf,bind_C8_flip_L_kr)
#
# Rule('catalyze_FADDANY_flip_L', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = 2, bDED2 = None) % C8(bf=2,flip = 4, state='A') % flip_L(bDED=4, state = 'A') >>
#      TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + RIP1(bscf=None, btraf=None, bub1=None, bub2=None, bub3=None,bDD = None,bRHIM=None,bMLKL=None, state='trunc') + FADD(bDD=None,bDED1 = None, bDED2 = None) + C8(bf=None,flip = None, state='I') + flip_L(bDED=None, state = 'I') , kc_c8_1)
#
#
# #RIP3 reactions to MLKL
#
# Rule('bind_FADDANY_proC8', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=None,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) + RIP3(bRHIM=None, bDD = None, state='unmod')
#      <> TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod'), bind_FADDANY_RIP3_kf, bind_FADDANY_RIP3_kr)
#
# Rule('C8_activation2', TRADD(brec = None, brip = 3, bDD1=None, bDD2=None) % RIP1(bscf=3, btraf=None, bub1=None, bub2=None, bub3=None,bDD = 1,bRHIM=5,bMLKL=None, state='deub') % FADD(bDD=1,bDED1 = None, bDED2 = None) % RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> TRADD(brec = None, brip = None, bDD1=None, bDD2=None) + FADD(bDD=None,bDED1 = None, bDED2 = None) + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod'), kc_c8_2)
#
# Parameter('bind_RIP1_RIP3po4_kf', 1e-2)
# Parameter('RIP1po4_RIP3po4_kf', 1e-3)
# Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf', 1e-3)
# Parameter('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr', 1e-6)
# Parameter('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc', 1)
#
# Rule('bind_FADDANYANYflip_LANYproC8ANY_RIP1unmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='unmod')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4'), bind_RIP1_RIP3po4_kf)
#
# Rule('Rip1_PO4lation', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'deub')% RIP3(bRHIM=5, bDD = None, state='po4')
#      >> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4'), RIP1po4_RIP3po4_kf)
#
# Rule('bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') + MLKL(bRHIM=None, state='unmod')
#      <> RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod'), bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kf,bind_RIP1po4_MLKLunmod_to_RIP1po4MLKLunmod_kr)
#
# Rule('catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive', RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = 1, state='po4') % MLKL(bRHIM=1, state='unmod')
#      >>  MLKL(bRHIM=None, state='active') + RIP1(bscf = None, bub1 = None, bub2 = None, bub3 = None, bDD=None, btraf=None, bMLKL = None, bRHIM = 5, state = 'po4')% RIP3(bRHIM=5, bDD = None, state='po4') , catalyze_RIP1po4MLKLunmod_to_RIP1po4_MLKLactive_kc)
#
# generate_equations(model)
#
# Observable('MLKLa_obs', MLKL(bRHIM=None, state='active'))
#
# fstpso = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
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
#
#
# # pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# #      3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
# #      4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
# #      4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
# #      3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
# #      1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
# #      3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
# #      5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
# #      6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
# #      1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
# #      9.631340201588598493e-01]
#
# # tnf = [2326, 233, 23, 2]
# # lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
# # tspan = np.linspace(0, 1440, 1440)
# tspan = np.linspace(0, 480, 481)
# sim = BngSimulator(model, tspan=tspan)
# # result = sim.run(method='ode', param_values=params_fst, initials={TNF(brec=None): tnf})
# res = sim.run(method='ode', param_values=fstpso)
# # result2 = sim.run(method='ode', param_values=params_ode_new, initials={TNF(brec=None): tnf})
# # df = result.dataframe
# df1 = res.dataframe
#
# # mlkl = [0,   170, 900, 4880,    9940,    10000]
# # y = [0, 1, 2, 4, 6, 8]
# color = ['b', 'g', 'm', 'r']
#
# plt.figure()
# # for n in range(0,4):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
# plt.plot(tspan, res.observables['MLKLa_obs'][:])
#     # plt.plot(tspan, df1.loc[n]['MLKLa_obs'].iloc[:],'--',c = color[n],lw = 1.5) #fppf
# plt.xlabel('Time [minutes]', fontsize=16)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# plt.title('Sensitivity of pMLKL to varying TNFa doses')
# plt.ylim(ymax = 10000)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# plt.show()
# # quit()
# #
# #
# #
# # pf = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 24000, 10000,
# #      3.497841355176434231e-05, 9.299353889254255087e-03, 1.222381727203981049e-02, 4.662733762385426803e-05,
# #      4.122411164261144054e-03, 6.539070372247510956e-06, 2.659930219941081947e-02, 4.278040527363371738e-05,
# #      4.349323568103394266e-02, 3.003477662783744624e-05, 7.944611372831474755e-03, 7.045871452352303610e-01,
# #      3.598099326725813344e-05, 5.181366206726087387e-02, 3.262502025743792904e-06, 6.553498738869752323e-03,
# #      1.007972408904296735e+00, 2.209541484408219299e-05, 4.529066165519779474e-02, 4.130360397594543542e+00,
# #      3.680177142525775213e+00, 2.283467806267297056e-06, 9.570718102817205031e-05, 3.109002812115645165e-03,
# #      5.730581943892969644e-05, 4.092670496970074456e-02, 1.400509228069881928e-01, 1.134783494309993535e-01,
# #      6.758187648940934267e-01, 1.741332876964058165e-06, 5.478726052303230058e-03, 4.769741203606630564e-01,
# #      1.193094491836479282e-01, 6.923045210935715359e-04, 6.980933857516307892e-03, 3.242773819946248566e-06,
# #      9.631340201588598493e-01]
# #
# # tnf = [2326, 233, 23, 2]
# # lab = ['100 ng/ml', '10 ng/ml', '1 ng/ml', '0.1 ng/ml']
# # # tspan = np.linspace(0, 1440, 1440)
# # tspan = np.linspace(0, 480, 481)
# # sim = BngSimulator(model, tspan=tspan)
# # # result = sim.run(method='ode', param_values=params_fst, initials={TNF(brec=None): tnf})
# # res = sim.run(method='ode', param_values=pf, initials={TNF(brec=None): tnf})
# # # result2 = sim.run(method='ode', param_values=params_ode_new, initials={TNF(brec=None): tnf})
# # # df = result.dataframe
# # df1 = res.dataframe
# #
# # # mlkl = [0,   170, 900, 4880,    9940,    10000]
# # # y = [0, 1, 2, 4, 6, 8]
# # color = ['b', 'g', 'm', 'r']
# # plt.figure()
# # for n in range(0,4):
# #     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
# #     plt.plot(tspan, df1.loc[n]['MLKLa_obs'].iloc[:],'--',c = color[n],lw = 1.5) #fppf
# # plt.xlabel('Time [minutes]', fontsize=16)
# # plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=16)
# # plt.title('Sensitivity of pMLKL to varying TNFa doses')
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# # plt.show()