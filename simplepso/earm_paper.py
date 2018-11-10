from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, MatchOnce, Annotation, ANY, WILD
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.util import alias_model_components
# from pysb.simulator.bng import BngSimulator

from pysb.simulator import ScipyOdeSimulator
# # params = [3000, 200, 100, 20000, 1000, 100000, 10000, 10000, 100000, 100000, 100000, 400000, 80000, 20000, 20000, 20000, 20000, 1000, 1000, 500000, 100000,
# # 3.598272649105753460e-06,
# # 2.948336214552251325e-02,
# # 5.243534941752143525e-06,
# # 7.363983998086260041e-07,
# # 7.169469721897811694e-05,
# # 3.065132538283446895e+00,
# # 8.716018090150592549e-06,
# # 4.649002393273928696e-02,
# # 2.873911246633154448e+01,
# # 3.295725878589613063e-05,
# # 1.359833144232308350e-03,
# # 3.699063453666581562e-05,
# # 3.800440699862902422e-03,
# # 2.622798799764678584e-02,
# # 2.830017840491271808e-01,
# # 1.837070394946599405e-02,
# # 1.778678619378387324e-02,
# # 8.068624234945922935e-06,
# # 1.636761396437235241e-02,
# # 2.016988085678836029e+01,
# # 6.942367850910414163e-07,
# # 1.054673463099062465e-03,
# # 8.853478387646071016e-09,
# # 1.800061947954972180e-02,
# # 8.218948383341494690e+00,
# # 1.894186085817541117e-05,
# # 4.245519200540074506e-02,
# # 1.344849305388241183e-04,
# # 5.301585056927489350e-04,
# # 4.730298950676735757e-07,
# # 2.302977783479264254e-02,
# # 1.833489135265260472e+00,
# # 4.038795961548112399e-05,
# # 3.009693815129029296e-03,
# # 1.367035499281706201e+00,
# # 1.343177919422906439e-05,
# # 8.001470801482343798e-02,
# # 1.837554044650650553e+01,
# # 2.890022537525762858e-06,
# # 7.791619197211129000e-03,
# # 2.042148353255123538e+01,
# # 3.725483653994006605e-08,
# # 2.680122584730700551e-03,
# # 6.997865372712475995e+00,
# # 5.858427447096542506e-01,
# # 1.446575502899243373e-02,
# # 1.403090241634387894e-07,
# # 4.111517438351722405e-03,
# # 2.951467435527249883e+00,
# # 1.425181025720966512e-06,
# # 6.686497819351074415e-03,
# # 2.269569788789402853e+00,
# # 4.284468271102515350e-08,
# # 8.617860196955739668e-03,
# # 1.137317294285738001e+01,
# # 6.792541998334622976e-08,
# # 2.429767496786963341e-02,
# # 4.546475422887918505e+00,
# # 4.461065758387060633e-07,
# # 3.933554122474381548e-02,
# # 3.017995061308956650e+01,
# # 5.198827713242158397e-06,
# # 1.692070209359504140e+00,
# # 2.935296145602478668e-05,
# # 5.775808957406885873e-03,
# # 1.241788359455894236e-06,
# # 3.510544357443208252e-01,
# # 1.701797289306908843e-05,
# # 1.533987163929844133e-01,
# # 2.529536518681317346e-05,
# # 4.279431074972147031e-01,
# # 1.025233714483479346e-05,
# # 2.123940377230496279e-01,
# # 3.841666882157645465e-06,
# # 6.010138621017479133e-04,
# # 7.955671510609168420e-06,
# # 4.194255415744318344e-03,
# # 4.516016714659239206e-05,
# # 1.429780748159576731e-01,
# # 5.741353589806094878e-06,
# # 7.571125470247778988e-02,
# # 9.777563214735014458e-04,
# # 9.605923385849690025e-04,
# # 7.493063684368412095e-03,
# # 7.102959730573846914e-03,
# # 3.887547068753295469e-03,
# # 9.009656196079472972e-04,
# # 2.558941179061214981e-03,
# # 4.520593402010276290e-02,
# # 1.659658998375355531e-03,
# # 4.914815570752861969e-03,
# # 1.610473564953765588e-04,
# # 1.278326285947190694e-02,
# # 7.365883129905032684e-05,
# # 1.632045181224880612e-02,
# # 2.264045575094220197e+01,
# # 2.939989952418117686e-06,
# # 1.538289472512236996e-05,
# # 1.524137952432992336e+02,
# # 1.065443737938143165e-04,
# # 3.222415859744079403e-02,
# # 8.185624660956458243e+00,
# # 1.078079038229597456e-03,
# # 6.813272208952928302e-03,
# # 7.501340191534067969e+01
# # ]

Model()

Monomer('L', ['bf'])
Monomer('R', ['bf'])
Monomer('DISC', ['bf'])
Monomer('flip', ['bf'])
Monomer('C8', ['bf', 'state'], {'state': ['pro', 'A']})
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


Parameter('L_0', 3000.0)
Parameter('R_0', 200.0)
Parameter('flip_0', 100.0)
Parameter('C8_0', 20000.0)
Parameter('BAR_0', 1000.0)
Parameter('bind_L_R_to_LR_kf', 4e-07)
Parameter('bind_L_R_to_LR_kr', 0.001)
Parameter('convert_LR_to_DISC_kc', 1e-05)
Parameter('bind_DISC_C8pro_to_DISCC8pro_kf', 1e-06)
Parameter('bind_DISC_C8pro_to_DISCC8pro_kr', 0.001)
Parameter('catalyze_DISCC8pro_to_DISC_C8A_kc', 1.0)
Parameter('bind_C8A_BidU_to_C8ABidU_kf', 1e-06)
Parameter('bind_C8A_BidU_to_C8ABidU_kr', 0.001)
Parameter('catalyze_C8ABidU_to_C8A_BidT_kc', 1.0)
Parameter('bind_DISC_flip_kf', 1e-06)
Parameter('bind_DISC_flip_kr', 0.001)
Parameter('bind_BAR_C8A_kf', 1e-06)
Parameter('bind_BAR_C8A_kr', 0.001)
Parameter('Apaf_0', 100000.0)
Parameter('C3_0', 10000.0)
Parameter('C6_0', 10000.0)
Parameter('C9_0', 100000.0)
Parameter('XIAP_0', 100000.0)
Parameter('PARP_0', 1000000.0)
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
Parameter('Bid_0', 40000.0)
Parameter('BclxL_0', 20000.0)
Parameter('Mcl1_0', 20000.0)
Parameter('Bcl2_0', 20000.0)
Parameter('Bad_0', 1000.0)
Parameter('Noxa_0', 1000.0)
Parameter('CytoC_0', 500000.0)
Parameter('Smac_0', 100000.0)
Parameter('Bax_0', 80000.0)
Parameter('Bak_0', 20000.0)
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


Rule('bind_L_R_to_LR', L(bf=None) + R(bf=None) <> L(bf=1) % R(bf=1), bind_L_R_to_LR_kf, bind_L_R_to_LR_kr)
Rule('convert_LR_to_DISC', L(bf=1) % R(bf=1) >> DISC(bf=None), convert_LR_to_DISC_kc)
Rule('bind_DISC_C8pro_to_DISCC8pro', DISC(bf=None) + C8(bf=None, state='pro') <>  DISC(bf=1) % C8(bf=1, state='pro'), bind_DISC_C8pro_to_DISCC8pro_kf, bind_DISC_C8pro_to_DISCC8pro_kr)
Rule('catalyze_DISCC8pro_to_DISC_C8A', DISC(bf=1) % C8(bf=1, state='pro') >> DISC(bf=None) + C8(bf=None, state='A'), catalyze_DISCC8pro_to_DISC_C8A_kc)
Rule('bind_C8A_BidU_to_C8ABidU', C8(bf=None, state='A') + Bid(bf=None, state='U') <> C8(bf=1, state='A') % Bid(bf=1, state='U'), bind_C8A_BidU_to_C8ABidU_kf, bind_C8A_BidU_to_C8ABidU_kr)
Rule('catalyze_C8ABidU_to_C8A_BidT', C8(bf=1, state='A') % Bid(bf=1, state='U') >> C8(bf=None, state='A') + Bid(bf=None, state='T'), catalyze_C8ABidU_to_C8A_BidT_kc)
Rule('bind_DISC_flip', DISC(bf=None) + flip(bf=None) <> DISC(bf=1) % flip(bf=1), bind_DISC_flip_kf, bind_DISC_flip_kr)
Rule('bind_BAR_C8A', BAR(bf=None) + C8(bf=None, state='A') <> BAR(bf=1) % C8(bf=1, state='A'), bind_BAR_C8A_kf, bind_BAR_C8A_kr)
Rule('equilibrate_SmacC_to_SmacA', Smac(bf=None, state='C') <> Smac(bf=None, state='A'), equilibrate_SmacC_to_SmacA_kf, equilibrate_SmacC_to_SmacA_kr)
Rule('equilibrate_CytoCC_to_CytoCA', CytoC(bf=None, state='C') <> CytoC(bf=None, state='A'), equilibrate_CytoCC_to_CytoCA_kf, equilibrate_CytoCC_to_CytoCA_kr)
Rule('bind_CytoCA_ApafI_to_CytoCAApafI', CytoC(bf=None, state='A') + Apaf(bf=None, state='I') <> CytoC(bf=1, state='A') % Apaf(bf=1, state='I'), bind_CytoCA_ApafI_to_CytoCAApafI_kf, bind_CytoCA_ApafI_to_CytoCAApafI_kr)
Rule('catalyze_CytoCAApafI_to_CytoCA_ApafA', CytoC(bf=1, state='A') % Apaf(bf=1, state='I') >> CytoC(bf=None, state='A') + Apaf(bf=None, state='A'), catalyze_CytoCAApafI_to_CytoCA_ApafA_kc)
Rule('convert_ApafA_C9_to_Apop', Apaf(bf=None, state='A') + C9(bf=None) <> Apop(bf=None), convert_ApafA_C9_to_Apop_kf, convert_ApafA_C9_to_Apop_kr)
Rule('bind_Apop_C3pro_to_ApopC3pro', Apop(bf=None) + C3(bf=None, state='pro') <> Apop(bf=1) % C3(bf=1, state='pro'), bind_Apop_C3pro_to_ApopC3pro_kf, bind_Apop_C3pro_to_ApopC3pro_kr)
Rule('catalyze_ApopC3pro_to_Apop_C3A', Apop(bf=1) % C3(bf=1, state='pro') >> Apop(bf=None) + C3(bf=None, state='A'), catalyze_ApopC3pro_to_Apop_C3A_kc)
Rule('bind_Apop_XIAP', Apop(bf=None) + XIAP(bf=None) <> Apop(bf=1) % XIAP(bf=1), bind_Apop_XIAP_kf, bind_Apop_XIAP_kr)
Rule('bind_SmacA_XIAP', Smac(bf=None, state='A') + XIAP(bf=None) <> Smac(bf=1, state='A') % XIAP(bf=1), bind_SmacA_XIAP_kf, bind_SmacA_XIAP_kr)
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

Initial(L(bf=None), L_0)
Initial(R(bf=None), R_0)
Initial(flip(bf=None), flip_0)
Initial(C8(bf=None, state='pro'), C8_0)
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


# Annotation(L, 'http://identifiers.org/uniprot/P50591', 'is')
# Annotation(R, 'http://identifiers.org/uniprot/O14763', 'is')
# Annotation(DISC, 'http://identifiers.org/obo.go/GO:0031264', 'is')
# Annotation(flip, 'http://identifiers.org/uniprot/O15519', 'is')
# Annotation(C8, 'http://identifiers.org/uniprot/Q14790', 'is')
# Annotation(BAR, 'http://identifiers.org/uniprot/Q9NZS9', 'is')
# Annotation(Bid, 'http://identifiers.org/uniprot/P55957', 'is')
# Annotation(Bax, 'http://identifiers.org/uniprot/Q07812', 'is')
# Annotation(Bak, 'http://identifiers.org/uniprot/Q16611', 'is')
# Annotation(Bcl2, 'http://identifiers.org/uniprot/P10415', 'is')
# Annotation(BclxL, 'http://identifiers.org/uniprot/Q07817', 'is')
# Annotation(Mcl1, 'http://identifiers.org/uniprot/Q07820', 'is')
# Annotation(Bad, 'http://identifiers.org/uniprot/Q92934', 'is')
# Annotation(Noxa, 'http://identifiers.org/uniprot/Q13794', 'is')
# Annotation(CytoC, 'http://identifiers.org/uniprot/P99999', 'is')
# Annotation(Smac, 'http://identifiers.org/uniprot/Q9NR28', 'is')
# Annotation(Apaf, 'http://identifiers.org/uniprot/O14727', 'is')
# Annotation(Apop, 'http://identifiers.org/obo.go/GO:0043293', 'is')
# Annotation(C3, 'http://identifiers.org/uniprot/P42574', 'is')
# Annotation(C6, 'http://identifiers.org/uniprot/P55212', 'is')
# Annotation(C9, 'http://identifiers.org/uniprot/P55211', 'is')
# Annotation(PARP, 'http://identifiers.org/uniprot/P09874', 'is')
# Annotation(XIAP, 'http://identifiers.org/uniprot/P98170', 'is')


# print(len(model.parameters))
# print(len(model.initial_conditions))
# print(len(model.parameters_rules()))
# print(model.parameters)
# # print(model.parameters_rules())
# quit()
# print(len(model.parameters))
# print(len(model.parameters_rules()))
# print(len(model.initial_conditions))
# print(len(model.reactions))

Observable('mBid', Bid(state='M'))
Observable('aSmac', Smac(state='A'))
Observable('cPARP', PARP(state='C'))

tspan = np.linspace(0, 20160, 20161)
sim = ScipyOdeSimulator(model, tspan=tspan)
result = sim.run()

plt.figure()
plt.subplot(131)
plt.plot(tspan, result.observables['mBid'][:])
plt.xlabel('Time [seconds]', fontsize=16)
plt.ylabel('mBid amount [molecules]', fontsize=16)

plt.subplot(132)
plt.plot(tspan, result.observables['aSmac'][:])
plt.xlabel('Time [seconds]', fontsize=16)
plt.ylabel('aSmac amount [molecules]', fontsize=16)

plt.subplot(133)
plt.plot(tspan, result.observables['cPARP'][:])
plt.xlabel('Time [seconds]', fontsize=16)
plt.ylabel('cParp amount [molecules]', fontsize=16)

plt.show()