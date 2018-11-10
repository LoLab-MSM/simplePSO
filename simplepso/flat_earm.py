from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, Annotation, ANY, WILD, MatchOnce
import numpy as np
# from pysb.simulator import ScipyOdeSimulator
from matplotlib.pyplot import savefig
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms


# import csv
# x, y = np.loadtxt('smac.csv', delimiter=',', unpack=True)
# plt.plot(x,y, label='Loaded from file!')
#
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Interesting Graph\nCheck it out')
# plt.legend()
# plt.show()
# quit()
# trans = mtransforms.blended_transform_factory(plt.transData, plt.transAxes)
from pysb.simulator.bng import BngSimulator
#
#
# def graph(formula, x_range):
#     x = np.array(x_range)
#     y = eval(formula)
#     plt.plot(x, y)
#     plt.show()
#
# graph('555*x', range(9810,9991,1))
#
# quit()

Model()

#Bid, bclxl, mcl1, bcl2, bad, noxa, cytoc, smac, bax, bak


# params_pf1 = [3000, 200, 100, 20000, 1000,
# 3.478091829880603936e-06, 2.848050023928675248e-02, 1.757031688876240377e-06, 7.891289845099509911e-07, 1.295624125685610351e-03,
# 4.071100091324983516e+00, 9.915286195393510930e-06, 4.666526196433289780e-02, 3.056785543623280077e+01, 3.285259655206213973e-05,
# 8.557276863348971469e-04, 3.725135839713434390e-05,  2.439080571244632914e-03, 100000, 10000, 10000, 100000, 100000, 1000000, 3.234183618443283587e-02, 2.780659771654789836e-01,
# 1.048097032643262316e-02, 2.674590322398063061e-02, 8.516572279160523586e-06, 1.737499372132231268e-02, 2.008375868368192840e+01,
# 6.935338096755502175e-07, 1.136499577852876025e-03, 1.193140532701628915e-08, 1.704916204628463550e-02, 8.742959394189499633e+00,
# 1.789906265823865676e-05, 4.105817161361738354e-02, 1.257025706723589745e-04, 4.584165277907347380e-04, 4.477444975401404352e-07,
# 2.339105449471683834e-02, 1.437047733704246433e+00, 3.904440867402459445e-05, 4.050606833131481337e-03, 1.457357024246903165e+00,
# 1.481127446411864950e-05, 7.931692096744334675e-02, 1.861627590626201112e+01, 4.159901453036069896e-06, 6.949718471014031292e-03,
# 2.168274357751361237e+01, 8.172850066101791072e-09, 2.114022553946774425e-03, 7.991962128089038941e+00,40000, 20000, 20000, 20000, 1000, 1000, 500000, 100000, 80000, 20000, 4.914667053129129770e-01,
# 1.488544919161798909e-02, 8.753418053135645940e-08, 3.553658885386833631e-03, 3.487054552750191938e+00, 1.410938017097035102e-06,
# 5.281833265220347087e-03, 1.838783995087956979e+00, 1.846462486144627379e-07, 7.181094070538411762e-03, 1.414912638138554613e+01,
# 7.918057623442218813e-08, 2.557938516386364103e-02, 4.140187150409025740e+00, 4.628013429646770746e-07, 3.971588715668186026e-02,
# 2.926220419742611867e+01, 4.252592366011532563e-06, 1.772014030804755613e+00, 2.925424231886606470e-05, 1.416385300463742294e-02,
# 2.743581529884081422e-07, 3.489965190145944418e-01, 1.680355524575344650e-05, 1.599775060612102506e-01, 2.758087801413124917e-05,
# 4.058891325205095102e-01, 8.260941261351661434e-06, 2.683939582833357318e-01, 4.421782231009156618e-06, 1.262196020574339042e-03,
# 7.524067995853102508e-06, 2.291182184871169356e-02, 4.353160388915986318e-05, 1.462061024589780267e-01, 8.398656517454044677e-06,
# 7.137140016734222492e-02, 9.489017023052876429e-04, 1.088051487629779473e-03, 7.317755992150802927e-03, 7.698521947424297705e-03,
# 3.636094217774524413e-03, 1.758735980323047985e-03, 2.771032714154535951e-03, 4.539781170042395120e-02, 1.471681236829634103e-03,
# 5.070857154442596527e-03, 1.168458118827437025e-04, 1.316367714276471187e-02, 7.466281118626521109e-05, 1.530637869157581293e-02,
# 1.792821754011748325e+01, 1.730684905406810431e-05, 2.403550248333405283e-04, 1.482867898155281807e+02, 1.830894482001253916e-04,
# 2.948362685351018156e-02, 1.289615954383603125e+01, 1.041500541329146547e-03, 6.126570720621737215e-03, 7.452755852831029415e+01
# ]

params_pf2 = [3000, 200, 100, 20000, 1000,
3.475694840532985580e-07, 3.548300852314466985e-02,2.993769354561528080e-05, 6.412609686981107681e-07,1.392429565586885375e-02,
1.188298631861059462e+01, 4.738004703994751267e-06,3.156199384538281583e-02, 2.741608902541684856e+01,3.287653509701260636e-05,
8.525977033942941027e-03, 3.466167069582748830e-05,3.120582321126550915e-05, 100000, 10000, 10000, 100000, 100000, 1000000,1.174253546462952630e-01,2.708012361695244508e-01,
1.268684099816612121e-01, 3.499021788426404100e-02,5.555670934978483816e-06, 2.942984692366612531e-02,2.392387912717449083e+01,
7.904885289574353146e-07, 7.333357127315982826e-03,2.636420074522533219e-08, 3.400040026174559055e-02,9.607307995627268227e+00,
3.073901208137210162e-05, 4.141228287946343428e-02,6.745923259556057497e-05, 1.063760597252086783e-04,4.312305721417325816e-07,
2.358218248790772478e-02, 3.175693413760052319e+00,4.873933390535256123e-05, 6.817971212285222314e-03,2.919132477195663267e+00,
1.268011373790858949e-05, 6.162573475851259447e-02,1.132569248208334045e+01, 1.279771417692554619e-05,6.538135947297156990e-03,
1.727851846936411917e+01, 1.053109568608831023e-08,8.861296709458745284e-04, 5.272472603765057109e+00,40000, 20000, 20000, 20000, 1000, 1000, 500000, 100000, 80000, 20000, 5.154641867216540607e-01,
6.441896996272443573e-04, 7.805008361443064053e-07,8.226857084331193684e-03, 6.026059192935764308e+00,1.899689602431782005e-06,
1.244225640348547469e-02, 8.037299324475483786e+00,6.635719118797360977e-07, 4.457206285807203440e-03,8.788509108434418238e+00,
8.489212323005825352e-07, 2.043063687052431091e-02,3.666843511206120176e+00, 4.423715634764409734e-08,4.670495589159086997e-02,
2.550152376250139241e+01, 1.318968790321588663e-05,1.701421499088885403e+00, 3.288929887385575482e-05,5.433309627888781648e-03,
7.206856981797517197e-06, 4.946822154967753793e-01,2.527524080008720759e-05, 1.319918901397579558e-01,2.759227569556548396e-05,
4.129378295972826463e-01, 1.586452552366658512e-05,6.707547008617840145e-01, 1.084802318669966126e-05,1.631362913727264774e-02,
4.775612143893501351e-06, 4.343257101518788882e-02,5.004563851242942521e-05, 1.160257281102752697e-01,1.025320434547356969e-05,
4.628442928798731648e-02, 1.299625861886722147e-03,3.991960067140984739e-03, 6.122810453011656197e-03,3.751155021679380454e-03,
3.717296409709610564e-03, 1.142863397251036141e-03,5.939539109287267360e-03, 5.234992009976399685e-02,3.619694981407414009e-04,
4.363312780119012191e-03, 7.719795717353994510e-05,1.700016590401536021e-02, 1.197967768941826422e-04,1.754148028673731263e-02,
1.876093532353690208e+01, 2.119704671964884918e-04,4.684862733038934815e-03, 1.756463003734903623e+02,5.035472991495496637e-04,
3.274130296611463958e-02, 2.879450876104883861e+01,1.010446810371324212e-03, 5.315125221211113200e-03,4.326872537562300636e+01
]

# paramsfst = [3000, 200, 100, 20000, 1000,
# 3.598272649105753460e-06,
# 2.948336214552251325e-02,
# 5.243534941752143525e-06,
# 7.363983998086260041e-07,
# 7.169469721897811694e-05,
# 3.065132538283446895e+00,
# 8.716018090150592549e-06,
# 4.649002393273928696e-02,
# 2.873911246633154448e+01,
# 3.295725878589613063e-05,
# 1.359833144232308350e-03,
# 3.699063453666581562e-05, #14th
# 3.800440699862902422e-03, 100000, 10000, 10000, 100000, 100000, 1000000,
# 2.622798799764678584e-02,
# 2.830017840491271808e-01,
# 1.837070394946599405e-02,
# 1.778678619378387324e-02,
# 8.068624234945922935e-06,
# 1.636761396437235241e-02,
# 2.016988085678836029e+01,
# 6.942367850910414163e-07,
# 1.054673463099062465e-03,
# 8.853478387646071016e-09,
# 1.800061947954972180e-02,
# 8.218948383341494690e+00,
# 1.894186085817541117e-05,
# 4.245519200540074506e-02,
# 1.344849305388241183e-04,
# 5.301585056927489350e-04,
# 4.730298950676735757e-07,
# 2.302977783479264254e-02,
# 1.833489135265260472e+00,
# 4.038795961548112399e-05,
# 3.009693815129029296e-03,
# 1.367035499281706201e+00,
# 1.343177919422906439e-05,
# 8.001470801482343798e-02,
# 1.837554044650650553e+01,
# 2.890022537525762858e-06,
# 7.791619197211129000e-03,
# 2.042148353255123538e+01,
# 3.725483653994006605e-08,
# 2.680122584730700551e-03, #32nd8
# 6.997865372712475995e+00,40000, 20000, 20000, 20000, 1000, 1000, 500000, 100000, 80000, 20000,
# 5.858427447096542506e-01,
# 1.446575502899243373e-02,
# 1.403090241634387894e-07,
# 4.111517438351722405e-03,
# 2.951467435527249883e+00,
# 1.425181025720966512e-06,
# 6.686497819351074415e-03,
# 2.269569788789402853e+00,
# 4.284468271102515350e-08,
# 8.617860196955739668e-03,
# 1.137317294285738001e+01,
# 6.792541998334622976e-08,
# 2.429767496786963341e-02,
# 4.546475422887918505e+00,
# 4.461065758387060633e-07,
# 3.933554122474381548e-02,
# 3.017995061308956650e+01,
# 5.198827713242158397e-06,
# 1.692070209359504140e+00,
# 2.935296145602478668e-05,
# 5.775808957406885873e-03,
# 1.241788359455894236e-06,
# 3.510544357443208252e-01,
# 1.701797289306908843e-05,
# 1.533987163929844133e-01,
# 2.529536518681317346e-05,
# 4.279431074972147031e-01,
# 1.025233714483479346e-05,
# 2.123940377230496279e-01,
# 3.841666882157645465e-06,
# 6.010138621017479133e-04,
# 7.955671510609168420e-06,
# 4.194255415744318344e-03,
# 4.516016714659239206e-05,
# 1.429780748159576731e-01,
# 5.741353589806094878e-06,
# 7.571125470247778988e-02,
# 9.777563214735014458e-04,
# 9.605923385849690025e-04,
# 7.493063684368412095e-03,
# 7.102959730573846914e-03,
# 3.887547068753295469e-03,
# 9.009656196079472972e-04,
# 2.558941179061214981e-03,
# 4.520593402010276290e-02,
# 1.659658998375355531e-03,
# 4.914815570752861969e-03,
# 1.610473564953765588e-04,
# 1.278326285947190694e-02,
# 7.365883129905032684e-05,
# 1.632045181224880612e-02,
# 2.264045575094220197e+01,
# 2.939989952418117686e-06,
# 1.538289472512236996e-05,
# 1.524137952432992336e+02,
# 1.065443737938143165e-04,
# 3.222415859744079403e-02,
# 8.185624660956458243e+00,
# 1.078079038229597456e-03,
# 6.813272208952928302e-03,
# 7.501340191534067969e+01
# ]

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

# Observable('mBid', Bid(state='M'))
# Observable('aSmac', Smac(state='A'))
# Observable('cPARP', PARP(state='C'))
# Observable('C8a', C8(state = 'A'))
# Observable('tBid', Bid(bf=None, state='T'))
# Observable('mBax', Bax(bf=None, state='M'))
# Observable('fBcl2',Bcl2(bf=None) )
# Observable('smac:xiap', Smac(bf=1, state='A') % XIAP(bf=1))
# Observable('C3a', C3(bf=None, state='A'))
# Observable('C3a:xiap',XIAP(bf=1) % C3(bf=1, state='A'))


Rule('bind_L_R_to_LR', L(bf=None) + R(bf=None) <> L(bf=1) % R(bf=1), bind_L_R_to_LR_kf, bind_L_R_to_LR_kr)
Rule('convert_LR_to_DISC', L(bf=1) % R(bf=1) >> DISC(bf=None), convert_LR_to_DISC_kc)
Rule('bind_DISC_C8pro_to_DISCC8pro', DISC(bf=None) + C8(bf=None, state='pro') <> DISC(bf=1) % C8(bf=1, state='pro'), bind_DISC_C8pro_to_DISCC8pro_kf, bind_DISC_C8pro_to_DISCC8pro_kr)
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

Observable('mBid', Bid(state='M'))
Observable('aSmac', Smac(state='A'))
Observable('mSmac', Smac(state='M'))
Observable('cSmac', Smac(state='C'))
Observable('cPARP', PARP(state='C'))
Observable('C8a', C8(state = 'A'))
Observable('tBid', Bid(bf=None, state='T'))
Observable('mBax', Bax(bf=None, state='M'))
Observable('fBcl2',Bcl2(bf=None))
Observable('smac_xiap', Smac(bf=1, state='A') % XIAP(bf=1))
Observable('C3a', C3(bf=None, state='A'))
Observable('C3a_xiap',XIAP(bf=1) % C3(bf=1, state='A'))


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

trail = [740853, 18521, 3000, 740]
tspan = np.linspace(0, 20160, 20161)
sim = BngSimulator(model, tspan=tspan, verbose=True)
# result1 = sim.run(param_values=params_pf1, initials={L(bf=None): trail}
result2 = sim.run(param_values=params_pf2)
a = [9720, 9756, 9797,9828, 9864, 9900]
print(result2.observables['mBid'][a])
# quit()
# df = sim_result.dataframe
#
# # df1 = sim_result.dataframe
# a = (df['aSmac'].iloc[:]).div(100000)
# b = (df['cPARP'].iloc[:]).div(1000000)
# c = (df['tBid'].iloc[:])
# print(c)
# quit()

# print(b)
# y = a.div(100000)
# quit()
# x = [9360/3600, 9540/3600, 9720/3600, 9900/3600, 10080/3600, 10260/3600]
x = [2.6, 2.65, 2.7, 2.75, 2.8, 2.85]


plt.figure()
# # plt.fill_between(x, 0, 100000, facecolor='green', alpha=.3, label = 'MOMP')
plt.plot(tspan, result2.observables['mBid'][:], label = 'aSmac')
# plt.plot(tspan, result2.observables['cSmac'][:], label = 'cSmac')
# plt.plot(tspan, result2.observables['mSmac'][:], label = 'mSmac')
# plt.plot(tspan/3600, a, label = 'MOMP (cytosolic Smac)')
# plt.plot(tspan/3600, b, label = 'EC substrate (cPARP)')
# plt.plot(tspan/3600, c, label = 'IC substrate (tBid)')
plt.xlabel('Time [hours]', fontsize=15)
plt.ylabel('fraction', fontsize=15)
# plt.ylim(ymin = 0, ymax = 1.0)
plt.legend(loc = 'best')
# savefig('icecmomp.pdf')
plt.show()
#
quit()
#
# color = ['b', 'm', 'r', 'g']
# plt.figure()
# for n in range(0,4):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan/3600, df.loc[n]['aSmac'].iloc[:],c = color[n],lw = 4) #fppf
# plt.xlabel('Time [hours]', fontsize=16)
# plt.ylabel('amount [molecules]', fontsize=16)
# plt.title('Sensitivity of Smac to varying TRAIL doses')
# plt.legend(['1000 ng/ml','250 ng/ml ', ' 40 ng/ml', '10ng/ml'] , title = 'TRAIL FPPF', loc=0, fontsize = 8)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# savefig('smac.pdf')
#
# plt.figure()
# for n in range(0,4):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan/3600, df.loc[n]['cPARP'].iloc[:],c = color[n],lw = 4) #fppf
# plt.xlabel('Time [hours]', fontsize=16)
# plt.ylabel('amount [molecules]', fontsize=16)
# plt.title('Sensitivity of cPARP to varying TRAIL doses')
# plt.legend(['1000 ng/ml','250 ng/ml ', ' 40 ng/ml', '10ng/ml'] , title = 'TRAIL FPPF', loc=0, fontsize = 8)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# savefig('cparp.pdf')
#
# plt.figure()
# for n in range(0,4):
#     # plt.plot(tspan, df.loc[n]['MLKLa_obs'].iloc[:], c = color[n],lw =1.5) #fst-pso
#     plt.plot(tspan/3600, df.loc[n]['tBid'].iloc[:],c = color[n],lw = 4) #fppf
# plt.xlabel('Time [hours]', fontsize=16)
# plt.ylabel('amount [molecules]', fontsize=16)
# plt.title('Sensitivity of tBID to varying TRAIL doses')
# plt.legend(['1000 ng/ml','250 ng/ml ', ' 40 ng/ml', '10ng/ml'] , title = 'TRAIL FPPF', loc=0, fontsize = 8)
# # plt.legend(['100 ng/ml ', '10 ng/ml', '1 ng/ml', '0.1 ng/ml'] , title = 'TNF FPPF', loc=0, fontsize = 5)
# # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# savefig('tBid.pdf')
# # plt.show()
# quit()

#
# plt.figure()
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['tBid'].iloc[:]).div(600), color = 'g', label = 'tBid', lw = 1)
# plt.plot(tspan/3600, (df['aSmac'].iloc[:]).div(100000),color = 'm',  label = 'aSmac', lw = 1)
# plt.plot(tspan/3600, (df['cPARP'].iloc[:]).div(1000000), color = 'b', label = 'cPARP', lw = 1)
# plt.plot(tspan/3600, (df['C8a'].iloc[:]).div(10000), color = 'y', label = 'C8a', lw = 1)
# plt.plot(tspan/3600, (df['C3a_xiap'].iloc[:]).div(2.5), color = 'teal', label = 'C3a:XIAP', lw = 1)
# plt.plot(tspan/3600, (df['smac_xiap'].iloc[:]).div(100000), color = 'chocolate', label = 'smac:XIAP', lw = 1)
# plt.xlabel('Time [hours]', fontsize=15)
# plt.ylabel('fraction', fontsize=15)
# plt.ylim(ymin = 0, ymax = 1.0)
# plt.legend(loc = 'best')
# # savefig('fullobs.pdf')
# plt.show()
# quit()

plt.figure(figsize=(20,10))
plt.subplot(3,2,1)
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['tBid'].iloc[:]).div(600), color = 'g', label = 'tBid', lw = 3)
plt.plot(tspan, result2.observables['tBid'][:], color = 'teal', label = 'tBid', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best',prop={'size': 20})

plt.subplot(3,2,2)
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['aSmac'].iloc[:]).div(100000),color = 'm',  label = 'aSmac', lw = 3)
plt.plot(tspan, result2.observables['aSmac'][:],color = 'teal',  label = 'aSmac', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best',prop={'size': 20})
#
plt.subplot(3,2,3)
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['cPARP'].iloc[:]).div(1000000), color = 'b', label = 'cPARP', lw = 3)
plt.plot(tspan, result2.observables['cPARP'][:], color = 'teal', label = 'cPARP', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'upper left',prop={'size': 20})
#
plt.subplot(3,2,4)
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['C8a'].iloc[:]).div(10000), color = 'y', label = 'C8a', lw = 3)
plt.plot(tspan, result2.observables['C8a'][:], color = 'teal', label = 'C8a', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.ylim(ymax = 20000)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best',prop={'size': 20})

plt.subplot(3,2,5)
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['C3a_xiap'].iloc[:]).div(100000), color = 'chocolate', label = 'C3a-xiap', lw = 3)
plt.plot(tspan, result2.observables['C3a_xiap'][:], color = 'teal', label = 'C3a-xiap', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best',prop={'size': 20})
#
plt.subplot(3,2,6)
# plt.fill_between(x, 0, 1, facecolor='red', alpha=.3)
# plt.plot(tspan/3600, (df['smac_xiap'].iloc[:]).div(250),color = 'k',  label = 'smac-xiap', lw = 3)
plt.plot(tspan, result2.observables['smac_xiap'][:],color = 'teal',  label = 'smac-xiap', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'upper left',prop={'size': 20})
# savefig('sepplots.pdf')
plt.show()
quit()
plt.figure(figsize=(10,15))
plt.subplot(5,1,1)
plt.fill_between(x, 0, 500, facecolor='red', alpha=.3, label = 'MOMP')
plt.plot(tspan/3600, result1.observables['tBid'][:],color = 'orange',  label = 'tBidPF', lw = 4)
# plt.plot(tspan/3600, result2.observables['tBid'][:],color = 'k',  label = 'tBidFP', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best')

plt.subplot(5,1,2)
plt.fill_between(x, 0, 6000, facecolor='red', alpha=.3, label = 'MOMP')
plt.plot(tspan/3600, result1.observables['mBax'][:], color = 'brown', label = 'mBaxPF', lw = 4)
# plt.plot(tspan/3600, result2.observables['mBax'][:], color = 'y', label = 'mBaxFP', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best')
#

plt.subplot(5,1,3)
plt.fill_between(x, 15000, 20000, facecolor='red', alpha=.3, label = 'MOMP')
# plt.plot(tspan/3600, result2.observables['smac_xiap'][:], label = 'smac:xiap')
plt.plot(tspan/3600, result1.observables['fBcl2'][:],color = 'purple',  label = 'fBcl2PF', lw = 4)
# plt.plot(tspan/3600, result2.observables['fBcl2'][:],color = 'b',  label = 'fBcl2FP', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best')
#

plt.subplot(5,1,4)
plt.fill_between(x, 0, 100000, facecolor='red', alpha=.3, label = 'MOMP')
plt.plot(tspan/3600, result1.observables['smac_xiap'][:], color = 'chocolate', label = 'smac:XIAPPF', lw = 4)
# plt.plot(tspan/3600, result2.observables['smac_xiap'][:], color = 'm', label = 'smac:XIAPFP', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best')

plt.subplot(5,1,5)
plt.fill_between(x, 0, 2.5, facecolor='red', alpha=.3, label = 'MOMP')
plt.plot(tspan/3600, result1.observables['C3a_xiap'][:], color = 'teal', label = 'C3a:XIAPPF', lw = 4)
# plt.plot(tspan/3600, result2.observables['C3a_xiap'][:], color = 'g', label = 'C3a:XIAPFP', lw = 4)
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('Amount [molecules]', fontsize=16)
plt.tick_params(labelsize=15)
plt.legend(loc = 'best')

# plt.figure(figsize=(20,10))
# plt.subplot(131)
# plt.fill_between(x, 0, 40000, facecolor='green', alpha=.3, label = 'MOMP')
# # plt.plot(tspan, result1.observables['mBid'][:], color = 'r')
# plt.plot(tspan/3600, result2.observables['mBid'][:], label = 'mBid')
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('mBid amount [molecules]', fontsize=20)
# plt.legend(loc = 'best', prop={'size': 20})
#
# plt.subplot(132)
# plt.fill_between(x, 0, 100000, facecolor='green', alpha=.3)
# # plt.plot(tspan, result1.observables['aSmac'][:], color = 'r')
# plt.plot(tspan/3600, result2.observables['aSmac'][:], label = 'aSmac')
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('aSmac amount [molecules]', fontsize=20)
# plt.legend(loc = 'best',prop={'size': 20})
#
# plt.subplot(133)
# plt.fill_between(x, 0, 1000000, facecolor='green', alpha=.3)
# # plt.plot(tspan, result1.observables['cPARP'][:], color = 'r')
# plt.plot(tspan/3600, result2.observables['cPARP'][:], label = 'cPARP')
# plt.xlabel('Time [hrs]', fontsize=20)
# plt.ylabel('cParp amount [molecules]', fontsize=20)
# plt.legend(loc = 'best', prop={'size': 20})
plt.show()