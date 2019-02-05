import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import StochKitSimulator
from necroptosismodule import model
TNF = model.monomers['TNF']

import datetime
now = datetime.datetime.now()
fstpso = [2326, 4800, 9000, 40000, 9000, 9000, 9000, 9000, 8030, 3900, 7226, 9000, 40000, 10000,
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

# Definitions
NUM_SSA_RUNS = 200 #How many times SSA will be ran

#Length of simulation
tspan = np.linspace(0, 3600, 3601)

#turn off graphs showing up
plt.ioff()

#Set Path to save files needed
path = r'\home\ildefog\ParticleSwarmOptimization' #path for diablo


#RUN THROUGH EACH AMOUNT OF TNF:
TNF_LOOP = [('0.1 ng/ml TNF', 2)] #, ('10 ng/ml TNF', 232), ('1 ng/ml TNF', 23), ('.1 ng/ml TNF', 2)]
#TNF_LOOP = [('100 ng/ml TNF', 2326), ('10 ng/ml TNF', 232), ('1 ng/ml TNF', 23), ('.1 ng/ml TNF', 2)]
for tnf_title, dose in TNF_LOOP:

    #RUN STOCHASTIC SIMULATION ALGORITHM (SSA)
    ssa_sim = BngSimulator(model, tspan=tspan, verbose=True)
    ssa_sim_res = ssa_sim.run(method = 'ssa', initials={TNF(brec=None): dose}, n_runs=NUM_SSA_RUNS)
   # ssa_sim_res.save(path + 'SSA_data_%dTNF' % dose)
    df = ssa_sim_res.dataframe #pandas dataframe organizes data

    #FOR EACH OBSERVABLE AVERAGE THE SSA RUNS AT EACH TIME POINT
    avg = df.groupby(level='time').mean()

    #RUN ODE SIMULATION
    ode_sim = ScipyOdeSimulator(model, tspan=tspan)
    ode_sim_res = ode_sim.run(initials={TNF(brec=None): dose})

    # PLOT STOCHASTIC SIMULATION ALGORITHM (SSA) WITH AVG SSA (YELLOW) AND ODE (BLACK)
    # Array: [(Observable name, number to start y axis at, number to end y axis at)]
    # obs_y_range = [('CI_k63_obs', 0, 520), ('CII_obs', 0, 410), ('MLKLa_obs', 0, 10300)]
    #
    # for obs, y1, y2 in obs_y_range:
    #     plt.figure()
    #     plt.title(tnf_title)
    #     plt.ylim(y1, y2)
    #     for _, run in df.groupby('simulation'):
    #             plt.plot(tspan / 60, run.loc[:, obs])
    #     plt.plot(tspan / 60, avg.loc[:, obs], 'gold', linewidth=3)
    #     plt.plot(tspan / 60, ode_sim_res.observables[obs], 'black', linewidth=3, linestyle='dashed')
    #     plt.xlabel("Time (in hr)", fontsize=15)
    #     plt.ylabel("Molecules/Cell", fontsize=15)
    #     plt.title('%s Trajectories' % obs, fontsize=18)
    #     name = 'uncal_run5_%d_SSA_%s_' % (dose, obs)
    #     ssa_name = name + str(now.strftime('%Y-%m-%d_%H%M.png'))
    #     plt.savefig(ssa_name, bbox_inches='tight')

    obs_y_range = [('CI_k63_obs'), ('CII_obs'), ('MLKLa_obs')]

    for obs in obs_y_range:
        plt.figure()
        plt.title(tnf_title)
        # plt.ylim(y1, y2)
        for _, run in df.groupby('simulation'):
                plt.plot(tspan / 60, run.loc[:, obs])
        # plt.plot(tspan / 60, avg.loc[:, obs], 'gold', label = 'ssa_avg',linewidth=3)
        # plt.plot(tspan / 60, ode_sim_res.observables[obs], 'black', linewidth=3,label = 'ode_avg', linestyle='dashed')
        plt.xlabel("Time (in hr)", fontsize=15)
        plt.ylabel("Molecules/Cell", fontsize=15)
        plt.title('%s Trajectories' % obs, fontsize=18)
        # plt.legend(framealpha=1, frameon=True,loc = 'best')
        name = 'uncal_run5_%d_SSA_%s_%d' % (dose, obs,NUM_SSA_RUNS)
        ssa_name = name + str(now.strftime('%Y-%m-%d_%H%M.png'))
        plt.savefig(ssa_name, bbox_inches='tight')
