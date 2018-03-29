import numpy as np
from pydream.core import run_dream
from pysb.integrate import Solver
from pydream.parameters import SampledParam
from scipy.stats import norm, uniform
from pydream.convergence import Gelman_Rubin
from jnk3_no_ask1 import model
import pandas as pd

# Initialize PySB solver

exp_data = pd.read_csv('../data/exp_data_arrestin_noarrestin.csv')

tspan = exp_data['Time (secs)'].values[:-1]
solver = Solver(model, tspan=tspan)

like_mkk4_arrestin_pjnk3 = norm(loc=exp_data['pTyr_arrestin_avg'].values[:-1] + np.finfo(float).eps,
                                scale=exp_data['pTyr_arrestin_std'].values[:-1])
like_mkk7_arrestin_pjnk3 = norm(loc=exp_data['pThr_arrestin_avg'].values[:-1] + np.finfo(float).eps,
                                scale=exp_data['pThr_arrestin_std'].values[:-1])

like_mkk4_noarrestin_pjnk3 = norm(loc=exp_data['pTyr_noarrestin_avg'].values[:-1] + np.finfo(float).eps,
                                scale=exp_data['pTyr_noarrestin_std'].values[:-1])
like_mkk7_noarrestin_pjnk3 = norm(loc=exp_data['pThr_noarrestin_avg'].values[:-1] + np.finfo(float).eps,
                                scale=exp_data['pThr_noarrestin_std'].values[:-1])

pysb_sampled_parameter_names = ['kr_Ask1_Arr', 'kcat_Ask1_Activation', 'kf_uMKK4_to_uMKK7',
                                'kf_uMKK7_to_uMKK4', 'kcat_pMKK4_JNK3', 'kcat_pMKK7_JNK3']

# Add PySB rate parameters to be sampled as unobserved random variables to DREAM with normal priors

idx_pars_calibrate = [3, 21, 23, 25, 27, 29, 32, 33, 34, 35, 36, 37,  39, 41]
rates_of_interest_mask = [i in idx_pars_calibrate for i, par in enumerate(model.parameters)]

# Initial conditions of mkk4, mkk7, uJNK3 respectively
arrestin_idx = [42]
param_values = np.array([p.value for p in model.parameters])

sampled_parameter_names = [SampledParam(norm, loc=np.log10(par), scale=1.5) for par in param_values[rates_of_interest_mask]]

# kd_ASK1_Arr = SampledParam(norm, loc=np.log10(model.parameters['kr_Ask1_Arr'].value /
#                                               model.parameters['kf_Ask1_Arr'].value), scale=1.5)
# kcat_Ask1_Activation = SampledParam(norm, loc=np.log10(model.parameters['kcat_Ask1_Activation'].value), scale=1.5)
# kf_uMKK4_to_uMKK7 = SampledParam(norm, loc=np.log10(model.parameters['kf_uMKK4_to_uMKK7'].value), scale=1.5)
# kf_uMKK7_to_uMKK4 = SampledParam(norm, loc=np.log10(model.parameters['kf_uMKK7_to_uMKK4'].value), scale=1.5)
# kcat_pMKK4_JNK3 = SampledParam(norm, loc=np.log10(model.parameters['kcat_pMKK4_JNK3'].value), scale=1.5)
# kcat_pMKK7_JNK3 = SampledParam(norm, loc=np.log10(model.parameters['kcat_pMKK7_JNK3'].value), scale=1.5)
#
# sampled_parameter_names = [kd_ASK1_Arr, kcat_Ask1_Activation, kf_uMKK4_to_uMKK7, kf_uMKK7_to_uMKK4,
#                            kcat_pMKK4_JNK3, kcat_pMKK7_JNK3]

nchains = 5
niterations = 1000


def likelihood(position):
    Y = np.copy(position)
    param_values[rates_of_interest_mask] = 10 ** Y

    pars = np.copy(param_values)

    solver.run(param_values=pars)
    logp_mkk4_arrestin = np.sum(like_mkk4_arrestin_pjnk3.logpdf(solver.yobs['pTyr_jnk3']))
    logp_mkk7_arrestin = np.sum(like_mkk7_arrestin_pjnk3.logpdf(solver.yobs['pThr_jnk3']))

    # e_mkk4 = np.sum((exp_data['pTyr_arrestin_avg'].values[:-1] - sim['pTyr_jnk3']) ** 2 /
    #                 (2 * exp_data['pTyr_arrestin_std'].values[:-1])) / len(exp_data['pTyr_arrestin_std'].values[:-1])
    # e_mkk7 = np.sum((exp_data['pThr_arrestin_avg'].values[:-1] - sim['pThr_jnk3']) ** 2 /
    #                 (2 * exp_data['pThr_arrestin_std'].values[:-1])) / len(exp_data['pThr_arrestin_std'].values[:-1])
    # error1 = e_mkk4 + e_mkk7

    # No arrestin experiments
    pars[arrestin_idx] = 0
    solver.run(param_values=pars)
    # e2_mkk4 = np.sum((exp_data['pTyr_noarrestin_avg'].values[:-1] - sim2['pTyr_jnk3']) ** 2 /
    #                 (2 * exp_data['pTyr_noarrestin_std'].values[:-1])) / len(exp_data['pTyr_noarrestin_std'].values[:-1])
    # e2_mkk7 = np.sum((exp_data['pThr_noarrestin_avg'].values[:-1] - sim2['pThr_jnk3']) ** 2 /
    #                 (2 * exp_data['pThr_noarrestin_std'].values[:-1])) / len(exp_data['pThr_noarrestin_std'].values[:-1])
    # error2 = e2_mkk4 + e2_mkk7
    # error = error1 + error2


    logp_mkk4_noarrestin = np.sum(like_mkk4_arrestin_pjnk3.logpdf(solver.yobs['pTyr_jnk3']))
    logp_mkk7_noarrestin = np.sum(like_mkk7_arrestin_pjnk3.logpdf(solver.yobs['pThr_jnk3']))

    return logp_mkk4_arrestin + logp_mkk7_arrestin + logp_mkk4_noarrestin + logp_mkk7_noarrestin


if __name__ == '__main__':

    # Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    converged = False
    total_iterations = niterations
    sampled_params, log_ps = run_dream(parameters=sampled_parameter_names, likelihood=likelihood,
                                       niterations=niterations, nchains=nchains, multitry=False,
                                       gamma_levels=4, adapt_gamma=True, history_thin=1,
                                       model_name='jnk3_dreamzs_5chain', verbose=True)

    # Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save('jnk3_dreamzs_5chain_sampled_params_chain_' + str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save('jnk3_dreamzs_5chain_logps_chain_' + str(chain)+'_'+str(total_iterations), log_ps[chain])

    #Check convergence and continue sampling if not converged

    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ',total_iterations,' GR = ',GR)
    np.savetxt('jnk3_dreamzs_5chain_GelmanRubin_iteration_'+str(total_iterations)+'.txt', GR)

    old_samples = sampled_params
    if np.any(GR>1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            sampled_params, log_ps = run_dream(parameters=sampled_parameter_names, likelihood=likelihood,
                                               niterations=niterations, nchains=nchains, start=starts, multitry=False, gamma_levels=4,
                                               adapt_gamma=True, history_thin=1, model_name='jnk3_dreamzs_5chain',
                                               verbose=True, restart=True)


            # Save sampling output (sampled parameter values and their corresponding logps).
            for chain in range(len(sampled_params)):
                np.save('jnk3_dreamzs_5chain_sampled_params_chain_' + str(chain)+'_'+str(total_iterations), sampled_params[chain])
                np.save('jnk3_dreamzs_5chain_logps_chain_' + str(chain)+'_'+str(total_iterations), log_ps[chain])

            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ',total_iterations,' GR = ',GR)
            np.savetxt('jnk3_dreamzs_5chain_GelmanRubin_iteration_' + str(total_iterations)+'.txt', GR)

            if np.all(GR<1.2):
                converged = True

    try:
        #Plot output
        import seaborn as sns
        from matplotlib import pyplot as plt
        total_iterations = len(old_samples[0])
        burnin = total_iterations/2
        samples = np.concatenate((old_samples[0][burnin:, :], old_samples[1][burnin:, :], old_samples[2][burnin:, :],
                                  old_samples[3][burnin:, :], old_samples[4][burnin:, :]))

        ndims = len(sampled_parameter_names)
        colors = sns.color_palette(n_colors=ndims)
        for dim in range(ndims):
            fig = plt.figure()
            sns.distplot(samples[:, dim], color=colors[dim], norm_hist=True)
            fig.savefig('PyDREAM_jnk3_dimension_'+str(dim))

    except ImportError:
        pass

else:

    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':likelihood, 'niterations':niterations, 'nchains':nchains, \
                  'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'jnk3_dreamzs_5chain', 'verbose':False}
