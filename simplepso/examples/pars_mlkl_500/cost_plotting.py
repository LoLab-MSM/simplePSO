#! /usr/bin/env/python
# _*_ coding: utf-8 _*_
import os
from heapq import nsmallest

try:
    import matplotlib

    # matplotlib.use('TkAgg')
    matplotlib.use('Agg')  # this is needed only if you run it in a cluster
    import matplotlib.pyplot as plt

    plot = True  # this is because in the run_example function you have the "if plot" expression
except ImportError:
    plot = False
    pass

import numpy as np

plt.ioff()

"""reads in the cost values and throws away those more than a std above the average. Plots if plot = True"""
def best10(path_cost, path_par, plot):
    cost_list = []
    for file in os.listdir(path_cost):
        if file.endswith(".npy"):
            cost = np.load(os.path.join(path_cost, file))
            # print("ok")
        else:
            print("no file with this extension")

        cost_list.append(cost)
    print(len(cost_list))
    if plot:
        """plot costs"""
        plt.hist(cost_list, bins=10)
        # plt.xlim(1.74, 1.95)
        plt.xlabel('Cost')
        plt.ylabel('Frequency')
        plt.title('Frequency of cost values')
        plt.show()
        plt.savefig('costs_all_mlkl.png', bbox_inches='tight')
        plt.clf()
        # quit()
    cost_array = np.array(cost_list)
    print(cost_array)
    std = cost_array.std()
    print(std)
    avg = cost_array.mean()
    print(avg)
    cost_selected = []
    # print(cost_list)
    for element in cost_list:
        if abs(element - avg) >= std:
            continue
        else:
            cost_selected.append(element)
    # if plot:
    #     plt.clf()
    #     plt.xlim(1.75, 1.85)
    #     plt.hist(cost_selected, bins=30)
    #     plt.xlabel('Cost')
    #     plt.ylabel('Frequency')
    #     plt.title('Frequency of the selected cost values')
    #     plt.show()
    #     plt.savefig("costs_selected" + '.png', bbox_inches='tight')

    """reads in the parameter sets"""
    print(cost_selected)
    # quit()
    par_list = []
    for file in os.listdir(path_par):
        if file.endswith(".npy"):
            pars = np.load(path_par + "/" + file)
            par_list.append(pars)
            # print("ahoy")
        else:
            print("no file with this extension")
    print(cost_list)
    par_cost = dict([(tuple(par_list[i]), cost_list[i]) for i in range(len(cost_list))])
    par_cost_selected = {k: v for k, v in par_cost.items() if v in cost_selected}

    """select the ten best parameter sets"""
    i = 0
    for _ in range(10):
        minimum = min(par_cost_selected, key=par_cost_selected.get)
        if i == 0:
            best10par = []
        i = i + 1
        best10par.append(minimum)
        del par_cost_selected[minimum]
    if plot:
        best10cost = nsmallest(10,cost_selected)
        # plt.xlim(1.74, 1.95)
        plt.hist(best10cost, bins=30)
        plt.xlabel('Cost')
        plt.ylabel('Frequency')
        plt.title('Frequency of the best 10 cost values')
        plt.show()
        plt.savefig("costs_best10n" + '.png', bbox_inches='tight')
        plt.clf()
    return best10par

# path_par_silvia = r"/Users/geenaildefonso/Projects/Necroptosis/optimizer.best"
# path_cost_silvia = r"/Users/geenaildefonso/Projects/Necroptosis/Cost_variance"
plot = True
path_par_geena = r"/Users/geenaildefonso/Projects/ParticleSwarmOptimization/simplepso/examples/pars_mlkl_500"
path_cost_geena = r"/Users/geenaildefonso/Projects/ParticleSwarmOptimization/simplepso/examples/cost_var_mlkl_500"


if '__main__' == __name__:
    """SILVIA"""
    # best10par = best10(path_cost_silvia, path_par_silvia, plot)
    # for i in range(len(best10par)):
    #     print( best10par[i])
    #     np.save('par_best_silvia_%s' % i, best10par[i])
    """GEENA"""
    best10par = best10(path_cost_geena, path_par_geena, plot)
    for i in range(len(best10par)):
        print(best10par[i])
        np.save('par_best_geenan_%s' % i, best10par[i])