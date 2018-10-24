from shared_code import *
from pylab import *
import seaborn as sns

import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def simulate(B, model, time_exp):
	fixed_parameters = map(lambda y: y.value, model.parameters[:SPECIES_PARAMETERS])
	new_parameters = list(fixed_parameters) + list(B)
	sim = BngSimulator(model, tspan=time_exp)
	#sim = ScipyOdeSimulator(model, tspan=time_exp)
	x = sim.run(param_values=new_parameters, method="ode")
	return x.all['MLKLobs']


if __name__ == '__main__':
	
	A = genfromtxt("mlkl_data_new.txt", delimiter="\t")
	TSM = A[1:,1:]
	TIME = A[0][1:]
	time_exp = TIME
	time_exp = append([0], time_exp)

	fig, ax = subplots(1,1, figsize=(10,6))

	model = create_model_and_stuff(time_exp)

	APRIORI = loadtxt("2018-08-27_16-47-39_solution") # use best solution from FST-PSO as mean
	# APRIORI = loadtxt("c:\\users\\aresi\\Desktop\\test.txt") # use best solution from FST-PSO as mean
	#2018-08-26_13-29-00_solution
	
	full = linspace(0, TIME[-1], 100)
	res = simulate(APRIORI, model, full)

	plot(full, res, "--", color="red", label="Simulation using $\\textbf{g}^*$")
	plot(TIME, TSM[0], "s", markersize=7, color="gray", label="Experimental data")
	xlabel("Time [minutes]")
	ylabel("Phoshorylated MLKL amount [molecules]")
	legend()
	tight_layout()
	sns.despine()
	savefig("Images/fitting-deterministic-newdata2-noprior.pdf")
	#show()