import matplotlib.pyplot as plt
import numpy as np
from pysb.examples.robertson import model
from pysb.integrate import odesolve
import pysb




t = np.linspace(0, 40,100)
obs_names = ['A_total', 'C_total']



solver = pysb.integrate.Solver(model, t, integrator='vode',rtol=1e-8, atol=1e-8)

solver.run()
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)

def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return np.vstack([recarray[name] for name in names]).T