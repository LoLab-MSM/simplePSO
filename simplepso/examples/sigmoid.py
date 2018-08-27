import numpy as np
import pylab
from scipy.optimize import *
# import curve_fit

xdata = [0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.]
ydata = [0., 0., 0., 0.10, 0.25, 0.5, 0.75, 0.90, 0.98, 1., 1., 1.,1.]

def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a*(x-b)))

popt, pcov = curve_fit(fsigmoid, xdata, ydata, method='dogbox',bounds=((0., 0.),(1., 1.)))

x = np.linspace(-1, 13, 1)
y = fsigmoid(x, *popt)

pylab.plot(xdata, ydata, 'o', label='data')
pylab.plot(x,y, label='fit')
pylab.ylim(0, 1.05)
pylab.legend(loc='best')
pylab.show()