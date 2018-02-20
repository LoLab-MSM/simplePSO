import numpy as np
from scipy import *
from scipy import interpolate
import matplotlib.pyplot as plt

x = np.array([0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.])
y = np.array([0., 0., 0., 0., 0., 0.5, 1., 1., 1., 1., 1., 1., 1.])
print(var(y))
quit()

tck = interpolate.splrep(x, y, s=0)
xnew = np.linspace(0, 12, num=41, endpoint=True)
ynew = interpolate.splev(xnew, tck, der = 0)

plt.figure()
plt.plot(x, y, 'x', xnew, ynew, 'b')
plt.legend(['Linear', 'Cubic Spline'])
plt.show()
