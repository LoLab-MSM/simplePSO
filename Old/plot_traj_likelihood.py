import pylab as plt
import numpy as np
min_val = []
for i in range(100):
    tmp = np.loadtxt('pso_%s.txt'%str(i))
    plt.semilogy(tmp)
    min_val.append(min(tmp))
plt.show()
plt.hist(min_val)
plt.show()
