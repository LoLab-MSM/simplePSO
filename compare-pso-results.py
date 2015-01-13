# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 22:38:07 2015

@author: pinojc
"""
import pylab as plt
import numpy as np
data = np.loadtxt("1.txt")
for i in range(2,100):
    tmp = np.loadtxt(str(i)+'.txt')
    data = np.column_stack((data,tmp))
print np.shape(data)
for i in range(105):
    plt.hist(data[i,:],25)
    plt.show()