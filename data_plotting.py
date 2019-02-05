from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate
import numpy as np

df1 = pd.read_csv("data_1_ms.csv")
df2 = pd.read_csv("data_23_ms.csv")
df3 = pd.read_csv("data_tr_231.csv")
print(df1.head())

df1.plot.bar(x='1_prot', rot=0, title='1 hour', figsize=(22,10), color = ['red', 'blue'], fontsize=16)
plt.savefig('1hr_cttr')
df2.plot.bar(x='23_prot', rot=0, title='23 hour', figsize=(22,10),color = ['red', 'blue'], fontsize=16)
plt.savefig('23hr_cttr')
df3.plot.bar(x='prot', rot=0, title='1 & 23 hr treated', figsize=(22,10),color = ['red', 'blue'], fontsize=16)
plt.savefig('1_23hr_tr')
plt.show()

